#!/usr/bin/env python3
"""
mutate_sd.py

Goal
----
Given a transcript (+1..ATG+N) and the list of opposite-arm positions that
pair to the SD (from sd_analysis), propose minimal mutations (single/double)
that open the SD at cold temperature(s). For each mutant we compute
p_unpaired_SD (u=6) with RNAplfold and rank mutants.

Inputs
------
- config/cspE_config.json      : SD window, ATG, temps, protected regions
- data/cspE_transcript_fold.fasta
- results/sd_partners.tsv      : SD partner mapping (for base-aware proposals)
- results/sd_mutation_mask.txt : 1-based positions to mutate (opposite arm)

Outputs (results/mutants/)
--------------------------
- mutants.csv   : id, n_mut, muts, p_unpaired at each temp, passes, notes
- topN.fasta    : top mutants (RNA + DNA) for quick copy/paste
- a folder of *temporary* plfold files is not kept (we only store scores)

Usage
-----
python3 scripts/mutate_sd.py \
  --config config/cspE_config.json \
  --temps 15 20 37 \
  --max-singles 100 \
  --max-doubles 200 \
  --top 25

Notes
-----
- We only mutate positions listed in sd_mutation_mask.txt.
- We never touch SD, AUG, or +1..+5 (verified again as safety).
- Base proposals are chosen to break pairing to the SD base.
"""

import argparse, pathlib, subprocess, sys, json, re, itertools, time
from typing import List, Dict, Tuple, Optional

# --------------------- utilities ---------------------

def run_cmd(cmd: List[str], stdin_bytes: Optional[bytes] = None) -> bytes:
    res = subprocess.run(cmd, input=stdin_bytes,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if res.returncode != 0:
        sys.stderr.write(res.stderr.decode() or f"Command failed: {' '.join(cmd)}\n")
        sys.exit(res.returncode)
    return res.stdout

def read_fasta_seq(path: pathlib.Path) -> str:
    txt = path.read_text().splitlines()
    seq = "".join(l for l in txt if not l.startswith(">")).strip()
    if not seq:
        sys.exit(f"[ERROR] No sequence in {path}")
    return seq  # DNA letters allowed; we convert to RNA for folding

def dna_to_rna(s: str) -> str:
    return s.upper().replace('T','U')

def rna_to_dna(s: str) -> str:
    return s.upper().replace('U','T')

def read_mask(path: pathlib.Path) -> List[int]:
    if not path.exists():
        sys.exit(f"[ERROR] mutation mask not found: {path}")
    s = path.read_text().strip()
    return [int(x) for x in s.split(",") if x.strip().isdigit()]

def parse_partners_tsv(path: pathlib.Path) -> Dict[int, Dict[str,str]]:
    """
    Return: map partner_pos -> { 'sd_base':..., 'partner_base':..., 'sd_pos':... }
    Only rows where partner_pos != '-' are used.
    """
    m = {}
    if not path.exists():
        sys.exit(f"[ERROR] sd_partners.tsv not found: {path}")
    for line in path.read_text().splitlines():
        if not line or line.startswith("pos"):  # header
            continue
        cols = line.strip().split("\t")
        sd_pos = int(cols[0]); sd_base = cols[1]
        partner_pos = cols[2]
        partner_base = cols[3]
        if partner_pos != '-':
            p = int(partner_pos)
            m[p] = {'sd_pos': sd_pos, 'sd_base': sd_base, 'partner_base': partner_base}
    return m

# ViennaRNA: run RNAplfold and return path to the *_lunp file
def run_rnaplfold_for_seq(rna_seq: str, temp_c: float, winsize: int, span: int, ulen: int, jobtag: str) -> pathlib.Path:
    before = set(p.name for p in pathlib.Path('.').glob('*_lunp'))
    data = f">{jobtag}\n{rna_seq}\n".encode()
    run_cmd(['RNAplfold','--temp',str(temp_c),'--winsize',str(winsize),'--span',str(span),'-u',str(ulen)], stdin_bytes=data)
    after = set(p.name for p in pathlib.Path('.').glob('*_lunp'))
    new = sorted(list(after - before))
    if not new:
        exist = sorted(list(after))
        if len(exist) == 1:
            new = exist
    if not new:
        sys.exit("[ERROR] RNAplfold didn’t produce a *_lunp file")
    return pathlib.Path(new[-1])

def read_p_unpaired(lunp: pathlib.Path, end_pos_1b: int, ulen: int) -> Optional[float]:
    with lunp.open() as fh:
        for line in fh:
            if not line.strip(): continue
            cols = line.split()
            if cols[0] == str(end_pos_1b):
                return float(cols[ulen])
    return None

# pairing rules helper: bases to A/U/G/C that AVOID strong pairing to sd_base
def anti_pairing_choices(sd_base: str, current_partner: str) -> List[str]:
    sd_base = sd_base.upper().replace('T','U')
    current_partner = current_partner.upper().replace('T','U')
    bases = ['A','U','G','C']
    avoid = set()
    # strong + wobble
    if sd_base == 'A': avoid = {'U'}
    elif sd_base == 'U': avoid = {'A'}
    elif sd_base == 'G': avoid = {'C','U'}
    elif sd_base == 'C': avoid = {'G'}
    allowed = [b for b in bases if b not in avoid]
    # don't propose same base (no-op)
    return [b for b in allowed if b != current_partner]

def has_upstream_aug(rna_seq: str, sd_end: int, atg_pos: int) -> bool:
    if atg_pos-1 <= sd_end: return False
    window = rna_seq[sd_end:atg_pos-1]
    return 'AUG' in window

def gc_content(s: str) -> float:
    s = s.upper()
    if not s: return 0.0
    return (s.count('G') + s.count('C')) / len(s)

# ------------------ main mutagenesis ------------------

def build_mutants(ref_rna: str,
                  mask_positions: List[int],
                  partner_map: Dict[int,Dict[str,str]],
                  max_singles: int,
                  max_doubles: int) -> List[Dict]:
    """
    Produce single and double mutants using anti-pairing choices at each
    partner position. Returns a list of dicts with fields:
      {'id','mutations':[(pos,ref,alt),...],'seq_rna'}
    """
    muts = []

    # singles
    for p in mask_positions:
        ref = ref_rna[p-1]
        sd_base = partner_map.get(p, {}).get('sd_base', None)
        partner_base = partner_map.get(p, {}).get('partner_base', ref)
        if sd_base is None:
            # fallback: if not in map, propose all ≠ ref
            choices = [b for b in "AUGC" if b != ref]
        else:
            choices = anti_pairing_choices(sd_base, partner_base)
        for alt in choices:
            s_list = list(ref_rna)
            s_list[p-1] = alt
            muts.append({
                'id': f"m{p}{ref}>{alt}",
                'mutations': [(p, ref, alt)],
                'seq_rna': "".join(s_list)
            })
    # cap
    if max_singles and len(muts) > max_singles:
        muts = muts[:max_singles]

    # best alt per position (by heuristic: prefer A over G over C over U to open)
    # We'll use this to make doubles non-explosive.
    best_alt_by_pos = {}
    for m in muts:
        pos, ref, alt = m['mutations'][0]
        # simple priority list to bias away from GC
        prio = {'A':3,'G':2,'C':1,'U':0}
        key = (pos, ref)
        existing = best_alt_by_pos.get(key)
        if existing is None or prio[alt] > prio[existing]:
            best_alt_by_pos[key] = alt

    # doubles
    doubles = []
    positions = [p for p in mask_positions]
    for p1, p2 in itertools.combinations(positions, 2):
        ref1 = ref_rna[p1-1]; ref2 = ref_rna[p2-1]
        alt1 = best_alt_by_pos.get((p1, ref1))
        alt2 = best_alt_by_pos.get((p2, ref2))
        if not alt1 or not alt2: continue
        s_list = list(ref_rna)
        s_list[p1-1] = alt1
        s_list[p2-1] = alt2
        doubles.append({
            'id': f"m{p1}{ref1}>{alt1}_m{p2}{ref2}>{alt2}",
            'mutations': [(p1,ref1,alt1),(p2,ref2,alt2)],
            'seq_rna': "".join(s_list)
        })
    if max_doubles and len(doubles) > max_doubles:
        doubles = doubles[:max_doubles]

    return muts + doubles

def main():
    ap = argparse.ArgumentParser(description="Mutagenesis around SD: propose single/double mutants and rank by p_unpaired_SD.")
    ap.add_argument('--config', required=True, help='JSON config (cspE_config.json)')
    ap.add_argument('--winsize', type=int, default=120)
    ap.add_argument('--span', type=int, default=80)
    ap.add_argument('--ulen', type=int, default=6)
    ap.add_argument('--temps', type=float, nargs='+', default=[15.0], help='Temperatures to evaluate (°C)')
    ap.add_argument('--max-singles', type=int, default=100)
    ap.add_argument('--max-doubles', type=int, default=200)
    ap.add_argument('--top', type=int, default=25)
    ap.add_argument('--outdir', default='results/mutants')
    args = ap.parse_args()

    cfg = json.loads(pathlib.Path(args.config).read_text())
    fasta = pathlib.Path(cfg['fasta'])
    sd_start, sd_end = int(cfg['sd_start']), int(cfg['sd_end'])
    atg_pos = int(cfg['atg_pos'])
    protected = cfg.get('protected_regions', [])
    outdir = pathlib.Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    seq_dna = read_fasta_seq(fasta)
    ref_rna  = dna_to_rna(seq_dna)
    mask_positions = read_mask(pathlib.Path('results/sd_mutation_mask.txt'))
    partner_map    = parse_partners_tsv(pathlib.Path('results/sd_partners.tsv'))

    # Safety: none of the mask positions may touch protected regions
    for m in mask_positions:
        for pr in protected:
            if pr['start'] <= m <= pr['end']:
                sys.exit(f"[ERROR] Proposed mutation position {m} overlaps protected region {pr['name']}")

    mutants = build_mutants(ref_rna, mask_positions, partner_map, args.max_singles, args.max_doubles)

    results_rows = []
    t0 = time.time()

    for k, mut in enumerate(mutants, 1):
        scores = {}
        passes = True
        notes = []

        # basic filters: spacing window shouldn't gain AUG (very unlikely since we don't mutate there)
        if has_upstream_aug(mut['seq_rna'], sd_end, atg_pos):
            passes = False
            notes.append("new_upstream_AUG")

        # optional: GC change of whole 5'UTR not more than +/-5%
        utr_ref = ref_rna[:atg_pos-1]
        utr_mut = mut['seq_rna'][:atg_pos-1]
        if abs(gc_content(utr_mut) - gc_content(utr_ref)) > 0.05:
            notes.append("ΔGC>5%")

        for T in args.temps:
            lunp = run_rnaplfold_for_seq(mut['seq_rna'], T, args.winsize, args.span, args.ulen, jobtag=f"mut{k}_{int(T)}C")
            pu = read_p_unpaired(lunp, sd_end, args.ulen)
            if pu is None:
                passes = False
                notes.append(f"no_pu@{T}")
                pu = float('nan')
            scores[T] = pu

        row = {
            'id': mut['id'],
            'n_mut': len(mut['mutations']),
            'muts': ";".join([f"{p}{r}>{a}" for (p,r,a) in mut['mutations']]),
            **{f"p_unpaired_{int(T)}C": scores[T] for T in args.temps},
            'passes': "yes" if passes else "no",
            'notes': ",".join(notes),
            'seq_rna': mut['seq_rna'],
            'seq_dna': rna_to_dna(mut['seq_rna'])
        }
        results_rows.append(row)

    # rank by primary temperature (first in list)
    primary = args.temps[0]
    results_rows.sort(key=lambda r: (-(r.get(f"p_unpaired_{int(primary)}C") or -1), r['n_mut']))

    # write CSV
    import csv
    csv_path = outdir / "mutants.csv"
    with csv_path.open("w", newline="") as fh:
        fieldnames = ['id','n_mut','muts'] + [f"p_unpaired_{int(T)}C" for T in args.temps] + ['passes','notes','seq_rna','seq_dna']
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in results_rows: w.writerow(r)

    # write topN FASTA (RNA + DNA)
    topN = results_rows[:args.top]
    fa = outdir / "topN.fasta"
    with fa.open("w") as fh:
        for r in topN:
            fh.write(f">{r['id']}_RNA  n_mut={r['n_mut']}  p15={r.get('p_unpaired_15C','')}\n{r['seq_rna']}\n")
            fh.write(f">{r['id']}_DNA\n{r['seq_dna']}\n")

    # write markdown summary
    md = outdir / "README_mutants.md"
    with md.open("w") as fh:
        fh.write("# Mutagenesis summary\n\n")
        fh.write(f"- total mutants evaluated: **{len(results_rows)}**\n")
        fh.write(f"- ranked by p_unpaired_SD at **{int(primary)} °C**\n")
        fh.write(f"- top written to: `{fa}`\n")
        fh.write(f"- full table: `{csv_path}`\n\n")
        fh.write("## Top candidates\n\n")
        fh.write("| id | n_mut | p_unpaired_15C | notes |\n|---|---:|---:|---|\n")
        for r in topN:
            p15 = r.get("p_unpaired_15C", "")
            fh.write(f"| {r['id']} | {r['n_mut']} | {p15:.4f} | {r['notes']} |\n")
        fh.write("\n")
    dt = time.time() - t0
    print(f"Wrote:\n  {csv_path}\n  {fa}\n  {md}\nDone in {dt:.1f}s.")
if __name__ == "__main__":
    main()
