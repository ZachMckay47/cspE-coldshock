#!/usr/bin/env python3
"""
sd_analysis.py

What it does (in plain English):
1) Runs RNAfold at a chosen temperature on your transcript (+1…ATG+N).
2) Reads the dot-bracket structure and tells you which bases are pairing
   with your Shine–Dalgarno (SD) sequence. Those opposite-arm bases are the
   ones to mutate to "open" the SD.
3) Optionally runs RNAplfold and reports p_unpaired for the SD (at length=6),
   which is our quantitative gate for "how open is the SD" at this temperature.
4) Writes tidy outputs in results/: a TSV table, a short Markdown report,
   and a simple mutation mask listing the partner positions to target.

Why this matters:
- If the SD is mostly paired (low p_unpaired), ribosomes struggle to bind,
  especially at cold temperatures. Opening the SD boosts translation.
"""

import argparse
import pathlib
import re
import subprocess
import sys
from typing import List, Tuple, Optional


# ---------------- helpers ----------------

def run_cmd(cmd: List[str], stdin_bytes: Optional[bytes] = None) -> bytes:
    """Run a shell command and return stdout or exit on error."""
    res = subprocess.run(
        cmd,
        input=stdin_bytes,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if res.returncode != 0:
        sys.stderr.write(res.stderr.decode() or f"Command failed: {' '.join(cmd)}\n")
        sys.exit(res.returncode)
    return res.stdout


def read_single_fasta_seq(path: pathlib.Path) -> str:
    """Return the single sequence (no header) from a FASTA file, T→U."""
    txt = path.read_text().splitlines()
    seq = "".join(l.strip() for l in txt if not l.startswith(">")).upper().replace("T", "U")
    if not seq:
        sys.exit(f"[ERROR] No sequence found in {path}")
    return seq


def run_rnafold(fasta_path: pathlib.Path, temp_c: float, out_dot: pathlib.Path) -> Tuple[str, str]:
    """Run RNAfold and return (sequence, structure)."""
    data = fasta_path.read_bytes()
    out = run_cmd(["RNAfold", "--temp", str(temp_c), "--noPS"], stdin_bytes=data)
    out_dot.write_bytes(out)
    lines = [l.strip() for l in out.decode().splitlines() if l.strip()]
    seq = next(l for l in lines if set(l) <= set("ACGTU"))
    struct = next(l for l in lines if any(c in l for c in "() .") or "(" in l or ")" in l).split()[0]
    return seq, struct


def partners_from_dot(struct: str) -> List[int]:
    """Return partner indices (0-based), or -1 if unpaired."""
    stack: List[int] = []
    partner = [-1] * len(struct)
    for i, ch in enumerate(struct):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            j = stack.pop()
            partner[i] = j
            partner[j] = i
    return partner


def run_rnaplfold(
    fasta_path: pathlib.Path,
    temp_c: float,
    winsize: int,
    span: int,
    ulen: int,
    outprefix: str,
) -> Optional[pathlib.Path]:
    """
    Run RNAplfold; whatever <something>_lunp it writes (depends on FASTA header),
    move/rename it to <outprefix>_lunp and return that path.
    """
    # Remember existing *_lunp files so we can detect the new one
    before = set(p.name for p in pathlib.Path(".").glob("*_lunp"))

    data = fasta_path.read_bytes()
    # --auto-id uses the first sequence header as an ID; we force a predictable prefix
    id_prefix = pathlib.Path(outprefix).stem or "seq"
    run_cmd(
        [
            "RNAplfold",
            "--temp", str(temp_c),
            "--winsize", str(winsize),
            "--span", str(span),
            "-u", str(ulen),
            "--auto-id",
            "--id-prefix", id_prefix,
        ],
        stdin_bytes=data,
    )

    # Find the new *_lunp file (the one that wasn't there before)
    after = {p.name for p in pathlib.Path(".").glob("*_lunp")}
    new_files = sorted(list(after - before))
    if not new_files:
        # fallback: if nothing changed, but there is exactly one *_lunp, use it
        existing = sorted(list(after))
        if len(existing) == 1:
            new_files = existing

    if new_files:
        src = pathlib.Path(new_files[-1])  # newest
        dest = pathlib.Path(f"{outprefix}_lunp")
        src.replace(dest)
        return dest

    return None


def read_pu_for_window_end(lunp_path: pathlib.Path, end_pos_1b: int, ulen: int) -> Optional[float]:
    """From <lunp> table, return P<ulen> at <end_pos_1b> (1-based)."""
    with lunp_path.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            cols = line.split()
            if cols[0] == str(end_pos_1b):
                # columns: pos P1 P2 ... P<ulen> ...
                return float(cols[ulen])
    return None


# ---- optional: simple SD autodetect (Purine-rich 6mer near ATG) ----
PWM = {  # toy RBS PWM for 6-nt SD (favor A/G)
    "A": [0.25, 0.20, 0.45, 0.20, 0.45, 0.20],
    "G": [0.45, 0.45, 0.20, 0.45, 0.20, 0.45],
    "C": [0.15, 0.15, 0.15, 0.15, 0.15, 0.15],
    "U": [0.15, 0.20, 0.20, 0.20, 0.20, 0.20],
}


def score_6mer(s: str) -> float:
    return sum(PWM.get(b, [0] * 6)[i] for i, b in enumerate(s))


def autodetect_sd(transcript: str, atg_pos_1b: int) -> Tuple[int, int, str]:
    """Search 6-nt window in [ATG-20 .. ATG-6] with highest PWM score."""
    start = max(1, atg_pos_1b - 20)
    end = max(1, atg_pos_1b - 6)  # inclusive end for SD start
    best = (0.0, -1, "NNNNNN")
    for s in range(start, end + 1):
        e = s + 5
        if e > len(transcript):
            break
        six = transcript[s - 1 : e].upper().replace("T", "U")
        if re.search(r"[^ACGU]", six):
            continue
        sc = score_6mer(six)
        if sc > best[0]:
            best = (sc, s, six)
    if best[1] == -1:
        raise SystemExit("[ERROR] Could not auto-detect SD—provide --sd-start/--sd-end.")
    return best[1], best[1] + 5, best[2]


# ---------------- CLI ----------------

def main():
    p = argparse.ArgumentParser(description="Map SD pairing partners and report SD accessibility.")
    p.add_argument("--fasta", required=True, help="Transcript FASTA (+1..ATG+N)")
    p.add_argument("--temp", type=float, default=15.0, help="Temperature (°C) for folding")
    p.add_argument("--sd-start", type=int, help="SD start (1-based, transcript coords)")
    p.add_argument("--sd-end", type=int, help="SD end   (1-based, transcript coords)")
    p.add_argument("--atg-pos", type=int, help="ATG position (1-based) if using autodetect")
    p.add_argument("--winsize", type=int, default=120)
    p.add_argument("--span", type=int, default=80)
    p.add_argument("--ulen", type=int, default=6)
    p.add_argument("--outdir", default="results")
    p.add_argument("--run-plfold", action="store_true")
    args = p.parse_args()

    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    fasta = pathlib.Path(args.fasta)
    seq = read_single_fasta_seq(fasta)  # U instead of T

    # Determine SD window
    if args.sd_start and args.sd_end:
        sd_start, sd_end = args.sd_start, args.sd_end
        sd_seq = seq[sd_start - 1 : sd_end]
    else:
        if not args.atg_pos:
            sys.exit("[ERROR] Provide --sd-start/--sd-end OR --atg-pos for autodetect.")
        sd_start, sd_end, sd_seq = autodetect_sd(seq, args.atg_pos)

    # Basic sanity checks
    if not (1 <= sd_start <= sd_end <= len(seq)):
        sys.exit(f"[ERROR] SD window {sd_start}-{sd_end} is out of bounds for sequence length {len(seq)}.")

    # 1) RNAfold
    dot_path = outdir / (f"{fasta.stem}_{int(args.temp)}C.dot")
    seq_fold, struct = run_rnafold(fasta, args.temp, dot_path)
    partners = partners_from_dot(struct)

    # 2) partner table for SD
    rows = []
    for pos1 in range(sd_start, sd_end + 1):
        i = pos1 - 1
        j = partners[i]
        rows.append(
            {
                "pos": pos1,
                "base": seq_fold[i],
                "partner_pos": (j + 1) if j >= 0 else "-",
                "partner_base": seq_fold[j] if j >= 0 else "-",
                "pair": f"{seq_fold[i]}{seq_fold[j]}" if j >= 0 else "-",
            }
        )

    # 3) RNAplfold (optional)
    pu = None
    lunp_path = None
    if args.run_plfold:
        lunp_path = run_rnaplfold(fasta, args.temp, args.winsize, args.span, args.ulen, str(outdir / fasta.stem))
        if lunp_path:
            pu = read_pu_for_window_end(lunp_path, sd_end, args.ulen)

    # 4) Write outputs
    tsv = outdir / "sd_partners.tsv"
    with tsv.open("w") as fh:
        fh.write("pos\tbase\tpartner_pos\tpartner_base\tpair\n")
        for r in rows:
            fh.write(f"{r['pos']}\t{r['base']}\t{r['partner_pos']}\t{r['partner_base']}\t{r['pair']}\n")

    # 4b) also write a simple mask of positions-to-mutate (opposite arm)
    mask_path = outdir / "sd_mutation_mask.txt"
    with mask_path.open("w") as fh:
        muts = [str(r["partner_pos"]) for r in rows if r["partner_pos"] != "-"]
        fh.write(",".join(muts) + "\n")

    md = outdir / "sd_summary.md"
    with md.open("w") as fh:
        fh.write("# SD Accessibility Summary\n\n")
        fh.write(f"- FASTA: `{fasta.name}`\n")
        fh.write(f"- Temperature: **{args.temp} °C**\n")
        fh.write(f"- SD (1-based): **{sd_start}–{sd_end}**  | SD seq: `{sd_seq}`\n")
        if args.atg_pos:
            fh.write(f"- ATG position (1-based): **{args.atg_pos}**\n")
        if pu is not None:
            fh.write(f"- p_unpaired_SD (u={args.ulen}): **{pu:.4f}**\n")
            fh.write(f"- Source: `{lunp_path.name}`\n")
        fh.write(f"- Mutation mask (partner_pos list): `{mask_path.name}`\n")
        fh.write("\n**Guidance:** mutate `partner_pos` bases to open the SD; do **not** change the SD itself, AUG, or +1..+5.\n")
        fh.write("\n## SD pairing partners (opposite arm = mutate here)\n\n")
        fh.write("| pos | base | partner_pos | partner_base | pair |\n|---:|:---:|---:|:---:|:---:|\n")
        for r in rows:
            fh.write(f"| {r['pos']} | {r['base']} | {r['partner_pos']} | {r['partner_base']} | {r['pair']} |\n")

    print("Wrote:")
    print(" ", tsv)
    print(" ", mask_path)
    print(" ", md)
    print(" ", dot_path)
    if pu is not None and lunp_path is not None:
        print(f" p_unpaired_SD = {pu:.4f}")
        print(" ", lunp_path)


if __name__ == "__main__":
    main()
