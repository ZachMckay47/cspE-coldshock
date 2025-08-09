#!/usr/bin/env python3
import json, pathlib, sys, re

CFG = pathlib.Path("config/cspE_config.json")
RES = pathlib.Path("results")
FASTA = None

def read_seq(path: pathlib.Path) -> str:
    txt = path.read_text().splitlines()
    return "".join(l for l in txt if not l.startswith(">")).upper().replace("T","U")

def length_ok(seq, start, end):
    return 1 <= start <= end <= len(seq)

def read_mask(p: pathlib.Path):
    if not p.exists(): return []
    s = p.read_text().strip()
    return [int(x) for x in s.split(",") if x.strip().isdigit()]

def read_p_unpaired(lunp: pathlib.Path, end_pos_1b: int, ulen: int = 6):
    if not lunp.exists(): return None
    with lunp.open() as fh:
        for line in fh:
            if not line.strip(): continue
            cols = line.split()
            if cols[0] == str(end_pos_1b):
                return float(cols[ulen])  # pos P1..P6 -> take P6
    return None

def region(name, s, e):
    return {"name":name, "start":s, "end":e}

def overlaps(a, b):
    return not (a["end"] < b["start"] or b["end"] < a["start"])

def main():
    if not CFG.exists():
        print("[ERROR] config/cspE_config.json not found.")
        sys.exit(1)
    cfg = json.loads(CFG.read_text())
    global FASTA
    FASTA = pathlib.Path(cfg["fasta"])
    if not FASTA.exists():
        print(f"[ERROR] FASTA not found: {FASTA}")
        sys.exit(1)

    seq = read_seq(FASTA)
    sd_start, sd_end, atg = cfg["sd_start"], cfg["sd_end"], cfg["atg_pos"]

    # 1) Bounds checks
    ok1 = length_ok(seq, sd_start, sd_end)
    ok2 = 1 <= atg <= len(seq)
    if not (ok1 and ok2):
        print("[FAIL] SD/ATG out of bounds for sequence length", len(seq))
        sys.exit(1)

    # 2) SD->AUG spacing (want ~7-10 nt)
    spacing = atg - sd_end
    print(f"[INFO] SD window: {sd_start}-{sd_end}; ATG at {atg}; spacing={spacing} nt")
    if not (7 <= spacing <= 10):
        print("[WARN] SD–AUG spacing is outside 7–10 nt. (This can still work, but note it.)")

    # 3) Required result files
    req = [
        RES/"sd_partners.tsv",
        RES/"sd_mutation_mask.txt",
        RES/"sd_summary.md",
        RES/"cspE_transcript_fold_15C.dot",
        RES/"cspE_transcript_fold_lunp"
    ]
    missing = [str(p) for p in req if not p.exists()]
    if missing:
        print("[FAIL] Missing results files:\n  " + "\n  ".join(missing))
        sys.exit(1)
    else:
        print("[OK] All baseline result files present.")

    # 4) SD mask sanity (no protected overlap)
    mask = read_mask(RES/"sd_mutation_mask.txt")
    print(f"[INFO] partner_pos (mutate here): {mask}")
    protected = cfg["protected_regions"]
    bad = []
    for m in mask:
        for pr in protected:
            if pr["start"] <= m <= pr["end"]:
                bad.append((m, pr["name"]))
    if bad:
        print("[FAIL] Some partner_pos collide with protected regions:", bad)
        sys.exit(1)
    else:
        print("[OK] partner_pos do not overlap protected regions.")

    # 5) p_unpaired_SD baseline @15C
    pu = read_p_unpaired(RES/"cspE_transcript_fold_lunp", sd_end, 6)
    if pu is None:
        print("[FAIL] Could not read p_unpaired_SD from _lunp.")
        sys.exit(1)
    print(f"[OK] p_unpaired_SD (u=6) @15°C = {pu:.4f}")

    # 6) Nuisance checks (no extra AUG between SD and ATG)
    between = seq[sd_end:atg-1] if atg-1 > sd_end else ""
    if "AUG" in between:
        print("[WARN] Found 'AUG' between SD and ATG; may create upstream start.")
    else:
        print("[OK] No upstream AUG between SD and ATG.")

    # 7) Emit summary
    print("\n=== Preflight summary ===")
    print(f"Length(seq) = {len(seq)}")
    print(f"SD = {sd_start}-{sd_end}  (len={sd_end-sd_start+1})")
    print(f"ATG pos = {atg}")
    print(f"Spacing = {spacing}")
    print(f"p_unpaired_SD@15C = {pu:.4f}")
    print("Protected:", ", ".join([f"{r['name']}({r['start']}-{r['end']})" for r in protected]))
    print("=========================")

if __name__ == "__main__":
    main()
