# cspE Cold-Shock: Opening the Shine–Dalgarno (SD) with RNAfold/RNAplfold

**One-line goal:** At **15 °C**, the mRNA of *E. coli* **cspE** tends to fold so the **Shine–Dalgarno (SD)** is hidden. We map which bases hide the SD and design minimal mutations (not touching the SD/AUG) to **open** it for better translation.

---

## Why this matters (plain English)

- A **gene** is like a recipe. The cell reads the recipe using **ribosomes**.
- To start reading, the ribosome needs to grab a small “handle” on the mRNA called the **Shine–Dalgarno (SD)** sequence (a 6-letter stretch ~6–9 bases before the start codon **AUG**).
- At **cold temperatures**, mRNAs form tighter **hairpins** that can **hide** the SD.
- If the SD is hidden, the ribosome can’t start easily → **low protein**.
- We use **RNAfold**/**RNAplfold** (ViennaRNA tools) to:
  1) find which bases are **pairing with** (hiding) the SD,
  2) measure how often the SD is **unpaired** (accessible).  
     > The number we track is **p_unpaired_SD** (probability the 6-nt SD is unpaired).  
     > Rules of thumb: **< 0.25 = hidden**, **≥ 0.30 = good**.

---

## What this repo gives you

1. A **reproducible pipeline** to:
   - run **RNAfold** at 15 °C,
   - list the **exact bases** that pair with the SD (the ones you’ll mutate),
   - run **RNAplfold** and record **p_unpaired_SD**.
2. Clean outputs in `results/` you can drop in slides or a report.

---

## Folder layout
IDEC/
├─ data/ # input sequences you provide
│ └─ cspE_transcript_fold.fasta # +1 … ATG+30 (176 nt for cspE)
├─ scripts/
│ ├─ sd_analysis.py # main analysis (fold + SD partners + p_unpaired_SD)
│ └─ preflight_checks.py # sanity checks before mutagenesis
├─ results/ # auto-generated outputs
│ ├─ sd_partners.tsv # SD base ↔ opposite-arm base (mutate these)
│ ├─ sd_mutation_mask.txt # comma list of positions to mutate
│ ├─ sd_summary.md # short report incl. p_unpaired_SD
│ ├─ cspE_transcript_fold_15C.dot # RNAfold dot-bracket used
│ └─ cspE_transcript_fold_lunp # RNAplfold unpaired probs (u=1..6…)
├─ config/
│ └─ cspE_config.json # SD window, ATG, protected regions, temps
├─ run.sh # one command to (re)run the analysis
└─ README.md


---

## Quick start (WSL Ubuntu)

1) **Install ViennaRNA** (if not already):
```bash
sudo apt update
sudo apt install -y vienna-rna      # provides RNAfold, RNAplfold

Put your transcript FASTA in data/ as cspE_transcript_fold.fasta.
For cspE we use +1 … ATG+30 (the whole 5′-UTR + 30 nt of CDS).

Run the analysis:

bash
Copy
Edit
./run.sh
You’ll get files in results/ and the terminal will print:

scss
Copy
Edit
p_unpaired_SD (u=6) at 15°C: 0.03...
What the outputs mean
results/sd_partners.tsv
Table showing each SD base and the partner_pos (the base it pairs with in the hairpin).
→ Mutate the partner_pos bases (opposite arm) to open the SD.
→ Do NOT mutate the SD itself, AUG, or +1..+5.

results/sd_mutation_mask.txt
Comma-separated list of the positions to mutate (opposite arm).

results/cspE_transcript_fold_lunp
RNAplfold table with unpaired probabilities. The P6 value at the SD end position is p_unpaired_SD.

results/sd_summary.md
Human-readable mini-report with the SD window, sequence, and (if run) p_unpaired_SD.

Preflight checks (so you don’t waste time)
Run:

bash
Copy
Edit
python3 scripts/preflight_checks.py
You should see:

SD window, ATG position, and spacing (7–10 nt is common),

All required files present,

partner_pos does not collide with protected regions,

A baseline p_unpaired_SD at 15 °C (example: 0.0316).

Optional multi-temp baseline:

Check 15/20/37 °C so you don’t break warm-temp behavior (script available on request).

What to mutate (and what not)
Mutate: the opposite-arm bases listed in sd_mutation_mask.txt (these are what pair with the SD).
Do not mutate:

SD (6 nt) itself,

AUG (start codon),

+1..+5 of the transcript,

promoter DNA (−35/−10/spacer) — not part of the transcript we fold.

Strategy: start with 1–2 single-nt changes on the opposite arm that break Watson–Crick pairing to the SD base (e.g., break GC→GA/AC, AU→AC/CC), re-run ./run.sh, and watch p_unpaired_SD. Stop when ≥ 0.30 at 15 °C.

Glossary (super short)
SD (Shine–Dalgarno): the “handle” ribosomes grab; ~6 bases long; purine-rich.

AUG: start codon where translation begins.

5′-UTR: the front part of mRNA before AUG; can fold into hairpins.

RNAfold: predicts one most-stable structure (shape).

RNAplfold: estimates probabilities across many shapes; here we read P6 at the SD end = p_unpaired_SD.

Troubleshooting
“RNAplfold not found” → sudo apt install -y vienna-rna

No _lunp file in results → run ./run.sh again; script now auto-captures whatever _lunp name ViennaRNA uses.

Wrong SD window → set it in config/cspE_config.json and re-run.

Working in the wrong (Windows) terminal → use the Ubuntu (WSL) terminal and ensure the VS Code window says WSL: Ubuntu (green badge).

