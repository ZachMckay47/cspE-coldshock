# SD Accessibility Summary

- FASTA: `cspE_transcript_fold.fasta`
- Temperature: **37.0 °C**
- SD (1-based): **134–139**  | SD seq: `AAGGUA`
- p_unpaired_SD (u=6): **0.1611**
- Source: `cspE_transcript_fold_lunp`
- Mutation mask (partner_pos list): `sd_mutation_mask.txt`

**Guidance:** mutate `partner_pos` bases to open the SD; do **not** change the SD itself, AUG, or +1..+5.

## SD pairing partners (opposite arm = mutate here)

| pos | base | partner_pos | partner_base | pair |
|---:|:---:|---:|:---:|:---:|
| 134 | A | - | - | - |
| 135 | A | 68 | U | AU |
| 136 | G | 67 | C | GC |
| 137 | G | 65 | C | GC |
| 138 | U | 64 | A | UA |
| 139 | A | 63 | U | AU |
