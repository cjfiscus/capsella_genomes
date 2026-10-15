# 002_alignment — methods and results

Inventory item B1 #2. Preliminary. Full detail, rules and caveats are in `README`.

## Methods
- **Genomes:** 9 assemblies (Cbp2-2, Co39, Cr145, CbpChinese, CbpRussian, CgMosher, CoMosher, CrMosher, CoRussian). Headers are trimmed to the first token and only the `n_chrom` longest sequences are kept (`000_prep_genomes.sh`).
- **Subgenomes:** each Cbp assembly is split into Co and Cr subgenomes.
  - Cbp2-2 uses the existing subgenome FASTAs (given, not inferred).
  - CbpChinese uses the authors' A/B sets, assigned A = Cr and B = Co on modest evidence.
  - CbpRussian is assigned by a minimap2 chromosome vote against Cbp2-2.
  - Orientation and homolog number go in `results/01_norm/<name>.chrmap.tsv`, and `-` strand chromosomes are reverse-complemented.
  - Source FASTAs are never modified and sequence names are kept (`001_chrom_normalize.sh`).
- **Pairs:** 30 pairs, all 15 within each of the Co set and the Cr/Cg set, with the diploid reference first.
- **Alignment and SV calling:** `minimap2 -ax asm10 --eqx` → `syri -F B` (SyRI 1.8.2) → plotsr per pair (`002_align_syri.sh`). asm10 was chosen over asm20 in pilots (more aligned, faster). TRANS/DUP counts are preset-dependent.
- **Summaries:** `003_summarize.sh` / `003a_summarize.py` give `sv_summary.tsv`, `pair_metrics.tsv`, `sv_size_filtered.tsv` and chained plotsr figures. `003b_inversion_attribution.py` clusters ≥ 1 Mb inversions per genome.
- **TE overlap:** merged TE intervals from the `001_edta` panEDTA annotation are intersected with SyRI regions (`004a_prep_te.py`, `004_te_overlap.py`). Only 16 of 30 pairs qualify, because CbpChinese and CoRussian were dropped from 001.

## Results
- All 30 pairs completed with no errors, and all align sensibly (syntenic fraction 0.46–0.99). The lowest are the Cbp-vs-CgMosher pairs, because CgMosher is a larger assembly.
- **Diploid pairs are closest:** Co39/CoMosher/CoRussian 0.8–1.2 SNP/kb, Cr145 vs CrMosher 4.2. Cbp2-2 and CbpRussian are similar to each other (Co 3.0, Cr 3.7 SNP/kb).
- **CbpChinese is the most diverged Cbp** in both subgenomes (Co 8–10 SNP/kb against the other Cbps, Cr about 9–10 against the other Cbps and 12–18 against diploids). Its subgenomes are also larger (Co 167 Mb, Cr 164 Mb vs 125–141 Mb for the others). A distinct lineage and an assembly or partitioning effect are both possible, so this is a lead only.
- **Large inversions (≥ 1 Mb)** are concentrated in CbpChinese and CgMosher. Cbp2-2 and CbpRussian have almost none. Whether these are real or scaffolding artifacts is not checked.
- **TE content of regions** (16 pairs): syntenic sequence is 9–13% TE in every pair, unaligned sequence is 20–80% TE, and highly diverged regions are 13–49% TE. Most HDR is therefore not annotated TE. The Mosher assemblies have the highest TE fractions in unaligned regions.

## Caveats
28 of the 30 per-pair plotsr figures have not been inspected. The `INV` rows do not reproduce an older "Cr SCF_7 11–19 Mb inversion" note. SyRI small-variant counts are not comparable to the dnadiff SV counts in the manuscript's Fig 1.
