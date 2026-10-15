# 003_diploid_baseline — methods and results

Inventory item B1 #9. Preliminary and descriptive. Details in `README`.

## Methods
No new alignment. `scripts/000_baseline.py` reads `../002_alignment/results/03_summary/pair_metrics.tsv` (minimap2 `asm10` + SyRI, 30 within-set pairs) and writes `results/pair_baseline.tsv` and `results/class_summary.tsv`. Per pair it reports:
- SNP/kb and indel/kb over syntenic reference bp
- large SVs (INV + TRANS + DUP rows) per Mb
- the minimum syntenic fraction of the two genomes

Pairs are classed as diploid–diploid, Cbp-subgenome–diploid, or Cbp–Cbp, separately for the Co set (Co39, CoMosher, CoRussian) and the Cr set (Cr145, CrMosher, CgMosher).

`scripts/001_plots.R` (R, ggplot2) draws `results/plots/snp_per_kb.{png,pdf}`: one point per pair, bar = class median, CgMosher pairs marked.

## Results
SNP/kb, median (range):

| set | diploid–diploid | Cbp–diploid | Cbp–Cbp |
|---|---|---|---|
| Co | 0.9 (0.8–1.2) | 5.1 (4.0–14.2) | 8.3 (3.0–9.4) |
| Cr | 8.9 (4.2–13.6) | 9.3 (7.0–18.1) | 9.2 (3.7–9.9) |

- **Co set:** Cbp's Co subgenome is about 4–5× more diverged from the Co diploids than they are from each other. CbpChinese drives the high end (10–14 SNP/kb).
- **Cr set:** the diploid baseline already spans 4–14 SNP/kb because it includes CgMosher (*C. grandiflora*, an outcrosser). Cbp-vs-diploid therefore falls inside the diploid range. Cr145 vs CrMosher (4.2) is the within-species baseline, and the Cg pairs are the more relevant comparison for the Cr progenitor.

## Caveats
n = 3 pairs per class, and pairs share genomes so are not independent. The denominator is syntenic bp only, so divergent regions are under-represented. The Co diploid accessions may be closely related. Assembly size differences (for example CoMosher 181 Mb vs Co39 135 Mb) affect SV counts.
