# Capsella genomes manuscript — analysis inventory & working list

**Companion to:** `PUBLICATION_PLAN.md` (writing/formatting/consistency blockers — not duplicated here)
**Purpose:** Part A = every analysis currently in the paper. Part B = analyses still needed, comparative
genomics section first.
**Date:** 2026-09-17

---

## Part A — Analyses already present in the manuscript

### A1. Genome assembly & annotation

| Analysis | Tool | Data | Location | Result |
|---|---|---|---|---|
| De novo assembly | Unstated (Weigel lab) | Cbp 721-48 (278.3 Mbp), Cr 879/721-48 (158.6 Mbp), Co 2007-03 (136.2 Mbp) | Results §1, Methods | Substantial contiguity improvement over prior assemblies |
| BUSCO completeness | BUSCO 5.2.2, brassicales_odb10 | All 3 assemblies | Results §1, Table 1/S2 | 99.4% (Cbp), 99% (Cr), ~99% (Co) |
| Whole-genome alignment/synteny | nucmer (MUMmer4) | Subgenome vs. progenitor, pairwise | Methods | Basis for SV calls below |
| Ortholog/homolog identification | GENESPACE 0.9.3 (OrthoFinder 2.5.4) | 3 genomes | Methods | Basis for CNV/GO analyses |
| Comparative assembly stats table | — | New assemblies vs. Kasianov 2017 (Cbp), Slotte 2013 (Cr) | Table 1, Table S2 | Stats-only comparison, no alignment |

*Note: assembly method itself is undocumented — flagged in `PUBLICATION_PLAN.md` B4/§2.1.*

### A2. Subgenome structural/gene-content divergence (Fig 1, 11 panels; Figs S1–S3)

| Analysis | Tool | Data | Location | Result |
|---|---|---|---|---|
| SV calling | dnadiff (MUMmer) | CbpCo vs. Co; CbpCg vs. Cr | Fig 1A–C, S1 | 144,111 total SVs; CbpCg has more SVs than CbpCo even size-normalized |
| SNP calling + mutation spectrum | GATK4 (HaplotypeCaller/GenotypeGVCFs), bwa-mem2, samtools, bcftools, vcftools | Subgenomes vs. progenitors | Fig 1D, Methods | 1,042,229 SNPs (CbpCo), 2,204,354 (CbpCg); comparison group for Fig 1D ambiguous in caption |
| SNP effect annotation + nonsyn:syn ratio | Ensembl VEP v95 | Same as above | Results §3 | 0.97 (CbpCo) vs. 0.83 (CbpCg) — raw ratio, not normalized dN/dS (see B2) |
| Gene CNV analysis | GENESPACE/OrthoFinder orthologs | Subgenomes vs. progenitors | Fig 1E–I | 5,605/20,582 genes show CNV; largest category is loss in both subgenomes (2,741) |
| GO enrichment on CNV sets | topGO 2.46.0, weight01, Fisher exact p<0.05 | CNV-gain and both-subgenome-loss gene sets | Table S3, S4 | 46 terms (gain), 111 terms (both-loss) |
| TE/repeat analysis | EDTA 2.0, custom library (A. lyrata, A. thaliana, Cbp2-2, Cr145, Co39) | Co39, CbpCo, CbpCg, Cr145 | Fig 1J–K, S2, S3 | TE copy number up in CbpCo vs. Co, down in CbpCr vs. Cr; >14 Mbp LTR/Copia loss in CbpCr |

### A3. Population structure & phenotype (Figs 2–3, S4–S9)

| Analysis | Tool | Data | Location | Result |
|---|---|---|---|---|
| RAD-seq genotyping | Stacks 2.60, bwa-mem2 | 1,441→1,251 samples | Methods | 30,668 SNPs |
| ADMIXTURE | ADMIXTURE 1.3.0, Evanno ΔK, pophelper | K=1–20 (Methods) / K=1–10 (Fig S4, inconsistent) | Fig 2A, S4, S5 | K=2 best (Cluster1 n=671, Cluster2 n=559, 21 admixed); K=3 map in Fig S5 |
| PCA of allele frequencies | adegenet 2.1.5 | All samples | Fig 2B | PC1 16%, PC2 7%, clean Cluster1/2 separation |
| Köppen climate classification | kgc 1.0.0.2 (R) | Sample collection sites | Fig 2C | Cluster1 in D/E climates, Cluster2 dry/temperate |
| XP-CLR selection scan | xpclr 1.1.2 | Cluster1 vs. Cluster2, 100kb windows | Fig S6, Table S5 | Top 1% windows enriched for stress/dev GO terms — rank cutoff only, no significance test (see B2) |
| Common garden phenotyping | Wilcoxon rank-sum, Fisher's exact (R) | Osnabrück 2015 garden | Fig 2D–G | Cluster2 earlier flowering, more branches, better germination; leaf-shape freq differs (χ²=30.318, p=4.2e-06) |
| Flow cytometry genome size | Undocumented method | By cluster | Fig 2H | Cluster2 reduced genome size vs. Cluster1 |
| Colonization/dispersal analysis | plink 1.9 Hamming distance | 4 focal regions (West Coast NA, other NA, Japan, GB/Ireland) | Fig 3, S8, S9 | Multi-lineage colonization narrative; distance metric/test unnamed in captions |

### A4. GWA / climate adaptation (Fig 4, S10–S11)

| Analysis | Tool | Data | Location | Result |
|---|---|---|---|---|
| Bioclimatic PCA | 19 WorldClim2 variables | Sample sites | Fig 4B, S10 | PC1 (temp, 33%) + PC2 (precip, 29%) = 62% variance |
| Phenotype-climate correlation | — | Flowering/branches ~ PC1; height ~ PC2 | Fig 4C, S11 | p<0.05; R² shown in-panel only |
| Genotype imputation | STITCH 1.6.6 | RAD SNPs | Methods | K chosen via INFO_SCORE |
| GWA | GEMMA 0.98.5, mixed linear model, centered kinship, Bonferroni α=0.05 | PC1, PC2 | Fig 4D–E | 45 regions (PC1), 90 regions (PC2), 24 shared; no QQ plots/λ reported (see B2) |
| GO enrichment on GWA hits | topGO 2.46.0 | Nearest genes to hits | Fig 4F–G, Table S7, S8 | 16 terms (PC1), 38 terms (PC2); "Enrichment" column is actually −log10(p), mislabeled |
| Candidate gene curation | Manual, vs. Arabidopsis homologs | Top GWA SNPs | Table S9 | REPRESSOR OF VERNALIZATION 1, ABA/osmotic receptor kinase, INCREASED DNA METHYLATION 3, ARF19, etc. |

*Caption/labeling ambiguities in the above (Fig 1D comparison group, Fig 2D–H test identity, Fig 3
metric, ADMIXTURE K range mismatch) are tracked in detail in `PUBLICATION_PLAN.md` §1.2 — not repeated
here.*

---

## Part B — Analyses needed

### B1. Comparative genomics section (NEW — main-text scope)

Comparing the paper's 3 genomes (Cbp 721-48, Cr 879/721-48, Co 2007-03) against 6 additional published
genomes: **Cbp-China**, **Cbp-Russia**, **Cg-Mosher**, **Co-Mosher**, **Cr-Mosher**, **Co-Russia** — 9
genomes total. All 6 external genomes are already downloaded.

This is planned as a **main-text section**, not a supplement — it should support a real headline claim,
most likely: **single vs. multiple independent origins of *C. bursa-pastoris***, and/or **shared vs.
convergent gene-loss/SV/TE patterns across independently formed Cbp lineages**. Having 3 independently
assembled Cbp genomes is the key leverage here — it turns the paper's existing subgenome CNV/SV/TE
results (Part A2) into a testable question about whether those patterns are identity-by-descent or
convergent evolution.

#### Alignment strategy: hybrid (GENESPACE + SyRI)

- **Gene-anchored layer:** extend GENESPACE (already in the pipeline, OrthoFinder backend) to all 9
  genomes. Backbone for synteny and orthology-based CNV/gene-loss comparison across lineages.
- **Base-level layer:** pairwise **SyRI + minimap2** — the natural upgrade from the dnadiff approach
  already used for Fig 1 (SyRI explicitly classifies SVs/indels rather than inferring them from delta
  alignments). Target ~15–20 biologically motivated pairs (each Cbp vs. its 2 donors, Cbp-vs-Cbp,
  diploid-vs-diploid) rather than all 36 pairwise combinations.
- **Not using** a Cactus/pangenome-graph approach as primary — heavier to build/interpret than this
  question requires. Revisit only if pairwise results prove insufficient for simultaneous cross-lineage
  homologous-locus comparison.

#### Checklist of analyses

| # | Analysis | Tool | Question answered | Priority | Reuses existing pipeline? |
|---|---|---|---|---|---|
| 1 | Extend GENESPACE orthology/synteny to 9 genomes | GENESPACE/OrthoFinder | Consistent subgenome structure/orthogroups across 3 independent Cbp assemblies? | High | Yes |
| 2 | Pairwise SyRI+minimap2 across ~15–20 key pairs | SyRI | SV/indel content per Cbp lineage; shared vs. convergent landscape | High | Upgrade of dnadiff step |
| 3 | Gene-loss/CNV set overlap across 3 Cbp genomes (Jaccard + permutation test) | Custom, on GENESPACE output | Single-origin (shared) vs. multi-origin (independent) test for the paper's own CNV result | High | Yes |
| 4 | Re-run EDTA on all 9 genomes with same TE library | EDTA 2.0 | Is subgenome TE gain/loss shared across Cbp lineages or lineage-specific? | High | Yes |
| 5 | Phylogenomics from 9-genome orthologs (4-fold sites or ASTRAL gene tree) | IQ-TREE/RAxML | Monophyletic Cbp (single origin) vs. recurrent origin? | High | No — new |
| 6 | Extend Ks/synonymous-divergence work to all 9 genomes | Existing Ks pipeline | Same subgenome-donor split time across lineages, or different (recurrent polyploidization)? | High | Partial |
| 7 | SV polarization using Co-Mosher/Co-Russia/Cr-Mosher as outgroups | SyRI output from #2 | Fixes reference-dependent ins/del asymmetry (`PUBLICATION_PLAN` §3.2) | High | Yes |
| 8 | Indel size-spectrum comparison at homologous loci, genic vs. intergenic, across 3 Cbp lineages | Derived from #2 + GENESPACE coords | Converts the "maybe indel signal" hunch into a testable claim | Med–High | No — new |
| 9 | Diploid-diploid baseline divergence (Co×3, Cr/Cg×2–3 accessions) | Reuses #2 output | What does normal within-species divergence look like, to judge if Cbp-vs-donor divergence is unusual | Med | Yes |
| 10 | Homeologous exchange (HE) scan replicated across all 3 Cbp genomes | Ancestry-HMM-style window scan | Shared (IBD) vs. independent HE tracts per lineage | Med | No — new |
| 11 | Whole-genome pairwise identity heatmap/dot-plot, 9 genomes | Output of #1/#2 | Summary figure | Low–Med | Yes |

Items 1, 3, 4, 6, 7, 9 are cheap reuse of pipelines already built for this paper. Items 5, 8, 10 are
genuinely new analytical work. Since this is main-text-bound, sequence 5 and 6 early (Phase 2) rather
than gating them behind a pilot — they're the analyses most likely to produce the headline finding.

#### Open items to confirm before/while starting

- Confirm annotation (GFF3) completeness/quality for all 6 external genomes — items 1, 3, 4, 6, 10 need
  real gene models, not just contigs.
- Confirm whether "Cbp — Russian group" is the same assembly already cited as `Cbp_Kasianov2017` in
  Table 1/S2, or a distinct one.
- Confirm whether the Russian Co accession is the same wild lineage the Introduction already treats as
  the maternal donor, or independent Co diversity — changes the interpretation of items 3/5/6.
- If assembly quality varies widely across the 6 external genomes, gene-loss/CNV comparisons need the
  same read-depth-validation control as B2 below, or loss patterns will be confounded with assembly
  quality rather than biology.

#### Suggested order of work

- **Phase 1** — verify annotation completeness for the 6 external genomes; resolve the
  Kasianov/Russian-Cbp identity question and the Russian-Co lineage question.
- **Phase 2** — GENESPACE 9-genome rerun (#1); EDTA rerun (#4); SyRI pairwise (#2); begin phylogenomics
  (#5) and Ks extension (#6) in parallel.
- **Phase 3** — gene-loss/CNV overlap test (#3); SV polarization (#7); diploid baseline (#9).
- **Phase 4** — indel size-spectrum analysis (#8); HE replication (#10) if the single/multi-origin
  result warrants it; summary figure (#11); write-up.

### B2. Other analyses needed (carried forward from `PUBLICATION_PLAN.md`)

Not comparative-genomics, but still open. See `PUBLICATION_PLAN.md` for full detail — listed here only
as checklist items:

- **Read-depth validation of gene loss** (§3.1) — highest scientific priority overall; map Cbp WGS reads
  to donor gene models called "lost." *Overlaps with B1 items 3/9 — do once, applies to both.*
- **Proper dN/dS or πN/πS with CIs** (§3.3) — replace the current raw nonsyn:syn ratio.
- **QQ plots + genomic inflation factor λ** for both GWAs (§3.4) — non-negotiable for GWA.
- **ADMIXTURE cross-validation error curve** (§3.5) — cheap, data already exists.
- **XP-CLR significance threshold via permutation** (§3.6) — currently a rank cutoff only.
- **Organellar genome assembly/analysis** (§2.5) — confirms C. orientalis as maternal donor.
- **Ks/synonymous divergence between subgenomes and donors** (§2.5) — *overlaps with B1 item 6, do once.*
- **Assembly-completeness / genome-size reconciliation** (§2.3) — k-mer spectrum, GenomeScope.
- **Expanded Table 1** with QV (Merqury), LAI (from existing EDTA output), contig N50/L50,
  chromosome-anchored fraction (§2.2) — *natural to expand alongside the 9-genome comparison table in B1.*

---

## Notes

- This document lists analyses only. Writing/formatting/consistency/deposition items remain in
  `PUBLICATION_PLAN.md` and are not duplicated here.
- Review this inventory against your own knowledge of the paper and the 9 genomes in hand — flag
  anything Part A missed or Part B got wrong before treating it as a fixed checklist.
