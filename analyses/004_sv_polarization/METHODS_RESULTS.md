# 004_sv_polarization — methods and results

Inventory item B1 #7. Preliminary and descriptive; no significance tests. Run details are in `README`.

## Methods

**Question.** For each structural variant between two genomes, which genome carries the change (an insertion on one lineage or a deletion on the other)? Pairwise SyRI calls are reference-dependent, so each event is polarized against outgroup genomes.

**Focal pairs.** The 30 within-set pairs from `002_alignment` (minimap2 `asm10` + SyRI 1.8.2): diploid × diploid, diploid × Cbp subgenome, and Cbp × Cbp. Each Cbp assembly is split into its Co and Cr subgenomes. Events are SyRI INS/DEL ≥ 50 bp and INV/TRANS/INVTR ≥ 50 bp. DUP/INVDP are counted but not polarized.

**Outgroups.** Diploid × diploid and diploid × Cbp pairs use the diploids of the other set (Cr/Cg for Co-set pairs, Co for Cr-set pairs). Cbp × Cbp pairs use the within-set diploids. Within-set Co diploids are not used as outgroups to Cbp-Co (about 0.9 SNP/kb among diploids vs about 5 SNP/kb to Cbp-Co, see `003_diploid_baseline`). This differs from the outgroups listed in the inventory for #7.

**New alignments** (asm10, `002_align_syri.sh` with `LIST`/`OUTROOT`/`PRESET` overrides): 9 Co-diploid × Cr/Cg-diploid pairs and 18 Cbp-subgenome × opposite-set diploid pairs. A Co39 vs Cr145 pilot gave 105 vs 92 Mb syntenic and 24 vs 38 Mb unaligned for asm10 vs asm20, so asm10 was used throughout.

**Polarization (`scripts/001_polarize.py`).** Each event is expressed in the coordinates of a frame genome X (A or B of the pair) and tested against the X-vs-outgroup alignment.
- *Carrier = X:* if at least 70% of the interval lies in an X-only indel of the outgroup alignment, the call is an insertion on X. If at least 70% lies in aligned blocks, the call is a deletion on the other genome.
- *Carrier = other genome:* an outgroup insertion of similar size (×0.5–2) within 50 bp gives a deletion on X. Aligned flanks with no indel give an insertion on the other genome.
- *INV/TRANS/INVTR:* an outgroup syntenic over the interval puts the change on the other genome. An outgroup with the same SV type puts it on X.
- HDR blocks are not counted as "outgroup has it".
- The consensus per frame is the unique plurality over the 3 outgroups. The two frames are then reconciled as `both`, `A_frame_only`, `B_frame_only`, `conflict` or `none`.
- Both frames are used because a single frame skewed calls about 4:1 toward the non-frame genome even for symmetric diploid pairs. Two-frame calls give about 1:1.

**Plots** (`scripts/002_plots.R`, R, ggplot2) are drawn from `results/polarization_summary.tsv`:
- Decided indels (`final` is ins/del) are assigned to the lineage of the carrying genome.
- `confidence == both` counts as a both-frame call, and anything else decided as single-frame.
- Lineage totals are summed over pairs, which share genomes and are not independent.

## Results

**Insertion vs deletion by lineage** (indels ≥ 50 bp, diploid × Cbp pairs, decided calls; `results/02_plots/indel_by_lineage.png`):

| lineage | deletions | insertions | ins:del |
|---|---|---|---|
| Cbp-Co | 4787 | 4341 | 0.91 |
| Co diploids | 4258 | 7451 | 1.75 |
| Cbp-Cr | 4360 | 3460 | 0.79 |
| Cr/Cg diploids | 5899 | 6802 | 1.15 |

- Both Cbp subgenomes are deletion-biased relative to their diploid counterparts.
- **Per pair** (`indel_ratio_per_pair.png`): the Cbp-Co ratio is bimodal. Two pairs are at about 0.35 and the rest at about 1.5, so the pooled 0.91 hides structure. Mosher pairs do not separate cleanly from the others. This is not yet explained.
- **Call quality** (`outcome_mix.png`): of 81,336 diploid × Cbp indel events, 38,103 stay ambiguous and 1,875 conflict. About 32% of all events (about 62% of decided calls) are single-frame. Diploid × diploid indels are the least resolved (11,224 of 18,236 ambiguous). INV/TRANS/INVTR are mostly ambiguous or single-frame in every class.

## Caveats
- The deletion bias is not tested for significance, by size class, per pair, excluding Mosher pairs, or normalized for assembly size.
- The diploid totals include the larger CoMosher and CgMosher assemblies (extra TE insertions).
- Filter on `confidence` before using individual calls. Conflicts (about 5%) are excluded.
- Cross-species Co × Cr/Cg pairs are not polarized (no outgroup).
- SV calls are unreviewed and TRANS/DUP counts depend on the alignment preset.
