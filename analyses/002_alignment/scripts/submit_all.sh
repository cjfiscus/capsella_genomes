#!/bin/bash
# Submit the full 002_alignment pipeline as a dependency chain.
# Run from the scripts/ directory: bash submit_all.sh
# Requires conda env syri_1.8.2 (see README) — run test_env.sh first on a new setup.
set -euo pipefail

mkdir -pv std

JOB0=$(sbatch --parsable 000_prep_genomes.sh)
echo "000_prep_genomes: $JOB0 (array 1-9)"

# 001 runs in two waves (col 5 of genomes.tsv is the ref_mode): 'cbp' and 'fixed_cbp' genomes
# are classified against the normalized Cbp2-2 subgenomes, so they must wait
# for the diploid/fixed wave to finish.
WAVE1=$(awk -F'\t' '$5!~/cbp/  {printf "%s%d", (n++?",":""), NR}' ../data/genomes.tsv)
WAVE2=$(awk -F'\t' '$5~/cbp/  {printf "%s%d", (n++?",":""), NR}' ../data/genomes.tsv)

JOB1A=$(sbatch --parsable --dependency=afterok:"$JOB0" --array="$WAVE1" 001_chrom_normalize.sh)
echo "001_chrom_normalize wave 1 (diploid/fixed): $JOB1A (array $WAVE1)"

JOB1=$(sbatch --parsable --dependency=afterok:"$JOB1A" --array="$WAVE2" 001_chrom_normalize.sh)
echo "001_chrom_normalize wave 2 (cbp): $JOB1 (array $WAVE2)"

JOB2=$(sbatch --parsable --dependency=afterok:"$JOB1" 002_align_syri.sh)
echo "002_align_syri: $JOB2 (array 1-30)"

JOB3=$(sbatch --parsable --dependency=afterok:"$JOB2" 003_summarize.sh)
echo "003_summarize: $JOB3"

echo ""
echo "Chain submitted. Watch with: squeue -u $USER"
