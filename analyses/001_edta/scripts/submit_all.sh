#!/bin/bash
# Submit the full 001_edta pipeline as a dependency chain.
# Run from the scripts/ directory: bash submit_all.sh
set -euo pipefail

mkdir -pv std

JOB0=$(sbatch --parsable 000_prep_genomes.sh)
echo "000_prep_genomes: $JOB0 (array 1-7)"

JOB1=$(sbatch --parsable --dependency=afterok:"$JOB0" 001_edta.sh)
echo "001_edta: $JOB1 (array 1-7)"

JOB2=$(sbatch --parsable --dependency=afterok:"$JOB1" 002_panlib_reanno.sh)
echo "002_panlib_reanno: $JOB2"

echo ""
echo "Chain submitted. Watch with: squeue -u $USER"
