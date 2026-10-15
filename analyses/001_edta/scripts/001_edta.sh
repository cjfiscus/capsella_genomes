#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --output=std/edta_%a_%A.stdout
#SBATCH --error=std/edta_%a_%A.stderr
#SBATCH --mail-user=fiscuscj@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --time=14-00:00:00
#SBATCH --job-name="edta"
#SBATCH -p koeniglab
#SBATCH --array=1-7

# software dependencies
## EDTA 2.3.0 (conda env EDTA_2.3.0)

source /opt/linux/rhel/8.x/x86_64/pkgs/miniconda3/py39_4.12.0/etc/profile.d/conda.sh
conda activate EDTA_2.3.0

# set -euo pipefail must come AFTER conda activate (see CLAUDE.md) —
# conda's activate functions are not strict-mode-safe and can silently
# leave PATH pointed at the base env under set -u/set -e.
set -euo pipefail
export PYTHONUNBUFFERED=1

# SET VARIABLES
LIST=../data/genomes.tsv
PREPDIR=$(realpath ../results/00_prep)
OUTBASE=$(realpath ../results/01_edta)
THREADS="${SLURM_CPUS_PER_TASK:?SLURM_CPUS_PER_TASK not set}"

#### PIPELINE #####
NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | cut -f1)
echo "[$NAME] starting EDTA (--sensitive 0)"

RUNDIR="$OUTBASE"/"$NAME"
mkdir -pv "$RUNDIR"
cd "$RUNDIR"

ln -sf "$PREPDIR"/"$NAME".fa "$NAME".fa
ln -sf "$PREPDIR"/"$NAME".cds.fa "$NAME".cds.fa

EDTA.pl \
    --genome "$NAME".fa \
    --cds "$NAME".cds.fa \
    --anno 1 \
    --sensitive 0 \
    --threads "$THREADS"

echo "[$NAME] EDTA complete"
