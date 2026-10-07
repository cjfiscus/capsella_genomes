#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --output=std/cgmosher_final_%j.stdout
#SBATCH --error=std/cgmosher_final_%j.stderr
#SBATCH --mail-user=fiscuscj@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --time=1-00:00:00
#SBATCH --job-name="cg_final"
#SBATCH -p koeniglab

# software dependencies
## EDTA 2.3.0 (conda env EDTA_2.3.0)
#
# Retry of panEDTA's per-genome "structural" reannotation for CgMosher only.
# In job 29303452 the CDS-cleanup RepeatMasker/BLAST stage of EDTA.pl
# --step final failed for CgMosher (empty intact.fa.cln2), so no
# CgMosher.fa.mod.EDTA.anno/TEanno.sum was produced even though panEDTA.sh
# exited 0. This runs the exact command panEDTA.sh uses (step 3, structural)
# without rebuilding the joint library. Inputs (genomes.list.panEDTA.TElib.fa,
# CgMosher.fa.mod.panEDTA.out, CDS) already exist in results/02_panlib.

source /opt/linux/rhel/8.x/x86_64/pkgs/miniconda3/py39_4.12.0/etc/profile.d/conda.sh
conda activate EDTA_2.3.0

# set -euo pipefail must come AFTER conda activate (see CLAUDE.md).
set -euo pipefail
export PYTHONUNBUFFERED=1

# SET VARIABLES
NAME=CgMosher
OUTDIR=$(realpath ../results/02_panlib)
THREADS="${SLURM_CPUS_PER_TASK:?SLURM_CPUS_PER_TASK not set}"

#### PIPELINE #####
cd "$OUTDIR"

for f in genomes.list.panEDTA.TElib.fa "$NAME".fa.mod.panEDTA.out "$NAME".cds.fa "$NAME".fa; do
    [ -s "$f" ] || { echo "ERROR: missing input $f" >&2; exit 1; }
done

# keep the failed attempt for reference; EDTA.pl would otherwise reuse its
# stale intermediates
if [ -d "$NAME".fa.mod.EDTA.final ]; then
    mv -v "$NAME".fa.mod.EDTA.final "$NAME".fa.mod.EDTA.final.failed_29303452
fi

EDTA.pl --genome "$NAME".fa -t "$THREADS" --step final --anno 1 \
    --curatedlib genomes.list.panEDTA.TElib.fa \
    --cds "$NAME".cds.fa \
    --rmout "$NAME".fa.mod.panEDTA.out

[ -s "$NAME".fa.mod.EDTA.anno/"$NAME".fa.mod.EDTA.TEanno.sum ] \
    || { echo "ERROR: no TEanno.sum produced for $NAME" >&2; exit 1; }

echo "[$NAME] final-step retry complete"
