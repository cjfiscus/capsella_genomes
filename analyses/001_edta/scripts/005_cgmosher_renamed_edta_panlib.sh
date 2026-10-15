#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --output=std/cgmosher_ren_%j.stdout
#SBATCH --error=std/cgmosher_ren_%j.stderr
#SBATCH --mail-user=fiscuscj@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --time=7-00:00:00
#SBATCH --job-name="cg_ren"
#SBATCH -p koeniglab

# software dependencies
## EDTA 2.3.0 (conda env EDTA_2.3.0)
#
# CgMosher with shortened scaffold IDs (see 004_cgmosher_rename.sh), run in
# results/03_cgmosher_renamed. Mirrors what 001_edta.sh + panEDTA.sh step 3 do
# for one genome, reusing the EXISTING joint library (no library rebuild):
#   1) standalone EDTA.pl   2) RepeatMasker with joint lib   3) --step final
# Each stage is skipped if its output exists, so the job can be resubmitted.

source /opt/linux/rhel/8.x/x86_64/pkgs/miniconda3/py39_4.12.0/etc/profile.d/conda.sh
conda activate EDTA_2.3.0

# set -euo pipefail must come AFTER conda activate (see CLAUDE.md).
set -euo pipefail
export PYTHONUNBUFFERED=1

# SET VARIABLES
NAME=CgMosher
RUNDIR=$(realpath ../results/03_cgmosher_renamed)
CDS=$(realpath ../results/00_prep/"$NAME".cds.fa)
LIB=$(realpath ../results/02_panlib/genomes.list.panEDTA.TElib.fa)
THREADS="${SLURM_CPUS_PER_TASK:?SLURM_CPUS_PER_TASK not set}"
RM_THREADS=$((THREADS / 4))   # as in panEDTA.sh

#### PIPELINE #####
cd "$RUNDIR"
for f in "$NAME".fa "$CDS" "$LIB"; do
    [ -s "$f" ] || { echo "ERROR: missing input $f" >&2; exit 1; }
done
MAXID=$(grep '^>' "$NAME".fa | awk '{l=length($0)-1; if(l>m)m=l} END{print m}')
[ "$MAXID" -le 30 ] || { echo "ERROR: scaffold IDs too long ($MAXID)" >&2; exit 1; }

ln -sf "$CDS" "$NAME".cds.fa

# 1) standalone EDTA (structural + homology, --sensitive 0 as in 001_edta.sh)
if [ ! -s "$NAME".fa.mod.EDTA.TElib.fa ]; then
    echo "[$NAME] standalone EDTA"
    EDTA.pl --genome "$NAME".fa --cds "$NAME".cds.fa --anno 1 --sensitive 0 \
        --threads "$THREADS"
fi
[ -s "$NAME".fa.mod.EDTA.TElib.fa ] || { echo "ERROR: standalone EDTA failed" >&2; exit 1; }

# 2) RepeatMasker with the joint library (panEDTA.sh step 3, homology)
if [ ! -s "$NAME".fa.mod.panEDTA.out ]; then
    echo "[$NAME] RepeatMasker with joint library"
    ln -sf "$NAME".fa.mod "$NAME".fa.mod.panEDTA
    RepeatMasker -e ncbi -pa "$RM_THREADS" -q -div 40 -lib "$LIB" -cutoff 225 \
        -gff "$NAME".fa.mod.panEDTA >/dev/null
fi
[ -s "$NAME".fa.mod.panEDTA.out ] || { echo "ERROR: no panEDTA.out" >&2; exit 1; }
perl -i -nle 's/\s+DNA\s+/\tDNA\/unknown\t/; print $_' "$NAME".fa.mod.panEDTA.out

# 3) structural reannotation with the joint library (panEDTA.sh step 3)
echo "[$NAME] EDTA --step final"
EDTA.pl --genome "$NAME".fa -t "$THREADS" --step final --anno 1 \
    --curatedlib "$LIB" --cds "$NAME".cds.fa \
    --rmout "$NAME".fa.mod.panEDTA.out

SUM="$NAME".fa.mod.EDTA.anno/"$NAME".fa.mod.EDTA.TEanno.sum
[ -s "$SUM" ] || { echo "ERROR: no TEanno.sum produced for $NAME" >&2; exit 1; }
echo "[$NAME] complete: $RUNDIR/$SUM"
