#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --output=std/panlib_%j.stdout
#SBATCH --error=std/panlib_%j.stderr
#SBATCH --mail-user=fiscuscj@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --time=7-00:00:00
#SBATCH --job-name="panlib"
#SBATCH -p koeniglab

# software dependencies
## EDTA 2.3.0 / panEDTA.sh (conda env EDTA_2.3.0)
#
# Runs panEDTA end-to-end: builds the joint (pan) TE library from all 7
# per-genome EDTA runs (results/01_edta/<name>/), then reannotates every
# genome with that joint library. This is a single sequential job because
# panEDTA.sh itself loops over genomes internally (not parallel array-safe) —
# see analyses/001_edta/README for why this isn't split into separate steps:
# EDTA.pl's --step final needs several exact intermediate files from the
# per-genome run (RM2.fa, EDTA.fa.stg1, intact.fa.cln, intact.raw.gff3) that
# panEDTA.sh symlinks in itself; reimplementing that by hand risks silently
# diverging from upstream.

source /opt/linux/rhel/8.x/x86_64/pkgs/miniconda3/py39_4.12.0/etc/profile.d/conda.sh
conda activate EDTA_2.3.0

# set -euo pipefail must come AFTER conda activate (see CLAUDE.md).
set -euo pipefail
export PYTHONUNBUFFERED=1

# SET VARIABLES
LIST=$(realpath ../data/genomes.tsv)
EDTADIR=$(realpath ../results/01_edta)
OUTDIR=$(realpath ../results/02_panlib)
THREADS="${SLURM_CPUS_PER_TASK:?SLURM_CPUS_PER_TASK not set}"

#### PIPELINE #####
mkdir -pv "$OUTDIR"
cd "$OUTDIR"

# build the 2-column genome list panEDTA expects: <genome.fa> <TAB> <cds.fa>
# pointing at each genome's own EDTA run dir, where EDTA.pl already wrote
# <name>.fa.mod.EDTA.* outputs that panEDTA.sh will detect and reuse.
> genomes.list
while IFS=$'\t' read -r NAME FASTA_SRC GFF_SRC SEQID_MAP; do
    printf '%s\t%s\n' \
        "$EDTADIR"/"$NAME"/"$NAME".fa \
        "$EDTADIR"/"$NAME"/"$NAME".cds.fa \
        >> genomes.list
done < "$LIST"

cat genomes.list

# -a 1: build the pan library AND reannotate every genome with it
# -o 1: allow panEDTA to overwrite/relink existing results in this dir (reruns)
# -f 3: keep candidate TEs with >=3 full-length copies in a genome
panEDTA.sh -g genomes.list -a 1 -o 1 -f 3 -t "$THREADS"

echo "panEDTA library: $OUTDIR/genomes.list.panEDTA.TElib.fa"
echo "Per-genome reannotation under $OUTDIR/<name>.fa.mod.EDTA.anno/"
