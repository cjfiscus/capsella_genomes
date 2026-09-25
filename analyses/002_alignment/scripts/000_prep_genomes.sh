#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=4G
#SBATCH --output=std/prep_%a_%A.stdout
#SBATCH --error=std/prep_%a_%A.stderr
#SBATCH --mail-user=fiscuscj@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH --time=0-02:00:00
#SBATCH --job-name="aln_prep"
#SBATCH -p koeniglab
#SBATCH --array=1-9

# software dependencies
## seqkit/2.4.0

module load seqkit/2.4.0

set -euo pipefail

# SET VARIABLES
LIST=../data/genomes.tsv
OUTDIR=../results/00_prep
MIN_CHR_LEN=10000000    # smallest kept sequence must be >= 10 Mb, else an unplaced scaffold slipped in

#### PIPELINE #####
mkdir -pv "$OUTDIR"

# parse genome list: name  fasta  n_chrom  group
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
NAME=$(echo "$LINE" | cut -f1)
FASTA_SRC=$(echo "$LINE" | cut -f2)
N_CHR=$(echo "$LINE" | cut -f3)

echo "[$NAME] fasta=$FASTA_SRC keep_longest=$N_CHR"

OUT_FA="$OUTDIR"/"$NAME".fa

# trim headers to first token and keep only the N longest sequences (the
# chromosomes; everything else is unplaced scaffolds / organelles). seqkit
# reads .gz transparently. The top-N list is built via files, not a
# `... | seqkit head` pipe: head exiting early SIGPIPEs the upstream seqkit
# and pipefail then fails the job.
TEMP_DIR="$TMPDIR"/prep_${SLURM_ARRAY_TASK_ID}   # $TMPDIR is the per-job scratch dir (see CLAUDE.md)
mkdir -pv "$TEMP_DIR"
trap 'rm -rf "$TEMP_DIR"' EXIT

seqkit fx2tab -n -l -i "$FASTA_SRC" | sort -k2,2nr > "$TEMP_DIR"/lengths.sorted.tsv
head -n "$N_CHR" "$TEMP_DIR"/lengths.sorted.tsv | cut -f1 > "$TEMP_DIR"/keep.txt
seqkit seq -i "$FASTA_SRC" | seqkit grep -f "$TEMP_DIR"/keep.txt > "$OUT_FA"

N_OUT=$(grep -c '^>' "$OUT_FA")
if [ "$N_OUT" -ne "$N_CHR" ]; then
    echo "[$NAME] ERROR: expected $N_CHR sequences, got $N_OUT (source has fewer than N_CHR sequences?)" >&2
    exit 1
fi

MIN_LEN=$(seqkit fx2tab -l -n "$OUT_FA" | cut -f2 | sort -n | sed -n 1p)
echo "[$NAME] kept $N_OUT sequences, smallest = $MIN_LEN bp"
if [ "$MIN_LEN" -lt "$MIN_CHR_LEN" ]; then
    echo "[$NAME] ERROR: smallest kept sequence ($MIN_LEN bp) < $MIN_CHR_LEN — n_chrom in genomes.tsv is probably wrong" >&2
    exit 1
fi

seqkit stats -T "$OUT_FA"
echo "[$NAME] prep complete: $OUT_FA"
