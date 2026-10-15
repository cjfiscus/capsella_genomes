#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=4G
#SBATCH --output=std/prep_%a_%A.stdout
#SBATCH --error=std/prep_%a_%A.stderr
#SBATCH --mail-user=fiscuscj@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH --time=0-04:00:00
#SBATCH --job-name="prep"
#SBATCH -p koeniglab
#SBATCH --array=1-7

# software dependencies
## seqkit/2.4.0
## gffread/0.12.7

module load seqkit/2.4.0
module load gffread/0.12.7

set -euo pipefail

# SET VARIABLES
LIST=../data/genomes.tsv
OUTDIR=../results/00_prep
# $TMPDIR is auto-provisioned per job by SLURM (points at
# /scratch/cfisc004/$SLURM_JOB_ID, already owned by us). The top-level
# /scratch/cfisc004 dir itself is root-owned — mkdir'ing a custom-named
# directory directly under it fails with Permission denied. Always use
# $TMPDIR (or a subdir under it), never hand-roll a /scratch/$USER/... path.
TEMP_DIR="$TMPDIR"/prep_${SLURM_ARRAY_TASK_ID}

#### PIPELINE #####
mkdir -pv "$OUTDIR"
mkdir -pv "$TEMP_DIR"
trap 'rm -rf "$TEMP_DIR"' EXIT

# parse genome list
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
NAME=$(echo "$LINE" | cut -f1)
FASTA_SRC=$(echo "$LINE" | cut -f2)
GFF_SRC=$(echo "$LINE" | cut -f3)
SEQID_MAP=$(echo "$LINE" | cut -f4)

echo "[$NAME] fasta=$FASTA_SRC gff=$GFF_SRC seqid_map=$SEQID_MAP"

OUT_FA="$OUTDIR"/"$NAME".fa
OUT_CDS="$OUTDIR"/"$NAME".cds.fa

# 1) trim FASTA headers to first token (seqkit reads .gz transparently)
seqkit seq -i "$FASTA_SRC" > "$OUT_FA"

N_SRC=$(seqkit stats -T "$FASTA_SRC" | tail -n1 | cut -f4)
N_OUT=$(seqkit stats -T "$OUT_FA" | tail -n1 | cut -f4)
echo "[$NAME] sequences: source=$N_SRC output=$N_OUT"
if [ "$N_SRC" != "$N_OUT" ]; then
    echo "[$NAME] ERROR: sequence count mismatch after header trimming" >&2
    exit 1
fi

# 2) stage GFF (decompress if needed)
GFF_PLAIN="$TEMP_DIR"/"$NAME".gff
if [[ "$GFF_SRC" == *.gz ]]; then
    zcat "$GFF_SRC" > "$GFF_PLAIN"
else
    cat "$GFF_SRC" > "$GFF_PLAIN"
fi

# 3) remap seqids in GFF col1 if a map is provided
GFF_FOR_CDS="$GFF_PLAIN"
if [ "$SEQID_MAP" != "-" ]; then
    if [ ! -s "$SEQID_MAP" ]; then
        echo "[$NAME] ERROR: seqid map $SEQID_MAP not found" >&2
        exit 1
    fi
    GFF_MAPPED="$TEMP_DIR"/"$NAME".mapped.gff
    awk -F'\t' -v OFS='\t' '
        NR==FNR { map[$1]=$2; next }
        /^#/    { print; next }
        {
            if ($1 in map) $1 = map[$1]
            print
        }
    ' "$SEQID_MAP" "$GFF_PLAIN" > "$GFF_MAPPED"
    GFF_FOR_CDS="$GFF_MAPPED"
fi

# 3b) keep only GFF features on sequences present in this FASTA (needed for
# the Cbp2-2 subgenomes, which share the full-genome GFF; no-op elsewhere)
seqkit seq -n -i "$OUT_FA" > "$TEMP_DIR"/"$NAME".seqids
GFF_FILT="$TEMP_DIR"/"$NAME".filt.gff
awk -F'\t' -v OFS='\t' 'NR==FNR { keep[$1]=1; next } /^#/ || ($1 in keep)' \
    "$TEMP_DIR"/"$NAME".seqids "$GFF_FOR_CDS" > "$GFF_FILT"
GFF_FOR_CDS="$GFF_FILT"

# 4) extract CDS sequences against the trimmed genome
gffread -x "$OUT_CDS" -g "$OUT_FA" "$GFF_FOR_CDS"

N_CDS=$(grep -c '^>' "$OUT_CDS" || true)
echo "[$NAME] CDS sequences extracted: $N_CDS"
if [ "$N_CDS" -eq 0 ]; then
    echo "[$NAME] ERROR: no CDS sequences extracted — check seqid match between GFF and FASTA" >&2
    exit 1
fi

echo "[$NAME] prep complete: $OUT_FA , $OUT_CDS"
