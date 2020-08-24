#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=16G
#SBATCH --output=std/gvcf%j.stdout
#SBATCH --error=std/gvcf%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="gvcf"
#SBATCH -p intel
#SBATCH --array=111

# software dependencies
# gatk 4.1.8.1

# SET VARIABLES
REFERENCE=/rhome/cfisc004/shared/GENOMES/CAPSELLA/FINISHED_v2/data/assemblies/Cr145_plus_manual_correction.fasta
SPP=Cr
RESULTS=/rhome/cfisc004/bigdata/projects/capsella_genomes/results
SEQLIST=../../data/Capsella_WGS_Sample_List.txt
#TEMP_DIR=./
TEMP_DIR=/scratch/cfisc004
THREADS=4

#### PIPELINE #####
# parse sample file
#NAME=$(head -n "$SLURM_ARRAY_TASK_ID" "$SEQLIST" | tail -n1 | cut -f2)
NAME=$(head -n "$SLURM_ARRAY_TASK_ID" "$SEQLIST" | tail -n1 | cut -f4)
BAMIN="$RESULTS"/alignments/"$SPP"/"$NAME".bam

# mk temp
TEMP_DIR="$TEMP_DIR"/"$NAME"_"$RANDOM"
mkdir -pv "$TEMP_DIR"
cd "$TEMP_DIR"

####
# mark duplicates
gatk --java-options "-Xmx32G" MarkDuplicatesSpark \
            -I "$BAMIN" \
            -O marked_duplicates.bam \
	    --conf 'spark.local.dir=./'

# haplotype caller
gatk --java-options "-Xmx32g" HaplotypeCaller  \
   -R "$REFERENCE" \
   -I marked_duplicates.bam \
   -O "$RESULTS"/gvcf/"$SPP"/"$NAME".g.vcf.gz \
   -ERC GVCF

# clean up temp
rm -r $TEMP_DIR
