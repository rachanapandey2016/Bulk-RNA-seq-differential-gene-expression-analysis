#!/bin/bash
#
#SBATCH --time=10:00:00
#SBATCH --ntasks=24
#SBATCH --mem=80g
#SBATCH -o star_alignment.%j.o
#SBATCH --tmp=120g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pande250@umn.edu 

#Load STAR
module load star/2.7.11b

genomedir=/scratch.global/pande250/GCD-8141/Project-3-rnaseq/star_indexed_genome

# Output index directory
outdir="/scratch.global/pande250/GCD-8141/Project-3-rnaseq/star_output"

# Run STAR 
# Loop through all trimmd fast.gz files in current directory
for fq in *.trimmed.fastq; do
    # Extract sample name (remove -trimmed.fastq.gz)
    sample=$(basename "$fq" .trimmed.fastq)

    echo "Running STAR for sample: $sample"

    STAR --runThreadN 24 \
         --genomeDir "$genomedir" \
         --readFilesIn "$fq" \
         --outFileNamePrefix "$outdir/${sample}_" \
         --outSAMtype BAM Unsorted \
         --outFilterType BySJout \
         --alignSJoverhangMin 8 \
         --outFilterMultimapNmax 20 \
         --alignSJDBoverhangMin 1 \
         --outFilterMismatchNmax 999 \
         --outFilterMismatchNoverReadLmax 0.04 \
         --alignIntronMin 20 \
         --alignIntronMax 1000000 \
         --alignMatesGapMax 1000000 \
         --quantMode TranscriptomeSAM \
         --outSAMattributes NH HI AS NM MD

    echo "Finished STAR for sample: $sample"

done


