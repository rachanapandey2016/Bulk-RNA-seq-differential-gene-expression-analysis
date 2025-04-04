#!/bin/bash
#
#SBATCH --time=8:00:00
#SBATCH --ntasks=20
#SBATCH --mem=50g
#SBATCH --tmp=120g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pande250@umn.edu 

# Activate conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq  # Replace with your actual environment name

# Output directory
OUTDIR="/scratch.global/pande250/GCD-8141/Project-3-rnaseq/fastp_output"

# Loop through all .fastq files in the current directory
for file in *.fastq; do
    # Get the sample name without extension
    sample=$(basename "$file" .fastq)

    echo "Processing sample: $sample"

    fastp \
        -i "$file" \
        -o "$OUTDIR/${sample}.trimmed.fastq" \
        -l 25 \
        -j "$OUTDIR/${sample}.fastp.json" \
        -h "$OUTDIR/${sample}.fastp.html"

    echo "Finished processing $sample"
done
