#!/bin/bash
#
#SBATCH --time=8:00:00
#SBATCH --ntasks=20
#SBATCH --mem=50g
#SBATCH --tmp=120g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pande250@umn.edu 

#loading fastqc and multiqc
module load fastqc/0.11.5

# Define output directory
OUTPUT_DIR="/scratch.global/pande250/GCD-8141/Project-3-rnaseq/fastqc_output"

# Run FastQC on all FASTQ files in the current directory
echo "Running FastQC..."
fastqc *.fastq -o "$OUTPUT_DIR"

# Move to output directory and run MultiQC there
echo "Running MultiQC..."
cd "$OUTPUT_DIR" || exit
multiqc . -o "$OUTPUT_DIR"
echo "✅ All done! Results saved in $OUTPUT_DIR"

