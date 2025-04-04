#!/bin/bash
#
#SBATCH --job-name=salmon_quant
#SBATCH --time=06:00:00
#SBATCH --ntasks=8
#SBATCH --mem=60G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pande250@umn.edu
#SBATCH -o salmon_quant.%j.out

#Load salmon from MSI
module load salmon/1.10.0

#Activate the conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq

#Define path to input directoriess
GENOME="/home/gcd8141/public/aazaidi/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
GTF="/home/gcd8141/public/aazaidi/gtf/gencode.v36.basic.annotation.gtf"
OUTDIR="/scratch.global/pande250/GCD-8141/Project-3-rnaseq/salmon_gtf"

mkdir -p "$OUTDIR"

#Salmon requires an annotated transcriptome. We can create transcriptome FASTA from GTF
gffread -w "${OUTDIR}/GRCh38_no_alt_analysis_set_gencode.v36.transcripts.fa" \
        -g "$GENOME" \
        "$GTF"

# this is optional but we will also create  a bed file for eQTL analysis needed by tensorQTL
gffread -g "$GENOME" \
        -G "$GTF" \
        -t "${OUTDIR}/gencode.v36.basic.annotation.fa" \
        -o "${OUTDIR}/gencode.v36.basic.annotation.bed" \
        --bed

# Run Salmon quantification 
TRANSCRIPTS="${OUTDIR}/GRCh38_no_alt_analysis_set_gencode.v36.transcripts.fa"
SALMON_OUT_DIR="/scratch.global/pande250/GCD-8141/Project-3-rnaseq/salmon_quant"
STAR_BAM_DIR="/scratch.global/pande250/GCD-8141/Project-3-rnaseq/star_output"

mkdir -p "$SALMON_OUT_DIR"

cd /scratch.global/pande250/GCD-8141/Project-3-rnaseq/star_output

for bam in *_Aligned.toTranscriptome.out.bam; do
    sample=$(basename "$bam" _Aligned.toTranscriptome.out.bam)
    echo "Quantifying $sample..."

    salmon quant \
      -t "$TRANSCRIPTS" \
      --libType A \
      -a "$bam" \
      -o "${SALMON_OUT_DIR}/${sample}.salmon_quant" \
      --gcBias --seqBias
    echo "Completed: $sample"

done
