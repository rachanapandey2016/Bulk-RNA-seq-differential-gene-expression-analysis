#!/bin/bash
#
#SBATCH --time=10:00:00
#SBATCH --ntasks=24
#SBATCH --mem=80g
#SBATCH --tmp=120g
#SBATCH -o star_indexed.%j.o
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pande250@umn.edu 

#Load star
module load star/2.7.11b

genomedir=/home/gcd8141/public/aazaidi/genome/
gtfdir=/home/gcd8141/public/aazaidi/gtf/

# Output index directory
outdir="/scratch.global/pande250/GCD-8141/Project-3-rnaseq/star_indexed_genome"

# Run STAR genomeGenerate
STAR --runThreadN 4 \
     --runMode genomeGenerate \
     --genomeDir "$outdir" \
     --genomeFastaFiles "${genomedir}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" \
     --sjdbGTFfile "${gtfdir}/gencode.v36.basic.annotation.gtf" \
     --sjdbOverhang 75
