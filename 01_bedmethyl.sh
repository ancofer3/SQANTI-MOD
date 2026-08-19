#!/bin/bash
#SBATCH --job-name=bedmethyl
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --output=logs/bedmethyl_%j.out
#SBATCH --error=logs/bedmethyl_%j.err
#SBATCH --time=04:00:00
#SBATCH --qos=short

bams_path="../iPSCs/data/bams_primary/*primary.bam"
beds_path="data/beds"
mkdir -p $beds_path
for bam in $bams_path; do
	name=$(basename "$bam" _primary.bam)
	../iPSCs/dist_modkit_v0.6.4_cd85862/modkit pileup $bam \
		data/beds/${name}.bed \
		--modified-bases m5C m6A 2OmeC pseU 2OmeU inosine 2OmeA \
		--reference ../genomes/hg38/hg38.fa \
		--log logs/pileup.txt
done
