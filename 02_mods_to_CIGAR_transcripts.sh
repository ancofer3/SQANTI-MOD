#!/bin/bash
#SBATCH --job-name=mods_to_cigar
#SBATCH --cpus-per-task=50
#SBATCH --mem=250G
#SBATCH --output=logs/mods_to_CIGAR_%j.out
#SBATCH --error=logs/mods_to_CIGAR_%j.err
#SBATCH --time=04:00:00
#SBATCH --qos=short

source config.sh
module load anaconda
conda activate jae_env

python -u saca_tags_transcript_models.py \
    --bed "$BED" \
    --bam "$BAM" \
    --gtf "$GTF" \
    --assoc "$ASSOC" \
    --out_tsv "$OUT_TSV" \
    --prob_lim "$PROB_LIM"

