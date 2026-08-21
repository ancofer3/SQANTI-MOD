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
    --bam_dir "$BAM_DIR" \
    --bam_suffix "$BAM_SUFFIX" \
    --bed_dir "$BED_OUTDIR" \
    --isoquant_dir "$ISOQUANT_DIR" \
    --out_dir "$TSV_OUTDIR" \
    --prob_lim "$PROB_LIM"