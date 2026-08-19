#!/bin/bash
#SBATCH --job-name=mods_to_cigar_transcripts
#SBATCH --cpus-per-task=50
#SBATCH --mem=250G
#SBATCH --output=logs/mods_to_CIGAR_transcripts_%j.out
#SBATCH --error=logs/mods_to_CIGAR_transcripts_%j.err
#SBATCH --time=04:00:00
#SBATCH --qos=short

module load anaconda
conda activate jae_env
python -u saca_tags_transcript_models.py
