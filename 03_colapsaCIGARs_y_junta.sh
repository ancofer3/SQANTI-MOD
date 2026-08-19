#!/bin/bash
#SBATCH --job-name=collapse_cigars
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=logs/collapse_cigars_mergea_%j.out
#SBATCH --error=logs/collapse_cigars_mergea_%j.err
#SBATCH --time=04:00:00
#SBATCH --qos=short
module load anaconda
conda activate jae_env
python -u resume_Cigars_y_mergea.py
