#!/bin/bash
#SBATCH --job-name=collapse_cigars
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=logs/collapse_%j.out
#SBATCH --error=logs/collapse_%j.err
#SBATCH --time=04:00:00
#SBATCH --qos=short

source config.sh
module load anaconda
conda activate jae_env

python -u collapse_to_tx_and_merge.py \
    --tsv_in "$OUT_TSV" \
    --sqanti_class "$SQANTI_CLASS" \
    --out_tsv "$CLASS_OUTDIR" \
    --min_cov "$MIN_TX_COV" \
    --min_occ "$MIN_TX_OCC"
