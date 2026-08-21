#!/bin/bash
# ==========================================
# SQANTI-MOD pipeline configuration file
# ==========================================

# 1. Routes to tools and reference data
MODKIT_BIN="../iPSCs/dist_modkit_v0.6.4_cd85862/modkit"
REF_GENOME="../genomes/hg38/hg38.fa"

# 2. Routes to input data
BAM_DIR="../iPSCs/data/bams_primary"
BAM_SUFFIX="_primary.bam" # Sufijo para limpiar el nombre del archivo
ISOQUANT_DIR="../iPSCs/reconstruction/isoquant"
SQANTI_DIR="../iPSCs" # Directorio padre donde buscar los classification.txt

# 3. Routes to output data
BED_OUTDIR="data/beds"
TSV_OUTDIR="data/tsvs_transcritos"
CLASS_OUTDIR="data/classifications_mod"
LOG_DIR="logs"

# 4. Filtering and mod-specific parameters
MODS="m5C m6A 2OmeC pseU 2OmeU inosine 2OmeA"
MIN_GENOMIC_COV=5
PROB_LIM=0.95
MIN_TX_COV=20
MIN_TX_OCC=0.2