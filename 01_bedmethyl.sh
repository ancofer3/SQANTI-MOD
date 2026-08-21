#!/bin/bash
#SBATCH --job-name=bedmethyl
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --output=logs/bedmethyl_%j.out
#SBATCH --error=logs/bedmethyl_%j.err
#SBATCH --time=04:00:00
#SBATCH --qos=short
source config.sh
mkdir -p $beds_path "BED_OUTDIR" "LOG_DIR"

for bam in "$BAM_DIR"/*"$BAM_SUFFIX"; do
	name=$(basename "$bam" "$BAM_SUFFIX")
	"$MODKIT_BIN" pileup $bam \
		"$BED_OUTDIR"/${name}.bed \
		--modified-bases $MODS\
		--reference "$REF_GENOME"\
		--log "$LOG_DIR"/pileup_${name}.txt
	# Conservamos sitios con cobertura minima genómica >= 5 para para quitar ruido
	awk -v min_cov="$MIN_GENOMIC_COV" '($10 >= min_cov) || /^#/' "$BED_OUTDIR/${name}.bed" > "$BED_OUTDIR/${name}_filtered.bed"
done
