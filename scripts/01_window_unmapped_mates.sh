#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# 01_window_unmapped_mates.sh
#
# For a given scaffold and range of coordinates, slide a window and:
#   1) collect read names in a pre-filtered BAM (e.g. mates mapped in region)
#   2) from the full BAM, select reads with these names whose mates are UNMAPPED
#   3) export them as FASTA
#
# This was used to generate:
#   - unmapped mates in the interesting region
#   - unmapped mates in random regions (negative controls)
###############################################################################

########################
# CONFIGURATION
########################

# Directory with BAM files
WORK_DIR=/home/svetlana/clu_scratch/202305101343data/Parus_major_WT_HiC/new_map

# Full Hi-C BAM
FULL_BAM=PM_WT.new.bam

# BAM containing reads whose mates map to various regions of interest
MATE_UNMAPPED_BAM=PM_WT.mate_unmapped.bam

# Scaffold to scan
SCAFFOLD=HiC_scaffold_6

# Windowing parameters (start, end, step, final coordinate)
WINDOW_START=58390000
WINDOW_END=58590000
WINDOW_STEP=100000
WINDOW_STOP=71990400   # stop when start >= WINDOW_STOP

# Output directory
OUT_DIR=Parus_major
mkdir -p "${OUT_DIR}"

########################
# MAIN LOOP
########################

st=${WINDOW_START}
en=${WINDOW_END}
fname=$(basename "${MATE_UNMAPPED_BAM}" .bam)

while [ "${st}" -lt "${WINDOW_STOP}" ]; do
  echo "region: ${SCAFFOLD}:${st}-${en}"

  # 1) Collect read names mapping in the current window in the pre-filtered BAM
  samtools view "${WORK_DIR}/${MATE_UNMAPPED_BAM}" "${SCAFFOLD}:${st}-${en}"         | cut -f1         | sort         | uniq         > "${OUT_DIR}/${fname}.${st}-${en}.reads"

  # 2) From the full BAM, select reads with these names whose mates are unmapped
  #    (-f 4: read itself unmapped; adapt flags if needed)
  samtools view         -N "${OUT_DIR}/${fname}.${st}-${en}.reads"         -f 4         -o "${OUT_DIR}/${fname}.${st}-${en}.bam"         "${WORK_DIR}/${FULL_BAM}"

  # 3) Export to FASTA
  samtools view "${OUT_DIR}/${fname}.${st}-${en}.bam"         | awk '{ printf(">%s\n%s\n", $1, $10) }'         > "${OUT_DIR}/${fname}.${st}-${en}.bam.fa"

  echo "created: ${OUT_DIR}/${fname}.${st}-${en}.bam.fa"

  # advance window
  st=$(( st + WINDOW_STEP ))
  en=$(( en + WINDOW_STEP ))
done
