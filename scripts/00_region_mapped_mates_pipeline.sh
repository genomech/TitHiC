#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# 00_region_mapped_mates_pipeline.sh
#
# Example pipeline for a single genomic interval:
#   - select reads in a given coordinate range
#   - keep a subset of mates according to SAM flags
#   - convert to FASTA
#   - run Jellyfish k-mer counting
#
# This was used as a positive control (mapped mates in the interesting region).
###############################################################################

########################
# CONFIGURATION
########################

# Genomic interval (SAM coordinates, 1-based).
REGION_START=522979400
REGION_END=535700607

# Number of lines to keep for the “normalized” subset
SUBSET_LINES=138055

# Samtools threads
THREADS=23

# Jellyfish parameters
JF_HASH_SIZE=60M
JF_KMER_LEN_MAIN=10   # e.g. 10-mers
JF_KMER_LEN_ALT=7     # optional second k-mer size

# Input SAM without header and separate header
SAM_BODY=Parus_major1.1_HiC.sam
SAM_HEADER=Parus_major1.1_HiC_header.sam

########################
# STEP 1: Extract reads in region
########################

echo "[1/9] Extracting reads in genomic interval ${REGION_START}-${REGION_END}..."
awk -v start="${REGION_START}" -v end="${REGION_END}" '
  ($4 > start && $4 < end) { print }
' "${SAM_BODY}" > positive_region_body.sam

########################
# STEP 2: Add header and convert to BAM
########################

echo "[2/9] Adding header and converting to BAM..."
cat "${SAM_HEADER}" positive_region_body.sam > positive_region.sam

samtools view -bh positive_region.sam -o positive_region.bam

########################
# STEP 3: Filter by SAM flags
########################
# The exact flag combination depends on your design.
# Here we keep reads that are:
#   - properly paired (0x2)
#   - second in pair (0x80)
#   - mate not on reverse strand (exclude 0x20)
# Adapt -f / -F if needed.

echo "[3/9] Filtering BAM by SAM flags (-f 130 -F 32)..."
samtools view       -bh       -f 130       -F 32       -@ "${THREADS}"       positive_region.bam       > positive_mapped.bam

########################
# STEP 4: Sort and index
########################

echo "[4/9] Sorting and indexing..."
samtools sort -@ "${THREADS}" -o positive_mapped_mates.bam positive_mapped.bam
samtools index positive_mapped_mates.bam

READ_COUNT=$(samtools view positive_mapped_mates.bam | wc -l)
echo "[INFO] Reads after filtering: ${READ_COUNT}"

########################
# STEP 5: Normalized subset
########################

echo "[5/9] Creating normalized subset of ${SUBSET_LINES} SAM lines..."
samtools view -h positive_mapped_mates.bam > positive_mapped_mates.sam
head -n "${SUBSET_LINES}" positive_mapped_mates.sam > positive_mapped_mates_N.sam

samtools view -bh positive_mapped_mates_N.sam -o positive_mapped_mates_N.bam

########################
# STEP 6: Convert BAM to FASTA
########################

echo "[6/9] Converting normalized subset to FASTA..."
samtools view positive_mapped_mates_N.bam       | awk '{ printf(">%s\n%s\n", $1, $10) }'       > reads_pos_N.fasta

echo "[6/9] Converting full set to FASTA..."
samtools view positive_mapped_mates.bam       | awk '{ printf(">%s\n%s\n", $1, $10) }'       > reads_pos_all.fasta

########################
# STEP 7: Jellyfish k-mer counting (main)
########################

echo "[7/9] Jellyfish count (k=${JF_KMER_LEN_MAIN})..."
jellyfish count       -m "${JF_KMER_LEN_MAIN}"       -s "${JF_HASH_SIZE}"       reads_pos_all.fasta       -o positive_mapped_mates.jf

echo "[8/9] Histogram and dump..."
jellyfish histo positive_mapped_mates.jf > positive_mapped_mates_histo.txt
jellyfish dump positive_mapped_mates.jf > positive_mapped_mates_dump.jf

########################
# STEP 9: Optional second k
########################

echo "[9/9] Optional Jellyfish run with k=${JF_KMER_LEN_ALT}..."
jellyfish count       -m "${JF_KMER_LEN_ALT}"       -s "${JF_HASH_SIZE}"       reads_pos_all.fasta       -o positive_mapped_mates_k${JF_KMER_LEN_ALT}.jf

echo "Done."
