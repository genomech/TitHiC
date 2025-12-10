#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# 02_kmer_count_and_dump.sh
#
# For every FASTA file in a directory:
#   - run Jellyfish k-mer counting (k=N)
#   - save histogram
#   - dump all k-mers with counts
#   - dump only high-copy k-mers with count >= MIN_COUNT
#
# Used for:
#   - comparing interesting vs random regions
#   - inspecting sequences repeated >= 80 times
###############################################################################

########################
# CONFIGURATION
########################

INPUT_DIR=/home/svetlana/clu_scratch/202305101343data/Parus_major_inv_HiC/new_map/kmers_regions
OUTPUT_DIR=kmer_counts

# k-mer length
K=80

# Jellyfish hash size (adjust to coverage)
HASH_SIZE=100M

# Minimum count threshold for high-copy k-mers
MIN_COUNT=80

mkdir -p "${OUTPUT_DIR}"

########################
# MAIN LOOP
########################

for file in "${INPUT_DIR}"/*.fa; do
  [ -e "$file" ] || continue

  base=$(basename "$file" .fa)

  jf_out="${OUTPUT_DIR}/${base}.mates.k${K}.jf"
  hist_out="${OUTPUT_DIR}/${base}.mates.k${K}.histo"
  dump_all_out="${OUTPUT_DIR}/${base}.mates.k${K}.all.fa"
  dump_high_out="${OUTPUT_DIR}/${base}.mates.k${K}.ge${MIN_COUNT}.fa"

  echo "Processing ${file} -> ${jf_out}"

  # k-mer counting
  jellyfish count         -m "${K}"         -s "${HASH_SIZE}"         "$file"         -o "${jf_out}"

  # histogram of k-mer frequencies
  jellyfish histo "${jf_out}" > "${hist_out}"

  # dump all k-mers with counts
  jellyfish dump -c "${jf_out}" > "${dump_all_out}"

  # dump only high-copy k-mers (count >= MIN_COUNT)
  jellyfish dump -L "${MIN_COUNT}" "${jf_out}" > "${dump_high_out}"
done

echo "Done."
