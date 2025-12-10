# Hi-C k-mer enrichment around interesting regions

This repository contains example scripts used to:

1. Extract read mates from Hi-C BAM files in a given genomic interval  
   (both mapped and unmapped mates).
2. Scan a scaffold with sliding windows and extract unmapped mates whose
   partners map inside each window (negative controls / random regions).
3. Count and inspect k-mer frequencies with Jellyfish.
4. Plot k-mer histograms for interesting vs random regions in R.

The code is meant as a **reproducible example** rather than a general-purpose
pipeline: you will need to adapt file paths, scaffold names, and coordinates
to your own data. In practice, we often rely on visual assesment of data using 
IGV, JuicerBox, UGENE, and simple text editors, and GNU core utilities 
(`awk`, `sort`, `uniq`, `head`, etc.) for data manipulation.

## Dependencies

- `samtools`
- `jellyfish`
- GNU core utilities (`awk`, `sort`, `uniq`, `head`, etc.)
- R with packages: `dplyr`, `purrr`

## Scripts

### `scripts/00_region_mapped_mates_pipeline.sh`

Example pipeline for a single genomic interval:

- selects reads in a given coordinate range from a SAM file
- restores SAM header and converts to BAM
- filters paired reads by SAM flags
- sorts and indexes BAM
- (optionally) subsamples reads
- converts to FASTA
- runs Jellyfish k-mer counting and produces histogram and dump

This script was used to generate a **positive control**:
mapped mates in the interesting region.

### `scripts/01_window_unmapped_mates.sh`

Scans a scaffold with sliding windows and, for each window:

1. Collects read names in a pre-filtered BAM (`mate_unmapped`).
2. From the full BAM, selects reads with these names whose **mates are unmapped**.
3. Exports them as FASTA.

This script was used both for the **target region** and for multiple
**random regions** (negative controls).

### `scripts/02_kmer_count_and_dump.sh`

For every FASTA file in a directory:

- runs `jellyfish count` (k-mer length configurable, e.g. k=80)
- saves a k-mer histogram
- dumps all k-mers with counts
- dumps only high-copy k-mers with count ≥ threshold (e.g. ≥80)

Used to:
- determine that the interesting region differs from negative controls
  around a specific k-mer frequency (~60 in our case),
- inspect sequences repeated ≥80 times.

### `analysis/01_kmer_histograms_IR1_IR2.R`

R script that:

- reads Jellyfish histograms for an interesting region,
  positive control (mapped mates), and multiple negative controls
  (random regions),
- computes median and ±1 SD envelopes across negative controls,
- plots log-scaled k-mer histograms:
  - black — unmapped mates in the interesting region,
  - blue — median and ±SD of random regions,
  - red — mapped mates in the interesting region.

Two example blocks are provided for two interesting regions (IR_1 and IR_2).

## Manual step: UGENE

To investigate highly repeated k-mers (e.g. ≥80 copies), we:

1. Dump high-copy k-mers from Jellyfish:
   ```bash
   jellyfish dump -L 80 int_reg_149.jf > int_reg_149.highcopy.fa
