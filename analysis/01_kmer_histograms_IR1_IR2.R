###############################################################################
# 01_kmer_histograms_IR1_IR2.R
#
# Plot k-mer histograms for:
#   - unmapped mates in the interesting region
#   - unmapped mates in random regions (negative controls)
#   - mapped mates in the interesting region (positive control)
#
# For each interesting region we:
#   1) read Jellyfish histograms for a single interesting region
#   2) read multiple negative-control histograms
#   3) compute median and ±1 SD across negatives
#   4) overlay curves on a log-scaled plot
###############################################################################

library(dplyr)
library(purrr)

# Helper: read a Jellyfish histogram (two columns: frequency, count)
read_histo <- function(path) {
  read.table(path, col.names = c("frequency", "count"))
}

# Helper: compute median and ±SD envelope across negative controls
summarise_negative_controls <- function(neg_list) {
  # neg_list: named list of data.frames with columns frequency, count
  neg_df <- bind_rows(neg_list, .id = "replicate")

  neg_df %>%
    group_by(frequency) %>%
    summarise(
      median_count = median(count),
      sd_count = sd(count),
      .groups = "drop"
    ) %>%
    mutate(
      lower = median_count - sd_count,
      upper = median_count + sd_count
    )
}

# Helper: base plotting function
plot_kmer_comparison <- function(
  interesting,
  pos_control,
  neg_summary,
  main_title = "",
  xlab = "k-mer count (frequency)",
  ylab = "Number of distinct k-mers"
) {
  # Base plot: interesting region (unmapped mates)
  plot(
    interesting$frequency,
    interesting$count,
    type = "l",
    log = "y",
    xlab = xlab,
    ylab = ylab,
    main = main_title,
    col = "black",
    lwd = 2
  )

  # Negative controls: median and ±SD
  lines(
    neg_summary$frequency,
    neg_summary$median_count,
    col = "blue",
    lty = 1
  )

  lines(
    neg_summary$frequency,
    pmax(neg_summary$lower, 1e-6),  # protect log scale
    col = "blue",
    lty = 2
  )

  lines(
    neg_summary$frequency,
    pmax(neg_summary$upper, 1e-6),
    col = "blue",
    lty = 2
  )

  # Positive control: mapped mates in interesting region
  lines(
    pos_control$frequency,
    pos_control$count,
    col = "red",
    lty = 2
  )

  legend(
    "topright",
    legend = c(
      "Unmapped mates in interesting region",
      "Unmapped mates in random regions (median ± SD)",
      "Mapped mates in interesting region"
    ),
    col = c("black", "blue", "red"),
    lty = c(1, 1, 2),
    cex = 0.8
  )
}

###############################################################################
# Example 1: Interesting region IR_1 (k = 100)
###############################################################################

# Files for IR_1 (adapt paths / names as needed)
interesting_IR1_file <- "PM_inv.mate_unmapped.IR_1.100.histo"
positive_IR1_file    <- "PM_inv.mate_mapped.IR_1.N.100.histo"

negative_IR1_files <- c(
  "PM_inv.IR_1.neg_1.100.histo",
  "PM_inv.IR_1.neg_2.100.histo",
  "PM_inv.IR_1.neg_3.100.histo",
  "PM_inv.IR_1.neg_4.100.histo",
  "PM_inv.IR_1.neg_5.100.histo",
  "PM_inv.IR_1.neg_6.100.histo",
  "PM_inv.IR_1.neg_7.100.histo",
  "PM_inv.IR_1.neg_8.100.histo",
  "PM_inv.IR_1.neg_9.100.histo",
  "PM_inv.IR_1.neg_10.100.histo"
)

# Read data
interesting_IR1    <- read_histo(interesting_IR1_file)
positive_control_1 <- read_histo(positive_IR1_file)

negative_IR1_list <- setNames(
  lapply(negative_IR1_files, read_histo),
  nm = basename(negative_IR1_files)
)

# Summarise negatives (median ± SD)
neg_summary_IR1 <- summarise_negative_controls(negative_IR1_list)

# Plot
plot_kmer_comparison(
  interesting = interesting_IR1,
  pos_control = positive_control_1,
  neg_summary = neg_summary_IR1,
  main_title = "IR_1: 100-mers, interesting vs random regions"
)

###############################################################################
# Example 2: Interesting region IR_2 (k = 150)
###############################################################################

interesting_IR2_file <- "PM_WT.mate_unmapped.IR_2.150.txt"
positive_IR2_file    <- "PM_WT.mate_mapped.IR_2.N.150.txt"

negative_IR2_files <- c(
  "PM_WT.IR_2.neg_1.noH.150.txt",
  "PM_WT.IR_2.neg_2.noH.150.txt",
  "PM_WT.IR_2.neg_3.noH.150.txt",
  "PM_WT.IR_2.neg_4.noH.150.txt",
  "PM_WT.IR_2.neg_5.noH.150.txt",
  "PM_WT.IR_2.neg_6.noH.150.txt",
  "PM_WT.IR_2.neg_7.noH.150.txt"
)

interesting_IR2    <- read_histo(interesting_IR2_file)
positive_control_2 <- read_histo(positive_IR2_file)

negative_IR2_list <- setNames(
  lapply(negative_IR2_files, read_histo),
  nm = basename(negative_IR2_files)
)

neg_summary_IR2 <- summarise_negative_controls(negative_IR2_list)

plot_kmer_comparison(
  interesting = interesting_IR2,
  pos_control = positive_control_2,
  neg_summary = neg_summary_IR2,
  main_title = "IR_2: 150-mers, interesting vs random regions"
)
