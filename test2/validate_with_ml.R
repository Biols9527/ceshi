#!/usr/bin/env Rscript

# =============================================================================
#
# SCRIPT: Methodological Cross-Validation for Ancestral Karyotype Reconstruction
#
# PURPOSE:
# This script provides an independent validation of ancestral chromosome number
# estimates by using a different algorithmic approach: Maximum Likelihood (ML).
# It compares the results from our parsimony-based method with ML estimates
# to assess the robustness of the conclusions.
#
# VERSION: 1.2 (Robust Data Merging Fix)
# DATE: 2025-07-03
# AUTHOR: Gemini AI
#
# KEY IMPROVEMENTS in v1.2:
# 1.  Fixed a `cor()` error by ensuring proper alignment of node labels between
#     parsimony and ML results before merging.
# 2.  Switched to an inner join (`all=FALSE`) to guarantee that only nodes with
#     results from BOTH methods are compared.
#
# =============================================================================

# --- 1. Load Libraries ---
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(phytools)
  library(ggplot2)
  library(ggrepel) # For better plot labels
})

# --- 2. Command Line Argument Parsing ---
option_list <- list(
  make_option(c("--parsimony_results"), type="character",
              help="Path to the ancestral_karyotypes.tsv file from the parsimony script (required for comparison)"),
  make_option(c("--karyotypes"), type="character", 
              help="Path to the CSV file of extant species chromosome counts (required)"),
  make_option(c("--phylogeny"), type="character", 
              help="Path to the Newick-formatted species tree (must be the same as used in the main pipeline)"),
  make_option(c("--output_dir"), type="character", default=".",
              help="Directory for the validation output [default: current directory]")
)
parser <- OptionParser(option_list=option_list, usage = "%prog --parsimony_results <path> --karyotypes <path> --phylogeny <path> [options]")
args <- parse_args(parser)

# --- 3. Utility Functions ---
log_message <- function(msg, level="INFO") {
  cat(sprintf("[%s] [%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, msg))
}

ensure_dir <- function(dir) {
  if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
}

# --- 4. Main Analysis Block ---
main <- function(args) {
  ensure_dir(args$output_dir)
  log_message("=== Methodological Cross-Validation (Parsimony vs. ML) (v1.2) Started ===")

  # --- Data Loading ---
  log_message("Loading parsimony results, karyotypes, and phylogeny...")
  if (is.null(args$parsimony_results) || !file.exists(args$parsimony_results)) stop("Parsimony results file not found.")
  if (is.null(args$karyotypes) || !file.exists(args$karyotypes)) stop("Karyotypes file not found.")
  if (is.null(args$phylogeny) || !file.exists(args$phylogeny)) stop("Phylogeny file not found.")

  parsimony_dt <- fread(args$parsimony_results)
  karyotypes <- fread(args$karyotypes)
  tree <- read.tree(args$phylogeny)

  # --- Data Preparation ---
  common_species <- intersect(tree$tip.label, karyotypes$species)
  if(length(common_species) < 3) stop("Fewer than 3 common species between tree and data.")
  tree <- keep.tip(tree, common_species)
  extant_counts_vec <- setNames(karyotypes[species %in% common_species, chromosome_count], karyotypes[species %in% common_species, species])
  extant_counts_vec <- extant_counts_vec[tree$tip.label]

  # --- ML Reconstruction ---
  log_message("Performing Maximum Likelihood ancestral state reconstruction...")
  ml_results <- anc.ML(tree, extant_counts_vec, model = "BM")

  # --- Comparison ---
  log_message("Comparing Parsimony and ML results...")
  parsimony_ancestors <- parsimony_dt[node_type == "ancestor"]
  
  # CORRECTED: Create a robust mapping from node ID to node label from the parsimony results
  # This ensures the labels are consistent before merging.
  n_tips <- length(tree$tip.label)
  node_map <- parsimony_ancestors[, .(node, node_label)]
  ml_ace_dt <- data.table(node = (n_tips + 1):(n_tips + tree$Nnode), ml_estimate = ml_results$ace)
  ml_dt <- merge(ml_ace_dt, node_map, by = "node")

  # Use an INNER JOIN (all=FALSE) to ensure we only compare nodes present in both results.
  comparison_dt <- merge(parsimony_ancestors, ml_dt, by = c("node", "node_label"), all = FALSE)
  comparison_dt[, `:=` (parsimony_estimate = estimated_count, ml_estimate = round(ml_estimate, 2))]
  
  if(nrow(comparison_dt) == 0) {
      log_message("Could not find any common ancestral nodes between Parsimony and ML results. Check node labels.", "ERROR")
      stop("Halting due to no common data for comparison.")
  }

  correlation <- cor(comparison_dt$parsimony_estimate, comparison_dt$ml_estimate, use="complete.obs")
  log_message(paste("Pearson correlation between methods:", round(correlation, 4)))

  # --- Output Generation ---
  output_tsv <- file.path(args$output_dir, "comparison_parsimony_vs_ml.tsv")
  fwrite(comparison_dt[, .(node_label, parsimony_estimate, ml_estimate, confidence)], output_tsv, sep="\t")
  log_message(paste("Comparison table saved to:", output_tsv))

  plot <- ggplot(comparison_dt, aes(x = parsimony_estimate, y = ml_estimate)) +
    geom_point(aes(color = confidence), size = 4, alpha = 0.8) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    geom_text_repel(aes(label = node_label), size = 2.5, max.overlaps = 15) +
    scale_color_gradient(low="#D55E00", high="#009E73", name="Parsimony Confidence") +
    labs(
      title = "Parsimony vs. Maximum Likelihood Reconstruction",
      subtitle = paste("Pearson Correlation:", round(correlation, 3)),
      x = "Ancestral Chromosome # (Parsimony Estimate)",
      y = "Ancestral Chromosome # (ML Estimate)"
    ) +
    theme_bw()

  output_pdf <- file.path(args$output_dir, "correlation_plot_parsimony_vs_ml.pdf")
  ggsave(output_pdf, plot = plot, width = 8, height = 7)
  log_message(paste("Correlation plot saved to:", output_pdf))

  log_message("=== Validation with ML complete. ===")
}

# --- Script Entry Point ---
if (!interactive()) {
  main(args)
}