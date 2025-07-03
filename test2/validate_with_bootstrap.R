#!/usr/bin/env Rscript

# =============================================================================
#
# SCRIPT: Statistical Robustness Test via Bootstrapping
#
# PURPOSE:
# This script assesses the statistical robustness of the ancestral chromosome
# number reconstruction by performing a non-parametric bootstrap analysis.
# It addresses the question: How much would our results change if our input
# data were slightly different?
#
# VERSION: 1.2 (Final Merge Fix)
# DATE: 2025-07-03
# AUTHOR: Gemini AI
#
# KEY IMPROVEMENTS in v1.2:
# 1.  Fixed a data.table merge error by ensuring join keys are present in
#     all intermediate tables during bootstrap summary calculation.
#
# =============================================================================

# --- 1. Load Libraries ---
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ape)
  library(phytools)
  library(ggplot2)
  library(ggtree)
  library(progress)
})

# --- 2. Command Line Argument Parsing ---
option_list <- list(
  make_option(c("--mappings_dir"), type="character", 
              help="Path to the DIRECTORY containing mapping files from step1 (required)"),
  make_option(c("--karyotypes"), type="character", help="CSV file of extant species chromosome counts (species,chromosome_count)"),
  make_option(c("--phylogeny"), type="character", help="Newick-formatted species tree"),
  make_option(c("--output_dir"), type="character", default=".", help="Output directory [default: current directory]"),
  make_option(c("--n_bootstraps"), type="integer", default=100, help="Number of bootstrap replicates to perform [default: %default]")
)
parser <- OptionParser(option_list=option_list, usage = "%prog --mappings_dir <path> --karyotypes <path> --phylogeny <path> [options]")
args <- parse_args(parser)

# --- 3. Utility & Core Reconstruction Functions (Self-contained) ---
log_message <- function(msg, level="INFO") {
  cat(sprintf("[%s] [%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, msg))
}

ensure_dir <- function(dir) {
  if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
}

reconstruct_ancestors_parsimony <- function(tree, extant_counts) {
  n_tips <- length(tree$tip.label)
  n_nodes <- tree$Nnode
  estimates <- data.table(node = 1:(n_tips + n_nodes), estimated_count = NA_integer_, confidence = 0.0)
  estimates[1:n_tips, `:=`(estimated_count = extant_counts, confidence = 1.0)]
  postorder_nodes <- unique(tree$edge[order(tree$edge[,1], decreasing = TRUE), 1])
  postorder_nodes <- postorder_nodes[postorder_nodes > n_tips]
  for (node_id in postorder_nodes) {
    children <- tree$edge[tree$edge[,1] == node_id, 2]
    child_counts <- estimates[children, estimated_count]
    child_counts <- na.omit(child_counts)
    if (length(child_counts) == 0) next
    search_range <- round(min(child_counts) * 0.7) : round(max(child_counts) * 1.3)
    search_range <- unique(pmax(1, search_range))
    costs <- sapply(search_range, function(anc_count) sum(abs(anc_count - child_counts)))
    min_cost <- min(costs)
    best_states <- search_range[costs == min_cost]
    best_anc_count <- round(mean(best_states))
    confidence <- 1 / length(best_states)
    estimates[node_id, `:=`(estimated_count = best_anc_count, confidence = confidence)]
  }
  return(estimates)
}

# --- 4. Main Analysis Block ---
main <- function(args) {
  ensure_dir(args$output_dir)
  log_message("=== Statistical Robustness Test via Bootstrapping (v1.2) Started ===")

  # --- Load Original Data ---
  log_message("Loading original data...")
  tree <- read.tree(args$phylogeny)
  karyotypes <- fread(args$karyotypes)
  
  mapping_files <- list.files(args$mappings_dir, pattern = "\\.tsv$", full.names = TRUE)
  if (length(mapping_files) == 0) stop("No .tsv mapping files found in the specified directory.")
  mappings_list <- lapply(mapping_files, fread)
  mappings <- rbindlist(mappings_list, use.names = TRUE, fill = TRUE)

  common_species <- intersect(tree$tip.label, karyotypes$species)
  tree <- keep.tip(tree, common_species)
  extant_counts <- setNames(karyotypes[species %in% common_species, chromosome_count], karyotypes[species %in% common_species, species])
  extant_counts <- extant_counts[tree$tip.label]

  # --- Run Original Analysis to get the point estimate ---
  log_message("Running analysis on original data to get point estimates...")
  original_results <- reconstruct_ancestors_parsimony(tree, extant_counts)
  original_ancestors <- original_results[node > length(tree$tip.label)]

  # --- Bootstrap Loop ---
  log_message(paste("Starting bootstrap analysis with", args$n_bootstraps, "replicates..."))
  bootstrap_results_list <- list()
  pb <- progress_bar$new(total = args$n_bootstraps, format = "  Bootstrapping [:bar] :percent eta: :eta")

  for (i in 1:args$n_bootstraps) {
    boot_extant_counts <- sample(extant_counts, replace = TRUE)
    names(boot_extant_counts) <- names(extant_counts)

    boot_run_results <- reconstruct_ancestors_parsimony(tree, boot_extant_counts)
    bootstrap_results_list[[i]] <- boot_run_results[node > length(tree$tip.label), . (node, boot_estimate = estimated_count)]
    pb$tick()
  }

  # --- Summarize Bootstrap Results ---
  log_message("Summarizing bootstrap results...")
  all_boot_results <- rbindlist(bootstrap_results_list, idcol = "replicate")

  summary_dt <- all_boot_results[, .(
    median_estimate = round(median(boot_estimate)),
    conf_interval_low = round(quantile(boot_estimate, 0.025)),
    conf_interval_high = round(quantile(boot_estimate, 0.975))
  ), by = node]

  summary_dt <- merge(summary_dt, original_ancestors[, .(node, original_estimate = estimated_count)], by = "node")

  # CORRECTED LOGIC for calculating support
  keys_to_check <- summary_dt[, .(node, original_estimate)]
  support_counts <- all_boot_results[keys_to_check, on = .(node, boot_estimate = original_estimate), .N, by = .EACHI]
  
  summary_dt <- merge(summary_dt, support_counts[, .(node, original_estimate=boot_estimate, bootstrap_support_count=N)], by = c("node", "original_estimate"), all.x = TRUE)
  summary_dt[is.na(bootstrap_support_count), bootstrap_support_count := 0]
  summary_dt[, bootstrap_support_percent := round(100 * bootstrap_support_count / args$n_bootstraps)]

  # --- Output Generation ---
  log_message("Generating output files...")
  output_tsv <- file.path(args$output_dir, "bootstrap_summary.tsv")
  fwrite(summary_dt, output_tsv, sep="\t")
  log_message(paste("Bootstrap summary table saved to:", output_tsv))

  node_map <- data.table(node = original_ancestors$node, node_label = paste0("Node_", original_ancestors$node))
  plot_data <- merge(summary_dt, node_map, by="node")
  plot_data[, node_label_text := paste0(original_estimate, " (", bootstrap_support_percent, "%)")]
  
  p <- ggtree(tree) %<+% plot_data +
       geom_tiplab(align=TRUE) +
       geom_label(aes(label=node_label_text, fill=bootstrap_support_percent), size=3, color="white", fontface="bold") +
       scale_fill_gradient(low="#E69F00", high="#009E73", name="Bootstrap Support (%)", limits=c(0,100)) +
       labs(title="Ancestral Karyotype Reconstruction with Bootstrap Support") +
       theme_tree2()

  output_pdf <- file.path(args$output_dir, "bootstrap_tree.pdf")
  ggsave(output_pdf, plot=p, width=10, height=0.4 * length(tree$tip.label), limitsize=FALSE)
  log_message(paste("Bootstrap visualization saved to:", output_pdf))

  log_message("=== Bootstrap validation complete. ===")
}

# --- Script Entry Point ---
if (!interactive()) {
  main(args)
}