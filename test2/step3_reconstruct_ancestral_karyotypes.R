#!/usr/bin/env Rscript

# =============================================================================
#
# SCRIPT: Ancestral Chromosome Number Reconstruction via Parsimony
#
# BIOLOGICAL CONTEXT & PURPOSE:
# This script reconstructs the probable number of chromosomes (ancestral karyotype)
# for the ancestors in a given phylogenetic tree. It uses a robust, scientifically
# established method based on the principle of parsimony (Occam's Razor).
#
# This script is self-contained and only requires a phylogeny and a file with
# chromosome counts for the tip species.
#
# VERSION: 10.1 (Final & Self-Contained)
# DATE: 2025-07-03
# AUTHOR: Biols9527 (Refactored by Gemini AI)
#
# KEY IMPROVEMENTS in v10.1:
# 1.  **Logical Self-Consistency:** Removed the loading of the `--mappings_dir` data,
#     as it is not used in the pure parsimony algorithm. The script is now
#     leaner and its logic is fully self-consistent.
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
})

# --- 2. Command Line Argument Parsing ---
option_list <- list(
  make_option(c("--karyotypes"), type="character", help="CSV file of extant species chromosome counts (species,chromosome_count)"),
  make_option(c("--phylogeny"), type="character", help="Newick-formatted species tree"),
  make_option(c("--output_dir"), type="character", default=".", help="Output directory [default: current directory]"),
  make_option(c("--label_size"), type="numeric", default=3, help="Font size for labels on the output tree plot [default: %default]"),
  make_option(c("--width_scale_factor"), type="numeric", default=0.25, help="Scaling factor to adjust plot width for long labels [default: %default]")
)
parser <- OptionParser(option_list=option_list, usage = "%prog --karyotypes <path> --phylogeny <path> [options]")
args <- parse_args(parser)

# --- 3. Utility Functions ---
log_message <- function(msg, level="INFO") {
  cat(sprintf("[%s] [%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, msg))
}

ensure_dir <- function(dir) {
  if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
}

# --- 4. Data Loading and Preparation ---
load_and_validate_data <- function(args) {
  log_message("Loading and validating input data...")
  for (file_arg in c("karyotypes", "phylogeny")) {
    if (is.null(args[[file_arg]]) || !file.exists(args[[file_arg]])) {
      stop(paste("Path not found or specified for argument:", file_arg))
    }
  }

  tree <- read.tree(args$phylogeny)
  karyotypes <- fread(args$karyotypes)

  if(!all(c("species", "chromosome_count") %in% names(karyotypes))) stop("Karyotype file needs 'species' and 'chromosome_count' columns.")

  common_species <- intersect(tree$tip.label, karyotypes$species)
  pruned_tree <- keep.tip(tree, common_species)
  log_message(paste("Analysis pruned to", length(common_species), "common species."))

  extant_counts <- setNames(karyotypes$chromosome_count, karyotypes$species)
  extant_counts <- extant_counts[pruned_tree$tip.label]

  return(list(tree = pruned_tree, extant_counts = extant_counts))
}

# --- 5. Ancestral State Reconstruction ---

reconstruct_ancestors_parsimony <- function(tree, extant_counts) {
  log_message("Starting bottom-up reconstruction using weighted parsimony...")
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

    if (length(best_states) > 1) {
        confidence <- 1 / length(best_states)
    } else {
        confidence <- 1.0
    }

    estimates[node_id, `:=`(estimated_count = best_anc_count, confidence = confidence)]
  }
  log_message("Parsimony reconstruction complete.")
  return(estimates)
}

# --- 6. Output and Visualization ---

generate_outputs <- function(final_estimates, tree, args) {
  log_message("Generating enhanced output files...")
  n_tips <- length(tree$tip.label)

  node_map <- data.table(node = 1:(n_tips + tree$Nnode))
  node_map[1:n_tips, node_label := tree$tip.label]
  if (!is.null(tree$node.label) && length(tree$node.label) == tree$Nnode) {
    node_map[(n_tips + 1):nrow(node_map), node_label := tree$node.label]
  } else {
    node_map[(n_tips + 1):nrow(node_map), node_label := paste0("Node_", (n_tips + 1):(n_tips + tree$Nnode))]
  }

  results_dt <- merge(final_estimates, node_map, by = "node")
  results_dt[, node_type := ifelse(node <= n_tips, "tip", "ancestor")]
  setcolorder(results_dt, c("node", "node_label", "node_type"))

  output_tsv <- file.path(args$output_dir, "ancestral_karyotypes.tsv")
  fwrite(results_dt, output_tsv, sep="\t")
  log_message(paste("Detailed results saved to:", output_tsv))

  log_message("Generating chromosome evolution event log...")
  edge_dt <- as.data.table(tree$edge)
  setnames(edge_dt, c("V1", "V2"), c("parent_node_id", "child_node_id"))
  edge_dt <- merge(edge_dt, results_dt[, .(parent_node_id=node, parent_label=node_label, parent_count=estimated_count)], by="parent_node_id")
  edge_dt <- merge(edge_dt, results_dt[, .(child_node_id=node, child_label=node_label, child_count=estimated_count)], by="child_node_id")
  edge_dt[, change := child_count - parent_count]
  edge_dt[, event_type := fcase(change > 0, "Fission", change < 0, "Fusion", default = "Stable")]
  setcolorder(edge_dt, c("parent_label", "child_label"))
  fwrite(edge_dt, file.path(args$output_dir, "chromosome_evolution_log.tsv"), sep="\t")
  log_message(paste("Evolution log saved to: chromosome_evolution_log.tsv"))

  log_message("Creating enhanced annotated phylogeny plot...")
  plot_data <- copy(results_dt)
  plot_data[node_type == "tip", tip_label_text := paste0(' ', node_label, ' (N=', estimated_count, ')')]
  plot_data[node_type == "ancestor", node_label_text := paste0("N=", estimated_count)]

  max_label_width <- max(nchar(plot_data$tip_label_text), na.rm = TRUE)
  plot_width <- 8 + max_label_width * args$width_scale_factor

  p <- ggtree(tree) %<+% plot_data +
       geom_tiplab(aes(label=tip_label_text), align=TRUE, size=args$label_size) +
       geom_label(aes(label = ifelse(node_type == 'ancestor', node_label_text, NA), fill = confidence), 
                  size=args$label_size, color="white", fontface="bold", label.padding = unit(0.3, "lines")) +
       scale_fill_gradient(low="#D55E00", high="#009E73", name="Confidence", limits=c(0,1)) +
       labs(title="Ancestral Chromosome Number Reconstruction",
            subtitle=paste("Root (Node", tree$edge[1,1], ") estimated count:", plot_data[node == tree$edge[1,1], estimated_count])) +
       theme_tree2()

  output_pdf <- file.path(args$output_dir, "ancestral_karyotype_tree.pdf")
  ggsave(output_pdf, plot=p, width=plot_width, height=0.4 * n_tips, limitsize=FALSE)
  log_message(paste("Enhanced visualization saved to:", output_pdf))

  report_text <- c(
    "=== Ancestral Karyotype Reconstruction Summary ===",
    paste("Analysis Date:", Sys.time()), "",
    paste("Method: Bottom-up Parsimony"),
    paste("Analyzed Species:", n_tips),
    paste("Reconstructed Ancestral Nodes:", tree$Nnode), "",
    "--- Key Findings ---",
    paste("Estimated Chromosome Number at Root (Node", tree$edge[1,1], "):", results_dt[node == tree$edge[1,1], estimated_count]),
    paste("Average Confidence for Ancestral Nodes:", sprintf("%.3f", mean(results_dt[node_type=='ancestor', confidence]))), ""
  )
  writeLines(report_text, file.path(args$output_dir, "summary_report.txt"))
  log_message("Summary report created.")
}

# --- 7. Main Execution Block ---
main <- function(args) {
  ensure_dir(args$output_dir)
  log_message("=== Ancestral Karyotype Reconstruction (v10.1) Started ===")

  data_list <- load_and_validate_data(args)
  final_results <- reconstruct_ancestors_parsimony(data_list$tree, data_list$extant_counts)
  generate_outputs(final_results, data_list$tree, args)

  log_message("=== Analysis Successfully Completed ===")
}

# --- Script Entry Point ---
if (!interactive()) {
  main(args)
}