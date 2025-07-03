#!/usr/bin/env Rscript

# =============================================================================
#
# SCRIPT: Step 4 - Reconstruct Ancestral Linkage Groups (ALGs)
#
# BIOLOGICAL CONTEXT & PURPOSE:
# This script represents the pinnacle of the reconstruction workflow. While
# previous steps estimated the *number* of ancestral chromosomes, this script
# aims to reconstruct their *content and order*. It takes the conserved synteny
# blocks shared among modern species and attempts to stitch them back together
# into the most likely ancestral arrangement.
#
# VERSION: 1.2 (Syntax Fix)
# DATE: 2025-07-03
# AUTHOR: Gemini AI
#
# KEY IMPROVEMENTS in v1.2:
# 1.  Fixed a critical syntax error (missing closing brace) in the main loop.
#
# =============================================================================

# --- 1. Load Libraries ---
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ape)
  library(phytools)
  library(igraph)
  library(progress)
})

# --- 2. Command Line Argument Parsing ---
option_list <- list(
  make_option(c("--synteny_blocks"), type="character", 
              help="Path to the consolidated synteny_blocks_all.tsv file from step2"),
  make_option(c("--ancestral_karyotype"), type="character", 
              help="Path to the ancestral_karyotypes.tsv file from step3"),
  make_option(c("--phylogeny"), type="character", 
              help="Path to the Newick-formatted species tree"),
  make_option(c("--target_ancestor_node"), type="character", 
              help="The label of the ancestral node to reconstruct (e.g., 'Node_68')"),
  make_option(c("--output_dir"), type="character", default=".", 
              help="Output directory [default: current directory]")
)
parser <- OptionParser(option_list=option_list, usage = "%prog [options]")
args <- parse_args(parser)

# --- 3. Utility Functions ---
log_message <- function(msg, level="INFO") {
  cat(sprintf("[%s] [%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, msg))
}

ensure_dir <- function(dir) {
  if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
}

# --- 4. Core Reconstruction Logic ---

#' Build the adjacency graph from synteny blocks of descendant species.
build_adjacency_graph <- function(synteny_blocks, descendant_species) {
  log_message("Building adjacency graph from descendant genomes...")

  adj_data <- synteny_blocks[species1 %in% descendant_species & species2 %in% descendant_species]
  log_message(paste("Found", nrow(adj_data), "synteny blocks within", length(descendant_species), "descendant species."))

  all_adjacencies <- list()
  for (species in descendant_species) {
    s1_blocks <- adj_data[species1 == species, .(chr = chr1, block_id, pos = start1, dir = direction)]
    s2_blocks <- adj_data[species2 == species, .(chr = chr2, block_id, pos = start2, dir = -direction)]
    species_blocks <- rbind(s1_blocks, s2_blocks)
    if(nrow(species_blocks) == 0) next

    setkey(species_blocks, chr, pos)

    chr_adjacencies <- species_blocks[, .(
      block1 = block_id[-.N],
      block2 = block_id[-1],
      dir1 = dir[-.N],
      dir2 = dir[-1]
    ), by = chr]

    all_adjacencies[[species]] <- chr_adjacencies
  }
  all_adjacencies_dt <- rbindlist(all_adjacencies)
  if(nrow(all_adjacencies_dt) == 0) stop("Could not find any block adjacencies in descendant species.")

  all_adjacencies_dt[, extremity1 := ifelse(dir1 == 1, paste0(block1, "_t"), paste0(block1, "_h"))]
  all_adjacencies_dt[, extremity2 := ifelse(dir2 == 1, paste0(block2, "_h"), paste0(block2, "_t"))]

  edge_list <- all_adjacencies_dt[, .N, by = .(extremity1, extremity2)]
  setnames(edge_list, "N", "weight")
  setorder(edge_list, -weight)

  log_message(paste("Graph has", length(unique(c(edge_list$extremity1, edge_list$extremity2))), "nodes and", nrow(edge_list), "edges."))
  return(edge_list)
}

#' Reconstruct ALGs using a greedy algorithm.
reconstruct_algs_greedy <- function(edge_list, all_blocks, num_chromosomes) {
  log_message("Reconstructing ALGs with a greedy algorithm...")
  
  contigs <- as.list(all_blocks)
  names(contigs) <- all_blocks
  
  block_to_contig <- setNames(all_blocks, all_blocks)
  
  contig_ends <- list()
  for (block in all_blocks) {
    contig_ends[[block]] <- c(paste0(block, "_h"), paste0(block, "_t"))
  }

  pb <- progress_bar$new(total = nrow(edge_list), format = "  Merging contigs [:bar] :percent")
  
  for (i in 1:nrow(edge_list)) {
    pb$tick()
    if (length(contigs) <= num_chromosomes) {
      log_message(paste("Target number of", num_chromosomes, "chromosomes reached."))
      break
    }

    edge <- edge_list[i]
    e1 <- edge$extremity1
    e2 <- edge$extremity2

    b1 <- sub("_[ht]$", "", e1)
    b2 <- sub("_[ht]$", "", e2)

    c1_name <- block_to_contig[b1]
    c2_name <- block_to_contig[b2]

    if (c1_name != c2_name && e1 %in% contig_ends[[c1_name]] && e2 %in% contig_ends[[c2_name]]) {
      c1 <- contigs[[c1_name]]
      c2 <- contigs[[c2_name]]

      if (e1 != contig_ends[[c1_name]][2]) c1 <- rev(c1)
      if (e2 != contig_ends[[c2_name]][1]) c2 <- rev(c2)

      new_contig <- c(c1, c2)
      new_contig_name <- new_contig[1]
      contigs[[new_contig_name]] <- new_contig
      
      block_to_contig[new_contig] <- new_contig_name
      contig_ends[[new_contig_name]] <- c(contig_ends[[c1_name]][1], contig_ends[[c2_name]][2])
      
      contigs[[c1_name]] <- NULL
      contigs[[c2_name]] <- NULL
      contig_ends[[c1_name]] <- NULL
      contig_ends[[c2_name]] <- NULL
    }
  } # CORRECTED: Added missing closing brace for the for-loop
  
  log_message(paste("Greedy merge resulted in", length(contigs), "ALGs."))
  return(contigs)
}


# --- 5. Main Execution Block ---
main <- function(args) {
  ensure_dir(args$output_dir)
  log_message("=== Step 4: Ancestral Linkage Group Reconstruction Started ===")

  # --- Load Data ---
  log_message("Loading all required data files...")
  synteny_blocks <- fread(args$synteny_blocks)
  ancestral_karyotype <- fread(args$ancestral_karyotype)
  tree <- read.tree(args$phylogeny)

  # --- Get Target Info ---
  target_karyotype <- ancestral_karyotype[node_label == args$target_ancestor_node]
  if (nrow(target_karyotype) == 0) stop("Target ancestor node label not found in karyotype file.")
  num_chromosomes <- target_karyotype$estimated_count
  target_node_id <- target_karyotype$node
  log_message(paste("Targeting ancestor", args$target_ancestor_node, "(Node ID:", target_node_id, ") with predicted N=", num_chromosomes))

  # Get descendant species for the target node
  descendant_tips_idx <- getDescendants(tree, target_node_id)
  descendant_tips_idx <- descendant_tips_idx[descendant_tips_idx <= length(tree$tip.label)]
  descendant_species <- tree$tip.label[descendant_tips_idx]

  # --- Build Graph ---
  edge_list <- build_adjacency_graph(synteny_blocks, descendant_species)

  # --- Reconstruct ALGs ---
  all_blocks_in_graph <- unique(sub("_[ht]$", "", c(edge_list$extremity1, edge_list$extremity2)))
  reconstructed_algs <- reconstruct_algs_greedy(edge_list, all_blocks_in_graph, num_chromosomes)

  # --- Format and Save Output ---
  log_message("Formatting and saving reconstruction results...")
  output_list <- lapply(1:length(reconstructed_algs), function(i) {
    data.table(ALG_ID = paste0("ALG_", i), block_id = reconstructed_algs[[i]], order = 1:length(reconstructed_algs[[i]]))
  })
  output_dt <- rbindlist(output_list)

  output_tsv <- file.path(args$output_dir, paste0("ancestor_", args$target_ancestor_node, "_algs.tsv"))
  fwrite(output_dt, output_tsv, sep="\t")
  log_message(paste("Reconstructed ALGs saved to:", output_tsv))

  # --- Summary Report ---
  summary_text <- c(
    paste("=== ALG Reconstruction Summary for Ancestor", args$target_ancestor_node, "==="),
    paste("Analysis Date:", Sys.time()), "",
    paste("Predicted Chromosome Number (from Step 3):", num_chromosomes),
    paste("Reconstructed Ancestral Linkage Groups (ALGs):", length(reconstructed_algs)),
    paste("Total Synteny Blocks Used in Reconstruction:", length(all_blocks_in_graph)),
    paste("Average Blocks per ALG:", round(mean(sapply(reconstructed_algs, length)), 2))
  )
  writeLines(summary_text, file.path(args$output_dir, paste0("summary_ancestor_", args$target_ancestor_node, ".txt")))

  log_message("=== Analysis Successfully Completed ===")
}

# --- Script Entry Point ---
if (!interactive()) {
  main(args)
}
