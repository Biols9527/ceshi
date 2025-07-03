#!/usr/bin/env Rscript

# =============================================================================
#
# SCRIPT: Enhanced Chromosome Mapping Type Analysis
#
# BIOLOGICAL CONTEXT & PURPOSE:
# This script serves as a crucial step in comparative genomics, specifically for
# understanding the macro-level evolution of chromosome structure. After an
# initial synteny analysis has identified conserved blocks of genes between
# different species, this script takes the summarized results (a "trace file")
# to classify the relationship between entire chromosomes.
#
# The core biological questions it answers are:
#   - Is chromosome A in Species 1 a direct, one-to-one equivalent of a
#     chromosome in Species 2? (Indicates chromosomal stability).
#   - Does chromosome A in Species 1 correspond to multiple chromosomes in
#     Species 2? (Suggests a chromosome fission event in the lineage of Species 2).
#   - Do multiple chromosomes in Species 1 all point to a single chromosome in
#     Species 2? (Suggests chromosome fusion events in the lineage of Species 2).
#
# By quantifying these 1:1, 1:n, n:1, and n:m relationships, we can build a
# detailed map of the chromosomal rearrangements (fissions, fusions) that
# have shaped the genomes of the species under study.
#
# VERSION: 2.1 (Corrected Plotting Logic)
# DATE: 2025-07-03
# AUTHOR: ShrivatsaMehan (Refactored by Gemini AI)
#
# =============================================================================

# --- 1. Load Libraries ---
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(RColorBrewer))

# --- 2. Command Line Argument Parsing ---
option_list <- list(
  make_option(c("--trace_file"), type="character", default=NULL,
              help="Path to the chromosome trace file (required)"),
  make_option(c("--output_file"), type="character", default=NULL,
              help="Base name for output files"),
  make_option(c("--min_count"), type="numeric", default=10,
              help="Minimum synteny block count to consider a link [default: %default]"),
  make_option(c("--species_focus"), type="character", default=NULL,
              help="Optional: Analyze only mappings involving this species"),
  make_option(c("--disable_plot_filtering"), action="store_true", default=FALSE,
              help="Disable the filtering of top 500 links for visualization")
)

parser <- OptionParser(option_list=option_list, usage = "%prog --trace_file <path> [options]")
args <- parse_args(parser)

# --- 3. Utility Functions ---
log_message <- function(msg, level="INFO") {
  cat(sprintf("[%s] [%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, msg))
}

# --- 4. Core Logic: Bidirectional Mapping Analysis ---
determine_bidirectional_mapping_type <- function(trace_data, bidir_data) {
  log_message("Analyzing bidirectional mapping patterns...")

  species_pairs <- unique(trace_data[, .(focal_species, other_species)])
  results <- list()

  for (i in 1:nrow(species_pairs)) {
    sp1 <- species_pairs$focal_species[i]
    sp2 <- species_pairs$other_species[i]

    chr_pairs <- trace_data[focal_species == sp1 & other_species == sp2,
                            .(chr1 = focal_chromosome, chr2 = other_chromosome)]

    if (nrow(chr_pairs) == 0) next

    setkey(bidir_data, focal_species, focal_chromosome, other_species)

    for (j in 1:nrow(chr_pairs)) {
      chr1 <- chr_pairs$chr1[j]
      chr2 <- chr_pairs$chr2[j]

      forward_targets <- bidir_data[.(sp1, chr1, sp2), unique(other_chromosome)]
      forward_count <- length(forward_targets)

      reverse_targets <- bidir_data[.(sp2, chr2, sp1), unique(other_chromosome)]
      reverse_count <- length(reverse_targets)

      if (forward_count == 1 && reverse_count == 1) {
        mapping_type <- "1:1"
      } else if (forward_count > 1 && reverse_count == 1) {
        mapping_type <- "1:n (Fission)"
      } else if (forward_count == 1 && reverse_count > 1) {
        mapping_type <- "n:1 (Fusion)"
      } else {
        mapping_type <- "n:m (Complex)"
      }

      pair_count <- trace_data[focal_species == sp1 & focal_chromosome == chr1 &
                               other_species == sp2 & other_chromosome == chr2, count]

      results[[length(results) + 1]] <- list(
        species_A = sp1,
        chromosome_A = chr1,
        species_B = sp2,
        chromosome_B = chr2,
        mapping_type_A_to_B = paste0("1:", forward_count),
        mapping_type_B_to_A = paste0("1:", reverse_count),
        bidirectional_mapping_type = mapping_type,
        count = if(length(pair_count) > 0) pair_count[1] else NA
      )
    }
  }

  if (length(results) > 0) {
    return(rbindlist(results))
  } else {
    log_message("No bidirectional mappings could be calculated.")
    return(NULL)
  }
}

# --- 5. Main Execution Block ---
main <- function(args) {
  if (is.null(args$trace_file)) {
    print_help(parser)
    stop("Required argument --trace_file is missing.", call. = FALSE)
  }
  if (is.null(args$output_file)) {
    args$output_file <- gsub("\\.tsv$|\\.csv$", "", basename(args$trace_file))
    args$output_file <- paste0(args$output_file, "_bidirectional_map")
  }

  log_message(paste("Loading chromosome trace data from:", args$trace_file))
  trace_data <- tryCatch({
    fread(args$trace_file)
  }, error = function(e) {
    stop(paste("Fatal Error: Could not read trace file:", e$message))
  })
  log_message(paste("Loaded", nrow(trace_data), "records."))

  required_cols <- c("focal_species", "focal_chromosome", "other_species", "other_chromosome", "count")
  if (!all(required_cols %in% colnames(trace_data))) {
    stop(paste("Input file is missing required columns:", paste(setdiff(required_cols, colnames(trace_data)), collapse=", ")))
  }

  original_rows <- nrow(trace_data)
  trace_data <- trace_data[count >= args$min_count]
  log_message(paste("Filtered", original_rows - nrow(trace_data), "records with count <", args$min_count, ". Remaining:", nrow(trace_data)))

  if (!is.null(args$species_focus)) {
    trace_data <- trace_data[focal_species == args$species_focus | other_species == args$species_focus]
    log_message(paste("Filtered to", nrow(trace_data), "records involving species:", args$species_focus))
  }
  if (nrow(trace_data) == 0) {
      log_message("No data remains after filtering. Halting analysis.")
      return()
  }

  log_message("Creating bidirectional representation for analysis...")
  reverse_dt <- copy(trace_data)
  setnames(reverse_dt,
           old = c("focal_species", "focal_chromosome", "other_species", "other_chromosome"),
           new = c("other_species", "other_chromosome", "focal_species", "focal_chromosome"))
  bidir_data <- rbindlist(list(trace_data, reverse_dt), use.names = TRUE)
  log_message(paste("Bidirectional data prepared with", nrow(bidir_data), "entries."))

  bidirectional_mappings <- determine_bidirectional_mapping_type(trace_data, bidir_data)

  if (is.null(bidirectional_mappings) || nrow(bidirectional_mappings) == 0) {
    log_message("Analysis did not yield any bidirectional mappings. Halting.")
    return()
  }

  log_message(paste("Generated", nrow(bidirectional_mappings), "bidirectional mapping records."))

  output_tsv <- paste0(args$output_file, ".tsv")
  fwrite(bidirectional_mappings, output_tsv, sep="\t")
  log_message(paste("Saved detailed results to:", output_tsv))

  output_json <- paste0(args$output_file, ".json")
  tryCatch({
    jsonlite::write_json(bidirectional_mappings, output_json, pretty = TRUE)
    log_message(paste("Saved JSON results to:", output_json))
  }, error = function(e) {
    log_message(paste("Warning: Could not write JSON output.", e$message))
  })

  output_pdf <- paste0(args$output_file, "_network.pdf")
  log_message("Generating chromosome relationship network visualization...")
  tryCatch({
    vis_data_unique <- bidirectional_mappings[!duplicated(pmin(paste(species_A, chromosome_A), paste(species_B, chromosome_B)))]
    plot_title <- "Chromosome Relationship Network"

    if (!args$disable_plot_filtering && nrow(vis_data_unique) > 500) {
      log_message("-> WARNING: Too many relationships to plot clearly (>500). Visualizing top 500 by count.", level="WARN")
      log_message("-> To see the full, unfiltered graph, re-run with the --disable_plot_filtering flag.")
      viz_data <- vis_data_unique[order(-count)][1:500]
      plot_title <- paste(plot_title, "(Top 500 Links by Count)")
    } else {
      viz_data <- vis_data_unique
    }

    edges <- data.table(
      from = paste(viz_data$species_A, viz_data$chromosome_A, sep=":"),
      to = paste(viz_data$species_B, viz_data$chromosome_B, sep=":"),
      mapping_type = viz_data$bidirectional_mapping_type,
      count = viz_data$count
    )
    g <- graph_from_data_frame(edges, directed = FALSE)

    V(g)$species <- sapply(strsplit(V(g)$name, ":"), `[`, 1)
    V(g)$chromosome <- sapply(strsplit(V(g)$name, ":"), `[`, 2)
    V(g)$degree <- degree(g)

    species_list <- unique(V(g)$species)
    palette <- colorRampPalette(brewer.pal(min(length(species_list), 9), "Set1"))(length(species_list))
    species_colors <- setNames(palette, species_list)
    V(g)$color <- species_colors[V(g)$species]

    V(g)$size <- 3 + sqrt(V(g)$degree) * 2

    mapping_colors <- c(
        "1:1" = "#377EB8", 
        "1:n (Fission)" = "#4DAF4A", 
        "n:1 (Fusion)" = "#FF7F00", 
        "n:m (Complex)" = "#E41A1C"
    )
    E(g)$color <- mapping_colors[E(g)$mapping_type]
    E(g)$width <- 1 + 4 * (E(g)$count - min(E(g)$count, na.rm=T)) / (max(E(g)$count, na.rm=T) - min(E(g)$count, na.rm=T))

    pdf(output_pdf, width = 14, height = 11)
    plot(g,
         vertex.label = V(g)$chromosome,
         vertex.label.cex = 0.8,
         vertex.label.color = "black",
         layout = layout_with_fr(g),
         main = plot_title)
    legend("bottomright", legend = names(mapping_colors), col = mapping_colors, lty = 1, lwd = 2, title = "Mapping Type", bty = "n")
    legend("topright", legend = names(species_colors), fill = species_colors, title = "Species", bty = "n")
    dev.off()
    log_message(paste("Saved network visualization to:", output_pdf))
  }, error = function(e) {
    log_message(paste("Error during visualization:", e$message), level="ERROR")
  })

  log_message("Analysis complete!")
}

if (!interactive()) {
  main(args)
}