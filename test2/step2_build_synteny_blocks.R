#!/usr/bin/env Rscript

# =============================================================================
#
# SCRIPT: Ultrafast Synteny Block Reconstruction from RBH Data
#
# BIOLOGICAL CONTEXT & PURPOSE:
# This script forms the foundational step of ancestral genome reconstruction.
# It takes Reciprocal Best Hit (RBH) gene pairs and chains them together to
# form "synteny blocks".
#
# ALorithmic Innovation:
# This version (v4.7) uses the explicit '-' delimiter for robust species name
# parsing, based on user-provided file naming conventions.
#
# VERSION: 4.7 (Final Parsing Fix)
# DATE: 2025-07-03
# AUTHOR: Biols9527 (Refactored by Gemini AI)
#
# =============================================================================

# --- 1. Load Libraries ---
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(parallel)
})

# --- 2. Command Line Argument Parsing ---
option_list <- list(
  make_option(c("--rbh_dir"), type="character",
              help="Directory containing RBH files (one per species pair)"),
  make_option(c("--output_dir"), type="character", default=".",
              help="Output directory for this module [default: current directory]"),
  make_option(c("--min_block_genes"), type="integer", default=5,
              help="Minimum number of genes required to form a synteny block [default: %default]"),
  make_option(c("--max_rank_gap"), type="integer", default=20,
              help="Maximum allowed gap in gene ranks within a block [default: %default]"),
  make_option(c("--threads"), type="integer", default=4,
              help="Number of parallel threads to use [default: %default]"),
  make_option(c("--has_header"), action="store_true", default=FALSE,
              help="Specify if the input RBH files have a header row")
)

parser <- OptionParser(option_list=option_list, usage = "%prog --rbh_dir <path> [options]")
args <- parse_args(parser)

# --- 3. Utility Functions ---
log_message <- function(msg, level="INFO") {
  cat(sprintf("[%s] [%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, msg))
}

ensure_dir <- function(dir) {
  if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
}

# --- 4. Core Logic: Ultrafast Synteny Block Identification ---

identify_synteny_blocks_fast <- function(rbh_dt, min_genes, max_gap) {
  rbh_dt <- copy(rbh_dt)

  if (nrow(rbh_dt) < min_genes) return(data.table())

  rbh_dt[, rank1 := rank(pos1, ties.method="first")]
  rbh_dt[, rank2 := rank(pos2, ties.method="first")]

  strands <- list(forward = 1, reverse = -1)
  all_blocks <- list()

  for (strand_name in names(strands)) {
    direction <- strands[[strand_name]]

    if (direction == -1) {
      rbh_dt[, rank2_proc := (max(rank2) + 1) - rank2]
    } else {
      rbh_dt[, rank2_proc := rank2]
    }

    setkey(rbh_dt, rank1)
    rbh_dt[, rank_diff := rank1 - rank2_proc]
    rbh_dt[, block_start := (rank_diff - shift(rank_diff, fill=rank_diff[1])) > max_gap]
    rbh_dt[1, block_start := TRUE]
    rbh_dt[, block_id_strand := cumsum(block_start)]

    rbh_dt[, gene_count := .N, by = block_id_strand]
    valid_blocks_dt <- rbh_dt[gene_count >= min_genes]

    if (nrow(valid_blocks_dt) > 0) {
      block_summary <- valid_blocks_dt[, .(
        gene_count = .N,
        start1 = min(pos1), end1 = max(pos1),
        start2 = min(pos2), end2 = max(pos2),
        span1 = max(pos1) - min(pos1) + 1,
        span2 = max(pos2) - min(pos2) + 1,
        synteny_score = cor(rank1, rank2, method="spearman")
      ), by = block_id_strand]

      block_summary[, direction := direction]
      all_blocks[[strand_name]] <- block_summary
    }
  }

  if (length(all_blocks) == 0) return(data.table())

  final_blocks <- rbindlist(all_blocks, use.names = TRUE)
  final_blocks[, block_id_strand := NULL]

  return(final_blocks)
}

process_rbh_file <- function(file_path, min_block_genes, max_rank_gap, has_header) {
  tryCatch({
    if (has_header) {
        rbh_data <- fread(file_path, header = TRUE)
        setnames(rbh_data, 1:4, c("pos1", "pos2", "chr1", "chr2"))
    } else {
        rbh_data <- fread(file_path, header = FALSE, col.names=c("pos1", "pos2", "chr1", "chr2"))
    }

    if (nrow(rbh_data) == 0) return(NULL)

    rbh_data <- rbh_data[!is.na(pos1) & !is.na(pos2) & chr1 != "" & chr2 != ""]
    rbh_data[, `:=`(pos1 = as.numeric(pos1), pos2 = as.numeric(pos2))]
    rbh_data <- na.omit(rbh_data)

    if (nrow(rbh_data) < min_block_genes) return(NULL)

    all_blocks_list <- rbh_data[, {
        identify_synteny_blocks_fast(.SD, min_block_genes, max_rank_gap)
      }, by = .(chr1, chr2)]

    return(all_blocks_list)

  }, error = function(e) {
    stop(sprintf("Failed to process file %s: %s", basename(file_path), e$message))
    return(NULL)
  })
}

# --- 5. Main Execution Block ---
main <- function(args) {
  ensure_dir(args$output_dir)
  log_message("=== Ultrafast Synteny Block Reconstruction (v4.7) ===")

  rbh_files <- list.files(args$rbh_dir, full.names=TRUE, pattern="\\.(tsv|csv|txt)$")
  if (length(rbh_files) == 0) {
    stop(paste("No RBH files found in directory:", args$rbh_dir))
  }
  log_message(paste("Found", length(rbh_files), "RBH files to process."))

  log_message(paste("Processing files in parallel using", args$threads, "threads..."))
  cl <- makeCluster(args$threads)
  
  clusterExport(cl, varlist=c("args", "identify_synteny_blocks_fast", "process_rbh_file"))
  clusterEvalQ(cl, library(data.table))

  results_list <- parLapply(cl, rbh_files, function(file_path) {
    fname_noext <- tools::file_path_sans_ext(basename(file_path))
    
    # FINAL CORRECTED LOGIC: Split specifically on "-" as per user's file naming convention.
    species_pair <- strsplit(fname_noext, "-")[[1]]
    
    if (length(species_pair) == 2) {
        sp1 <- species_pair[1]
        sp2 <- species_pair[2]
    } else {
        log_message(paste("Could not parse species pair from filename (expected 'speciesA-speciesB.ext'):", basename(file_path)), "WARN")
        return(NULL) # Skip files with unexpected naming
    }

    blocks <- process_rbh_file(file_path, args$min_block_genes, args$max_rank_gap, args$has_header)

    if (!is.null(blocks) && nrow(blocks) > 0) {
      blocks[, `:=`(species1 = sp1, species2 = sp2)]
      return(blocks)
    }
    return(NULL)
  })

  stopCluster(cl)

  all_synteny_blocks <- rbindlist(Filter(Negate(is.null), results_list))

  if (nrow(all_synteny_blocks) == 0) {
    log_message("Analysis complete, but no synteny blocks were identified across all files.", "WARN")
    return()
  }

  all_synteny_blocks[, block_id := .I]
  setcolorder(all_synteny_blocks, c("block_id", "species1", "chr1", "species2", "chr2"))

  log_message(paste("Successfully identified", nrow(all_synteny_blocks), "synteny blocks."))

  output_file <- file.path(args$output_dir, "synteny_blocks_all.tsv")
  fwrite(all_synteny_blocks, output_file, sep="\t")
  log_message(paste("All synteny blocks saved to:", output_file))

  summary_text <- c(
    "=== Synteny Block Analysis Summary ===",
    paste("Timestamp:", Sys.time()),
    paste("Total RBH files processed:", length(rbh_files)),
    paste("Total synteny blocks found:", nrow(all_synteny_blocks)),
    paste("Average genes per block:", sprintf("%.2f", mean(all_synteny_blocks$gene_count))),
    paste("Average synteny score (Spearman's rho):", sprintf("%.3f", mean(all_synteny_blocks$synteny_score, na.rm=TRUE)))
  )
  writeLines(summary_text, file.path(args$output_dir, "analysis_summary.txt"))

  log_message("Analysis complete.")
}

# --- Script Entry Point ---
if (!interactive()) {
  main(args)
}