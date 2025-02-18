# Script Name: ProBLASTr.R
# Description: An R script for protein BLAST analysis and homology verification
#
# Copyright (C) [2025] [University of Cologne]
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#
# Contact Information:
# Lee Weinand
# [Your Email]
# Univeristy of Cologne
#
# Version: 0.8
# Last Updated: [2025.02.18]
#
# This script is part of [Project Name]
# Project Repository: [URL to your repository]
# Documentation: [URL to documentation if available]
#
# Dependencies:
# R version 4.4.1 (2024-06-14)
# Platform: x86_64-pc-linux-gnu
# Running under: Ubuntu 24.10
#
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
# LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.28.so;  LAPACK version 3.12.0
#
# locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
# [3] LC_TIME=en_IE.UTF-8        LC_COLLATE=en_US.UTF-8
# [5] LC_MONETARY=en_IE.UTF-8    LC_MESSAGES=en_US.UTF-8
# [7] LC_PAPER=en_IE.UTF-8       LC_NAME=C
# [9] LC_ADDRESS=C               LC_TELEPHONE=C
# [11] LC_MEASUREMENT=en_IE.UTF-8 LC_IDENTIFICATION=C
#
# time zone: Europe/Berlin
# tzcode source: system (glibc)
#
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets
# [6] methods   base
#
# other attached packages:
# [1] lubridate_1.9.4 forcats_1.0.0   stringr_1.5.1
# [4] dplyr_1.1.4     purrr_1.0.4     tidyr_1.3.1
# [7] tibble_3.2.1    ggplot2_3.5.1   tidyverse_2.0.0
# [10] readr_2.1.5
#
# loaded via a namespace (and not attached):
# [1] crayon_1.5.3      vctrs_0.6.5       cli_3.6.3
# [4] rlang_1.1.5       stringi_1.8.4     generics_0.1.3
# [7] glue_1.8.0        bit_4.5.0.1       colorspace_2.1-1
# [10] hms_1.1.3         scales_1.3.0      grid_4.4.1
# [13] munsell_0.5.1     tzdb_0.4.0        lifecycle_1.0.4
# [16] compiler_4.4.1    timechange_0.3.0  pkgconfig_2.0.3
# [19] rstudioapi_0.17.1 R6_2.5.1          tidyselect_1.2.1
# [22] vroom_1.6.5       pillar_1.10.1     parallel_4.4.1
# [25] magrittr_2.0.3    tools_4.4.1       withr_3.0.2
# [28] bit64_4.6.0-1     gtable_0.3.6
#
# Usage:
# Provide a brief example of how to use this script
#
# Example:
# Rscript your_script_name.R [arguments]
#
# Notes:
# Any additional information that users should know about this script
# --------------------------------------------------------------------------

# Load required libraries
library(tidyverse)
library(GenomicRanges)
library(Biostrings)
library(GenomicFeatures)
library(rtracklayer)
library(Rsamtools)
library(stringi)

# Error handling configuration
options(error = NULL)

#' Create BLAST database and run tBLASTn search
#'
#' This function creates a BLAST database from genome files and performs tBLASTn searches
#' using query protein sequences against these databases.
#'
#' @param genome_dir Directory containing genome FASTA files ending with ".fasta"
#' @param query_dir Directory containing query protein FASTA files ending with ".fasta"
#' @param output_dir Directory for output files
#' @param evalue E-value threshold for BLAST searches
#' @param dbtype Type of BLAST database ("nucl" or "prot")
#' @param gen_file_pattern Pattern to match genome files
#' @param qur_file_pattern Pattern to match query files
#' @return DataFrame containing paths to query and BLAST hit files
create_blastdb_and_run_tblastn <- function(genome_dir,
                                           query_dir,
                                           output_dir,
                                           evalue = 1e-5,
                                           dbtype = "nucl",
                                           gen_file_pattern = "\\.fasta$",
                                           qur_file_pattern = "\\.fasta$") {
  # Input validation
  if (!(dir.exists(genome_dir))) {
    stop(paste0("Genome directory '", genome_dir, "' not found!"))
  }
  if (!(dir.exists(query_dir))) {
    stop(paste0("Query directory '", query_dir, "' not found!"))
  }
  
  # Create directories for BLAST database
  blast_database_dir <- file.path(output_dir, "blast_database")
  if (!dir.exists(blast_database_dir)) {
    dir.create(blast_database_dir, recursive = TRUE)
  }
  
  # List and validate genome files
  genome_files <- list.files(genome_dir, pattern = gen_file_pattern, full.names = TRUE)
  if (length(genome_files) == 0) {
    stop("No genome files found!")
  }
  
  message("Creating BLAST databases...")
  db_paths <- character(length(genome_files))
  
  # Create BLAST database for each genome file
  for (i in seq_along(genome_files)) {
    genome_file <- genome_files[i]
    db_name <- file.path(blast_database_dir,
                         paste0(tools::file_path_sans_ext(basename(genome_file)), "_db"))
    
    system2(
      "makeblastdb",
      args = c(
        "-in",
        genome_file,
        "-dbtype",
        dbtype,
        "-out",
        db_name,
        "-parse_seqids"
      )
    )
    db_paths[i] <- db_name
  }
  
  # List and validate query files
  query_files <- list.files(query_dir, pattern = qur_file_pattern, full.names = TRUE)
  if (length(query_files) == 0) {
    stop("No query files found!")
  }
  
  message("Running BLAST searches...")
  
  # Create directory for BLAST hits
  blast_hits_dir <- file.path(output_dir, "blast_hits")
  if (!dir.exists(blast_hits_dir)) {
    dir.create(blast_hits_dir, recursive = TRUE)
  }
  
  # Initialize index for tracking BLAST queries and results
  index_blast <- data.frame(query_path = character(0),
                            blast_hit_file_path = character(0))
  
  # Run BLAST for each query against each database
  for (query in query_files) {
    query_name <- tools::file_path_sans_ext(basename(query))
    for (db in unlist(db_paths)) {
      cat(paste0("Query: '", query, "'\n", "Db: '", db, "'\n"))
      
      db_name <- basename(db)
      output_file <- file.path(blast_hits_dir,
                               paste0(query_name, "_in_", db_name, ".tsv"))
      
      # Run tBLASTn
      system2(
        "tblastn",
        args = c(
          "-query",
          query,
          "-db",
          db,
          "-out",
          output_file,
          "-evalue",
          evalue,
          "-outfmt",
          "\"6 qseqid sseqid percent_identity length mismatch gapopen qstart qend sstart send evalue bitscore\""
        )
      )
      
      # Handle empty results
      if (file.exists(output_file) && file.size(output_file) == 0) {
        file.remove(output_file)
      }
      
      # Record successful BLAST runs
      if (file.exists(output_file) && file.size(output_file) > 0) {
        index_blast <- rbind(
          index_blast,
          data.frame(
            query_path = query,
            blast_hit_file_path = output_file
          )
        )
      }
    }
  }
  return(index_blast)
}


#' Find genes from BLAST hits
#'
#' This function identifies genes in the genomic regions corresponding to BLAST hits
#' 
#' @param blast_hit DataFrame containing a BLAST hit with the column sseqid
#' @param gff_file Directory containing Gene Annotation gff files ending with ".gff3"
#' @param min_seq_len Minimus Sequence length for a Sequence in a BLAST hit
#' @return DataFrame containing BLAST hit and Gene names from the corresponfing gff_file if it has at least one row otherwise NULL
find_gene_from_blast_hit <- function(blast_hit, gff_file, min_seq_len) {
  # Filter BLAST hits based on overlap and quality criteria
  filtered_blast_hit <- blast_hit |>
    mutate(low_start = pmin(sstart, send),
           high_end = pmax(sstart, send)) |>
    group_by(sseqid) |>
    arrange(low_start) |>
    mutate(
      prev_end = lag(high_end, default = head(high_end, 1)),
      distance_to_prev = low_start - prev_end,
      is_overlap = low_start <= prev_end
    ) |>
    filter((distance_to_prev <= min_seq_len |
              is.na(distance_to_prev)) &
             !is_overlap) |>
    filter(bitscore >= 50) |>
    ungroup() |>
    select(-prev_end,-distance_to_prev,-is_overlap,
           low_start,
           high_end)
  
  # Match BLAST hits with gene annotations
  result <- inner_join(filtered_blast_hit,
                       gff_file,
                       by = "sseqid",
                       relationship = "many-to-many") |>
    select(-contains("idk")) |>
    mutate(seq_in_gene = ifelse(sstart >= gstart &
                                  send <= gend, TRUE, FALSE)) |>
    filter(seq_in_gene == TRUE) |>
    select(-seq_in_gene)
  
  if (nrow(result) > 0) {
    return(result)
  } else {
    return(NULL)
  }
}

#' Convert gene sequence to protein
#'
#' This function extracts and translates gene sequences to protein sequences
#'
#' @param genome_dir Directory containing genome FASTA files ending with ".fasta"
#' @param gff_file Directory containing gene annotation gff files ending with ".gff3"
#' @param transcript_id Transcript ID corresponding to a BLAST hit
#' @param output_dir Directory for output files
#' @param gene_hit_path Path to BLAST gene hit files
#' @return Path to the protein sequence
gene_seq_to_prot <- function(genome_file,
                             gff_file,
                             transcript_id,
                             output_dir,
                             gene_hit_path) {
  message("Processing gene to protein conversion...")
  message("genome_file = ", genome_file)
  message("gff_file = ", gff_file)
  message("transcript_id = ", transcript_id)
  message("output_dir = ", output_dir)
  
  # Read and filter GFF file
  gff <- import(gff_file)
  gff <- gff[gff$type == "CDS"]
  cds <- gff[as.character(gff$Parent) == transcript_id]
  
  # Sort CDS based on strand
  if (unique(strand(cds)) == "-") {
    cds <- sort(cds, decreasing = TRUE)
  } else {
    cds <- sort(cds)
  }
  
  # Extract sequences
  fa <- FaFile(genome_file)
  cds_seqs <- sapply(seq_along(cds), function(i) {
    chr <- as.character(seqnames(cds[i]))
    start <- start(cds[i])
    end <- end(cds[i])
    
    seq <- scanFa(fa, param = GRanges(
      seqnames = chr,
      ranges = IRanges(start = start, end = end)
    ))
    
    if (as.character(strand(cds[i])) == "-") {
      seq <- reverseComplement(seq)
    }
    
    return(as.character(seq))
  })
  
  # Join CDS sequences
  full_transcript <- DNAString(paste(cds_seqs, collapse = ""))
  
  # Find and translate longest ORF
  find_and_translate_longest_ORF <- function(sequence) {
    atg_positions <- start(matchPattern("ATG", sequence))
    longest_protein <- ""
    
    for (start_pos in atg_positions) {
      current_seq <- subseq(sequence, start_pos)
      
      if (length(current_seq) >= 3) {
        protein <- translate(current_seq)
        stop_pos <- min(c(start(matchPattern(
          "*", protein
        )), length(protein) + 1))
        current_protein <- substr(as.character(protein), 1, stop_pos - 1)
        
        if (nchar(current_protein) > nchar(longest_protein)) {
          longest_protein <- current_protein
        }
      }
    }
    return(longest_protein)
  }
  
  # Get protein sequence and gene name
  protein <- find_and_translate_longest_ORF(full_transcript)
  gene_name <- unique(as.character(mcols(cds)$Parent))
  
  # Create output directories
  protein_output_dir_path <- file.path(output_dir, "protein")
  transcript_output_dir_path <- file.path(output_dir, "transcript")
  for (dir in c(protein_output_dir_path, transcript_output_dir_path)) {
    if (!dir.exists(dir)) {
      dir.create(dir)
    }
  }
  
  # Write protein sequence
  protein_header <- paste0(
    ">",
    gene_name,
    " Gene hit: ",
    basename(gene_hit_path),
    " Genome: ",
    basename(genome_file),
    " gff_file: ",
    basename(gff_file),
    " - Longest Open Reading Frame"
  )
  protein_out_path <- file.path(protein_output_dir_path,
                                paste0(gene_name, "_protein.fasta"))
  writeLines(c(protein_header, protein), protein_out_path)
  
  # Write transcript sequence
  transcript_header <- paste0(
    ">",
    gene_name,
    " Gene hit: ",
    basename(gene_hit_path),
    " Genome: ",
    basename(genome_file),
    " gff_file: ",
    basename(gff_file),
    " - Full Transcript Sequence"
  )
  writeLines(
    c(transcript_header, as.character(full_transcript)),
    file.path(
      transcript_output_dir_path,
      paste0(gene_name, "_transcript.fasta")
    )
  )
  
  return(protein_out_path)
}


#' Verify homology between tomato and potato sequences
#'
#' This function performs detailed homology analysis between tomato and potato sequences.
#' It creates BLAST databases, performs comparisons, and generates comprehensive reports
#' of the relationships between sequences.
#'
#' @param output_dir Directory for output files
#' @param tomato_proteom Path to tomato proteome reference file
#' @param potato_prot Directory containing potato protein sequences
#' @param tomato_prot Directory containing tomato protein sequences
#' @param evalue E-value threshold for BLAST searches
#' @param index_de DataFrame containing file relationship information
#' @return List containing detailed comparison and summary statistics
verify_tomato_potato_homology <- function(output_dir,
                                          tomato_proteom = "tomato/UP000004994_4081.fasta",
                                          potato_prot,
                                          tomato_prot,
                                          evalue = 1e-5,
                                          index_de) {
  # Validate input files and directories
  for (file in c(tomato_proteom)) {
    if (!file.exists(file)) {
      stop("File not found: ", file)
    }
  }
  for (dir in c(potato_prot, tomato_prot)) {
    if (!dir.exists(dir)) {
      stop("Directory not found: ", dir)
    }
  }
  
  # Create output directory structure
  homologue_verifi_path <- file.path(output_dir, "homologue_verification")
  tomato_proteom_db_path <- file.path(homologue_verifi_path, "blastp_database")
  potato_prot_blast_hit_path <- file.path(homologue_verifi_path, "potato_prot_blast_hit")
  tomato_prot_blast_hit_path <- file.path(homologue_verifi_path, "tomato_prot_blast_hit")
  
  for (dir in c(
    homologue_verifi_path,
    tomato_proteom_db_path,
    potato_prot_blast_hit_path,
    tomato_prot_blast_hit_path
  )) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
  }
  
  # Create BLAST database from tomato proteome
  message("Creating BLAST database from tomato proteome...")
  proteom_db_path <- file.path(tomato_proteom_db_path,
                               paste0(tools::file_path_sans_ext(basename(tomato_proteom)), "_db"))
  system2(
    "makeblastdb",
    args = c(
      "-in",
      tomato_proteom,
      "-dbtype",
      "prot",
      "-out",
      proteom_db_path
    )
  )
  
  # Helper function to perform BLAST searches
  blastp <- function(query, db, out, evalue = 1e-5) {
    message("Running BLAST for query: ", basename(query))
    tryCatch({
      system2(
        "blastp",
        args = c(
          "-query",
          query,
          "-db",
          db,
          "-out",
          out,
          "-evalue",
          evalue,
          "-outfmt",
          "\"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\""
        )
      )
      
      if (file.exists(out)) {
        if (file.size(out) == 0) {
          message("No BLAST hits found for: ", basename(query))
          file.remove(out)
          return(FALSE)
        }
        return(TRUE)
      }
      return(FALSE)
    }, error = function(e) {
      warning("BLAST failed for query: ", query, "\nError: ", e$message)
      return(FALSE)
    })
  }
  
  # Process potato proteins
  message("\nProcessing potato proteins...")
  potato_results <- data.frame()
  potato_files <- list.files(potato_prot, pattern = "_protein.fasta$", full.names = TRUE)
  
  for (protein_seq in potato_files) {
    output_file <- file.path(
      potato_prot_blast_hit_path,
      paste0(
        tools::file_path_sans_ext(basename(protein_seq)),
        "_in_",
        basename(tools::file_path_sans_ext(tomato_proteom)),
        ".tsv"
      )
    )
    
    if (blastp(protein_seq, proteom_db_path, output_file, evalue)) {
      blast_results <- read.table(
        output_file,
        sep = "\t",
        col.names = c(
          "qseqid",
          "sseqid",
          "pident",
          "length",
          "mismatch",
          "gapopen",
          "qstart",
          "qend",
          "sstart",
          "send",
          "evalue",
          "bitscore"
        )
      )
      blast_results$source_file_potato <- basename(protein_seq)
      potato_results <- rbind(potato_results, blast_results)
    }
  }
  
  # Process tomato proteins
  message("\nProcessing tomato proteins...")
  tomato_results <- data.frame()
  tomato_files <- list.files(tomato_prot, pattern = "\\.fasta$", full.names = TRUE)
  
  for (protein_seq in tomato_files) {
    output_file <- file.path(
      tomato_prot_blast_hit_path,
      paste0(
        tools::file_path_sans_ext(basename(protein_seq)),
        "_in_",
        basename(tools::file_path_sans_ext(tomato_proteom)),
        ".tsv"
      )
    )
    
    if (blastp(protein_seq, proteom_db_path, output_file, evalue)) {
      blast_results <- read.table(
        output_file,
        sep = "\t",
        col.names = c(
          "qseqid",
          "sseqid",
          "pident",
          "length",
          "mismatch",
          "gapopen",
          "qstart",
          "qend",
          "sstart",
          "send",
          "evalue",
          "bitscore"
        )
      )
      blast_results$source_file_tomato <- basename(protein_seq)
      tomato_results <- rbind(tomato_results, blast_results)
    }
  }
  
  # Generate comparison results
  if (nrow(potato_results) > 0 && nrow(tomato_results) > 0) {
    # Get best hits for each query
    potato_best_hits <- potato_results |>
      group_by(source_file_potato) |>
      slice_max(order_by = bitscore, n = 1)
    
    tomato_best_hits <- tomato_results |>
      group_by(source_file_tomato) |>
      slice_max(order_by = bitscore, n = 1)
    
    # Compare hits and join with index_de information
    comparison <- potato_best_hits |>
      left_join(tomato_best_hits,
                by = "sseqid",
                suffix = c("_potato", "_tomato")) |>
      left_join(
        select(index_de, blast_hit_file_path, gene_hit_path, query_path),
        by = c("source_file_potato" = "blast_hit_file_path")
      ) |>
      mutate(
        match_quality = case_when(
          is.na(bitscore_tomato) ~ "No Match",
          pident_potato >= 90 &
            pident_potato / pident_tomato >= 0.9 ~ "Strong Homolog",
          pident_potato >= 70 &
            pident_potato / pident_tomato >= 0.7 ~ "Moderate Homolog",
          pident_potato >= 50 &
            pident_potato / pident_tomato >= 0.5 ~ "Weak Homolog",
          TRUE ~ "Poor Match"
        ),
        similarity_ratio = pident_potato / pident_tomato,
        coverage_ratio = length_potato / length_tomato
      )
    
    # Generate summary statistics
    summary_stats <- comparison |>
      group_by(match_quality) |>
      summarise(
        count = n(),
        avg_similarity = mean(similarity_ratio, na.rm = TRUE),
        avg_coverage = mean(coverage_ratio, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Write results to files
    write_csv(
      comparison,
      file.path(
        homologue_verifi_path,
        "homology_comparison_detailed.csv"
      )
    )
    write_csv(
      summary_stats,
      file.path(homologue_verifi_path, "homology_summary_stats.csv")
    )
    
    return(list(
      detailed_comparison = comparison,
      summary_statistics = summary_stats
    ))
  } else {
    warning("No BLAST results found for comparison")
    return(NULL)
  }
}

#' Main execution function for ProBLASTr pipeline
#'
#' This function runs the entire ProBLASTr workflow, from initial BLAST searches 
#' through gene finding, protein sequence translation and homology verification.
#'
#' @param genome_dir Directory containing genome files
#' @param query_dir Directory containing query files
#' @param gff_file Directory containing gene annotation gff files ending with ".gff3"
#' @param output_dir Directory for output files
#' @param evalue E-value threshold for BLAST searches
#' @param min_seq_len Minimum sequence length to consider
#' @param tomato_proteom Path to tomato proteome reference file
main <- function(genome_dir,
                 query_dir,
                 gff_dir,
                 output_dir,
                 evalue = 1e-5,
                 min_seq_len = 50000,
                 tomato_proteom = "tomato/UP000004994_4081.fasta") {
  # Check required parameters
  if (missing(genome_dir) ||
      missing(query_dir) ||
      missing(gff_dir) ||
      missing(output_dir) ||
      missing(tomato_proteom)) {
    stop("Required parameter missing!")
  }
  
  # Read and validate index file
  if (!file.exists("index.csv")) {
    stop("The file index.csv does not exist.")
  }
  
  index_from_file <- read_csv(file = "index.csv", show_col_types = FALSE) |>
    na.omit()
  
  if (!("gff_file_path" %in% colnames(index_from_file))) {
    stop("index.csv must contain a column named 'gff_file_path'")
  }
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Initialize empty DataFrame for file paths
  index <- data.frame(
    blast_hit_file_path = character(0),
    gff_file_path = character(0),
    gene_hit_path = character(0)
  )
  
  # Run initial BLAST searches
  index_blast <- create_blastdb_and_run_tblastn(genome_dir, query_dir, output_dir, evalue)
  
  # Process BLAST hits
  blast_hits <- list.files(path = file.path(output_dir, "blast_hits"),
                           pattern = "\\.tsv$")
  
  # Create gene hits directory
  gene_hit_dir_path <- file.path(output_dir, "gene_hits")
  if (!dir.exists(gene_hit_dir_path)) {
    dir.create(gene_hit_dir_path, recursive = TRUE)
  }
  
  # Process each BLAST hit
  for (blast_hit in blast_hits) {
    message("Searching gff3 files for annotated genes for blast hit: '",
            blast_hit,
            "'...")
    
    blast_hit_file_path <- file.path(output_dir, "blast_hits", blast_hit)
    
    # Read and process BLAST results
    blast_hit_table <- read.table(
      file = blast_hit_file_path,
      sep = "\t",
      col.names = c(
        "qseqid",
        "sseqid",
        "percent_identity",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore"
      )
    ) |>
      mutate(sseqid = as.character(sseqid))
    
    # Find matching GFF files
    blast_hit_unique_sseqid <- paste(unique(blast_hit_table$sseqid), collapse = "|")
    matching_gff_file <- system2("egrep",
                                 args = c(
                                   "-l",
                                   shQuote(blast_hit_unique_sseqid),
                                   file.path(gff_dir, "*")
                                 ),
                                 stdout = TRUE)
    
    # Process each matching GFF file
    for (gff_file_hit in basename(matching_gff_file)) {
      gff_file_path <- file.path(gff_dir, gff_file_hit)
      
      # Read and process GFF file
      gff_file_table <- read.table(
        file = gff_file_path,
        sep = "\t",
        col.names = c(
          "sseqid",
          "source",
          "type",
          "gstart",
          "gend",
          "score",
          "strand",
          "phase",
          "gname"
        )
      ) |>
        mutate(sseqid = as.character(sseqid)) |>
        filter(type %in% c("mRNA", "transcript"))
      
      # Find genes from BLAST hits
      result <- find_gene_from_blast_hit(
        blast_hit = blast_hit_table,
        gff_file = gff_file_table,
        min_seq_len = min_seq_len
      )
      
      if (!is.null(result) &&
          is.data.frame(result) && nrow(result) > 0) {
        gene_hit_file_path <- file.path(
          gene_hit_dir_path,
          paste0(
            tools::file_path_sans_ext(blast_hit),
            "_with_",
            tools::file_path_sans_ext(gff_file_hit),
            "_gene_hit.csv"
          )
        )
        
        write_csv(result, file = gene_hit_file_path)
        
        # Update index
        index <- rbind(
          index,
          data.frame(
            blast_hit_file_path = blast_hit_file_path,
            gff_file_path = gff_file_path,
            gene_hit_path = gene_hit_file_path
          )
        )
      }
    }
  }
  
  # Save backup of current index
  write_csv(index, file = "index.csv.bak")
  
  # Update file paths in index
  # Update file paths in index
  index_from_file$gff_file_path <- file.path(gff_dir, index_from_file$gff_file_path)
  index_from_file$genome_file_path <- file.path(genome_dir, index_from_file$genome_file_path)
  
  # Join indices to create complete relationship map
  index_de <- left_join(index, index_from_file, by = "gff_file_path") |>
    na.omit() |>
    left_join(index_blast, by = "blast_hit_file_path")
  
  # Set up for protein sequence generation
  protein_output_dir_path <- file.path(output_dir, "protein")
  
  # Initialize tracking for gene to protein conversion
  gene_to_prot_index <- data.frame(
    gene_hit_path = character(0),
    transcript_id = character(0),
    protein_out_path = character(0)
  )
  
  # Process each gene hit for protein sequence generation
  message("\nConverting genes to proteins...")
  for (i in seq_along(index_de$gene_hit_path)) {
    cat(sprintf(
      "\nProcessing gene hit %d of %d\n",
      i,
      length(index_de$gene_hit_path)
    ))
    
    gene_hit_path <- index_de$gene_hit_path[i]
    gff_file_path <- index_de$gff_file_path[i]
    genome_file_path <- index_de$genome_file_path[i]
    
    # Clean output_dir from paths for consistency
    gff_file_path <- sub(
      pattern = paste0("^", output_dir, "[/\\]"),
      replacement = "",
      gff_file_path
    )
    genome_file_path <- sub(
      pattern = paste0("^", output_dir, "[/\\]"),
      replacement = "",
      genome_file_path
    )
    
    # Read and process gene hit information
    gene_hit_file_path_df <- read_csv(file = gene_hit_path, show_col_types = FALSE)
    transcript_id_raw <- gene_hit_file_path_df$gname
    transcript_ids <- unique(stri_extract_first(transcript_id_raw, regex = "(?<=ID=)[^;]+"))
    
    # Process each transcript
    for (transcript_id in transcript_ids) {
      message(sprintf("Processing transcript: %s", transcript_id))
      
      tryCatch({
        protein_out_path <- gene_seq_to_prot(
          genome_file = genome_file_path,
          gff_file = gff_file_path,
          transcript_id = transcript_id,
          output_dir = output_dir,
          gene_hit_path = gene_hit_path
        )
        
        # Record successful protein generation
        gene_to_prot_index <- rbind(
          gene_to_prot_index,
          data.frame(
            gene_hit_path = gene_hit_path,
            transcript_id = transcript_id,
            protein_out_path = protein_out_path
          )
        )
      }, error = function(e) {
        warning(sprintf(
          "Error processing transcript %s: %s",
          transcript_id,
          e$message
        ))
      })
    }
  }
  
  # Update final index with protein information
  message("\nUpdating final index with protein information...")
  tryCatch({
    index_de <- index_de |>
      left_join(gene_to_prot_index, by = "gene_hit_path")
  }, error = function(e) {
    warning(sprintf("Error joining protein information: %s", e$message))
  })
  
  # Save complete index
  message("\nSaving final index...")
  write_csv(index_de, file = file.path(output_dir, "final_index.csv"))
  
  # Perform homology verification
  message("\nPerforming homology verification...")
  tryCatch({
    results <- verify_tomato_potato_homology(
      output_dir = output_dir,
      tomato_proteom = tomato_proteom,
      potato_prot = protein_output_dir_path,
      tomato_prot = query_dir,
      evalue = evalue,
      index_de = index_de
    )
    
    message("\nHomology verification completed successfully")
    return(results)
    
  }, error = function(e) {
    warning(sprintf("Error in homology verification: %s", e$message))
    return(NULL)
  })
}

# Example usage of the main function
if (T) {
  # Set to TRUE to run
  results <- main(
    genome_dir = "genomes",
    query_dir = "meiotic_genes_protein_fasta",
    gff_dir = "gff_files",
    output_dir = "out_2024_02_18",
    evalue = 1e-5
  )
}