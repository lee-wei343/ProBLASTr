#!/usr/bin/env Rscript
#
# Script Name: ProBLASTr.R
# Description: A pipeline for protein BLAST analysis to identify homologous genes
#              across genomic datasets
#
# Copyright (C) [2025] [University of Cologne]
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Version: 1.1
# Last Updated: 2025-02-28
#
# Contact Information:
# Lee Weinand
# University of Cologne
#
# Dependencies:
# R version 4.4.1 or higher
# Required packages: tidyverse, GenomicRanges, Biostrings, GenomicFeatures, 
#                    rtracklayer, Rsamtools, stringi
# External tools: BLAST+ suite (blastn, tblastn, blastp, makeblastdb)


# Libraries ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(Biostrings)
  library(GenomicFeatures)
  library(rtracklayer)
  library(Rsamtools)
  library(stringi)
  library(msa)
})

check_dependencies <- function() {
  required_packages <- c(
    "tidyverse", "GenomicRanges", "Biostrings", "GenomicFeatures",
    "rtracklayer", "Rsamtools", "stringi", "msa"
  )
  
  missing_packages <- required_packages[!requireNamespace(required_packages, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    stop(
      "The following required packages are missing: ", 
      paste(missing_packages, collapse = ", "), 
      "\nPlease install them using: install.packages() or BiocManager::install()"
    )
  }
}

# Global Config -----------------------------------------------------------
# Set up error logging
setup_error_logging <- function() {
  log_file <- "problaster_error_log.txt"
  
  error_handler <- function() {
    cat(format(Sys.time()), 
        geterrmessage(), 
        "\n", 
        file = log_file, 
        append = TRUE)
  }
  
  options(error = error_handler)
  
  message("Error logging initialized. Errors will be saved to: ", log_file)
}


# BLAST Database and Search Functions -------------------------------------
#' Create BLAST database and run tBLASTn search
#'
#' This function creates a BLAST database from genome files and performs tBLASTn searches
#' using query protein sequences against these databases.
#'
#' @param genome_dir Directory containing genome FASTA files
#' @param query_dir Directory containing query protein FASTA files
#' @param output_dir Directory for output files
#' @param evalue E-value threshold for BLAST searches (default: 1e-5)
#' @param dbtype Type of BLAST database ("nucl" or "prot", default: "nucl")
#' @param gen_file_pattern Pattern to match genome files (default: "\\.fasta$")
#' @param qur_file_pattern Pattern to match query files (default: "\\.fasta$")
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
  
  # Create directories
  blast_database_dir <- file.path(output_dir, "blast_databases")
  blast_hits_dir <- file.path(output_dir, "blast_search_results")
  for (dir in c(blast_database_dir, blast_hits_dir)) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
  }
  
  # List and validate genome files
  genome_files <- list.files(genome_dir, pattern = gen_file_pattern, full.names = TRUE)
  if (length(genome_files) == 0) {
    stop("No genome files found matching the pattern: ", gen_file_pattern)
  }
  
  message("Creating BLAST databases from ", length(genome_files), " genome files...")
  db_paths <- character(length(genome_files))
  
  # Create BLAST database for each genome file
  for (i in seq_along(genome_files)) {
    genome_file <- genome_files[i]
    genome_name <- tools::file_path_sans_ext(basename(genome_file))
    db_name <- file.path(blast_database_dir, paste0(genome_name, "_blast_db"))
    
    message("  Creating database for: ", genome_name)
    
    system2(
      "makeblastdb",
      args = c(
        "-in", genome_file,
        "-dbtype", dbtype,
        "-out", db_name
      )
    )
    
    # Store the database path
    db_paths[i] <- db_name
  }
  
  # List and validate query files
  query_files <- list.files(query_dir, pattern = qur_file_pattern, full.names = TRUE)
  if (length(query_files) == 0) {
    stop("No query files found matching the pattern: ", qur_file_pattern)
  }
  
  message("Running BLAST searches with ", length(query_files), " query files...")
  
  # Initialize results dataframe
  blast_results <- data.frame(
    query_path = character(0),
    query_name = character(0),
    target_genome = character(0),
    blast_hit_file_path = character(0),
    stringsAsFactors = FALSE
  )
  
  # Run BLAST for each query against each database
  for (query in query_files) {
    query_name <- tools::file_path_sans_ext(basename(query))
    message("  Processing query: ", query_name)
    
    for (db in db_paths) {
      db_basename <- basename(db)
      genome_name <- sub("_blast_db$", "", db_basename)
      
      output_file <- file.path(blast_hits_dir,
                               paste0(query_name, "_vs_", genome_name, "_blast_hits.tsv"))
      
      # Run tBLASTn
      system2(
        "tblastn",
        args = c(
          "-query", query,
          "-db", db,
          "-out", output_file,
          "-evalue", evalue,
          "-outfmt", "\"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\""
        )
      )
      
      # Handle empty files
      if (file.exists(output_file) && file.size(output_file) == 0) {
        message("    No hits found for ", query_name, " in ", genome_name)
        file.remove(output_file)
      } else if (file.exists(output_file) && file.size(output_file) > 0) {
        # Record successful BLAST run
        blast_results <- rbind(
          blast_results,
          data.frame(
            query_path = query,
            query_name = query_name,
            target_genome = genome_name,
            blast_hit_file_path = output_file,
            stringsAsFactors = FALSE
          )
        )
        message("    Hits found for ", query_name, " in ", genome_name)
      }
    }
  }
  
  message("BLAST searches completed. Found hits in ", nrow(blast_results), " searches.")
  
  # Save the results to a file for reference
  blast_results_file <- file.path(output_dir, "blast_search_mapping.csv")
  write_csv(blast_results, blast_results_file)
  message("Saved BLAST search results to: ", blast_results_file)
  
  return(blast_results)
}


# Gene Identification Functions -------------------------------------------
#' Find genes corresponding to BLAST hits
#'
#' This function identifies genes in genomic regions matching BLAST hits
#' 
#' @param blast_hit DataFrame containing BLAST hits with required columns
#' @param gff_file DataFrame with gene annotation data
#' @return ?
find_gene_from_blast_hit <- function(blast_hit, gff_file) {
  
  # Initialize empty data frame to store results
  df <- data.frame()
  
  # Process each blast hit
  for (i in 1:nrow(blast_hit)) {
    blast_row <- blast_hit[i, ]
    
    # Get the blast hit coordinates
    blast_sstart <- blast_row$sstart
    blast_send <- blast_row$send
    
    # Filter GFF entries that overlap with the blast hit
    if (blast_sstart < blast_send) {
      # Forward strand - find GFF entries that overlap the blast hit
      gff_file_in_range <- gff_file |> 
        filter(gend >= blast_sstart & gstart <= blast_send)
    } else {
      # Reverse strand - find GFF entries that overlap the blast hit
      gff_file_in_range <- gff_file |> 
        filter(gend >= blast_send & gstart <= blast_sstart)
    }
    
    # Append the filtered GFF entries to the results
    if (nrow(gff_file_in_range) > 0) {
      df <- rbind(df, gff_file_in_range)
    }
  }
  
  # Only perform the join if we found matching GFF entries
  if (nrow(df) > 0) {
    result <- inner_join(blast_hit, df, 
                         by = "sseqid", 
                         relationship = "many-to-many")
    return(result)
  } else {
    return(NULL)
  }
}



# Sequence Processing Functions -------------------------------------------
#' Extract and translate gene sequences to protein
#'
#' This function extracts CDS from genome using gene coordinates and translates to protein
#'
#' @param genome_file Path to the genome FASTA file
#' @param gff_file Path to the gene annotation GFF file
#' @param transcript_id Transcript ID corresponding to a BLAST hit
#' @param output_dir Directory for output files
#' @param gene_hit_path Path to the BLAST gene hit file (for reference)
#' @return Path to the generated protein sequence file
gene_seq_to_prot <- function(genome_file,
                             gff_file,
                             transcript_id,
                             output_dir,
                             gene_hit_path) {
  message("Processing gene to protein conversion for transcript: ", transcript_id)
  
  # Read and filter GFF file to get CDS regions
  gff <- import(gff_file)
  gff <- gff[gff$type == "CDS"]
  cds <- gff[as.character(gff$Parent) == transcript_id]
  
  if (length(cds) == 0) {
    stop("No CDS features found for transcript ID: ", transcript_id)
  }
  
  # Sort CDS based on strand
  if (unique(strand(cds)) == "-") {
    cds <- sort(cds, decreasing = TRUE)
  } else {
    cds <- sort(cds)
  }
  
  # Extract sequences from genome
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
  
  # Join CDS sequences to form complete transcript
  full_transcript <- DNAString(paste(cds_seqs, collapse = ""))
  
  # Translate the transcript to protein
  protein <- translate(
    full_transcript,
    genetic.code = getGeneticCode("1"),  # Standard genetic code
    if.fuzzy.codon = "solve"             # Handle ambiguous codons
  )
  
  # Extract the protein sequence up to the first stop codon
  protein <- as.character(protein)
  stop_pos <- regexpr("\\*", protein)
  if (stop_pos > 0) {
    protein <- substr(protein, 1, stop_pos - 1)
  }
  
  # Get gene name from CDS
  gene_name <- unique(as.character(mcols(cds)$Parent))
  
  # Create output directories
  protein_output_dir_path <- file.path(output_dir, "protein_sequences")
  transcript_output_dir_path <- file.path(output_dir, "transcript_sequences")
  for (dir in c(protein_output_dir_path, transcript_output_dir_path)) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
  }
  
  # Create meaningful file names
  genome_name <- tools::file_path_sans_ext(basename(genome_file))
  
  # Write protein sequence file
  protein_header <- paste0(
    ">",
    gene_name,
    " | Source: ",
    genome_name,
    " | Reference: ",
    basename(gene_hit_path)
  )
  
  protein_out_path <- file.path(
    protein_output_dir_path,
    paste0(gene_name, "_", genome_name, "_protein.fasta")
  )
  
  writeLines(c(protein_header, protein), protein_out_path)
  
  # Write transcript sequence file
  transcript_header <- paste0(
    ">",
    gene_name,
    " | Source: ",
    genome_name,
    " | Reference: ",
    basename(gene_hit_path),
    " | CDS Transcript"
  )
  
  transcript_out_path <- file.path(
    transcript_output_dir_path,
    paste0(gene_name, "_", genome_name, "_transcript.fasta")
  )
  
  writeLines(
    c(transcript_header, as.character(full_transcript)),
    transcript_out_path
  )
  
  message("Generated protein sequence: ", basename(protein_out_path))
  return(protein_out_path)
}

# Sequence Alignment Function ---------------------------------------------
#' Align translated protein sequences with original query sequences
#'
#' @param protein_file Path to the translated protein sequence file
#' @param query_file Path to the original query protein file
#' @param output_dir Directory for output files
#' @return Path to the alignment file
align_multiple_sequences <- function(reference_file, 
                                     comparison_file1, 
                                     comparison_file2,
                                     output_dir) {
  
  # Create directory for alignments
  align_dir <- file.path(output_dir, "sequence_alignments")
  if (!dir.exists(align_dir)) {
    dir.create(align_dir, recursive = TRUE)
  }
  
  # Generate descriptive output file name
  reference_basename <- tools::file_path_sans_ext(basename(reference_file))
  comp1_basename <- tools::file_path_sans_ext(basename(comparison_file1))
  comp2_basename <- tools::file_path_sans_ext(basename(comparison_file2))
  
  alignment_name <- paste0(reference_basename, "_vs_", comp1_basename, "_vs_", comp2_basename, ".aln")
  alignment_file <- file.path(align_dir, alignment_name)
  
  message("Aligning sequences: ", reference_basename, " with ", comp1_basename, " and ", comp2_basename)
  
  # Read sequences
  tryCatch({
    reference_seq <- readAAStringSet(reference_file)
    comp1_seq <- readAAStringSet(comparison_file1)
    comp2_seq <- readAAStringSet(comparison_file2)
    
    # For simplicity, we'll use just the first sequence from each file
    # If you need specific sequence selection logic, you can adapt that part from the original function
    if (length(reference_seq) > 1) reference_seq <- reference_seq[1]
    if (length(comp1_seq) > 1) comp1_seq <- comp1_seq[1]
    if (length(comp2_seq) > 1) comp2_seq <- comp2_seq[1]
    
    # Add descriptions for clarity
    names(reference_seq) <- paste0(names(reference_seq), " [REFERENCE]")
    names(comp1_seq) <- paste0(names(comp1_seq), " [COMPARISON1]")
    names(comp2_seq) <- paste0(names(comp2_seq), " [COMPARISON2]")
    
    # Combine sequences
    combined_seqs <- c(reference_seq, comp1_seq, comp2_seq)
    
    # Perform alignment using MUSCLE algorithm
    library(msa)
    alignment <- msa(combined_seqs, method = "Muscle")
    
    # Write alignment as FASTA file
    writeXStringSet(
      as(alignment, "AAStringSet"),
      filepath = alignment_file,
      format = "fasta"
    )
    
    # Check if alignment file was created successfully
    if (!file.exists(alignment_file) || file.size(alignment_file) == 0) {
      warning("Alignment failed for multiple sequences")
      return(NA)
    }
    
    message("Alignment completed: ", basename(alignment_file))
    return(alignment_file)
    
  }, error = function(e) {
    warning("Error performing alignment: ", e$message)
    return(NA)
  })
}


# Main Pipeline Function --------------------------------------------------

#' Main execution function for ProBLASTr pipeline
#'
#' Runs the entire workflow from BLAST searches through protein sequence generation
#'
#' @param genome_dir Directory containing genome files
#' @param query_dir Directory containing query protein files
#' @param gff_dir Directory containing gene annotation GFF files
#' @param output_dir Directory for output files
#' @param evalue E-value threshold for BLAST searches (default: 1e-5)
#' @param generate_alignments Whether to generate sequence alignments (default: TRUE)
#' @return Data frame with complete pipeline results
run_problaster_pipeline <- function(genome_dir,
                                    query_dir,
                                    gff_dir,
                                    output_dir,
                                    evalue = 1e-5,
                                    generate_alignments = TRUE,
                                    reference_file) {
  
  # Initialize error logging
  setup_error_logging()
  
  # Check dependencies
  check_dependencies()
  
  # Start timestamp for the run
  start_time <- Sys.time()
  message("Starting ProBLASTr pipeline at: ", start_time)
  message("------------------------------------------------------")
  
  # Check required parameters
  if (missing(genome_dir) || 
      missing(query_dir) || 
      missing(gff_dir) || 
      missing(output_dir)) {
    stop("Required parameter missing! Please provide all required directories.")
  }
  
  # Create main output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  }
  
  # Read and validate index file mapping genomes to GFF files
  index_file_path <- "index.csv"
  if (!file.exists(index_file_path)) {
    stop("Required index file 'index.csv' does not exist. This file should map genome files to GFF files.")
  }
  
  message("Reading genome-to-GFF mapping from: ", index_file_path)
  genome_to_gff_mapping <- read_csv(file = index_file_path, show_col_types = FALSE) |>
    na.omit()
  
  if (!all(c("gff_file_path", "genome_file_path") %in% colnames(genome_to_gff_mapping))) {
    stop("index.csv must contain columns named 'gff_file_path' and 'genome_file_path'")
  }
  
  message("Found mappings for ", nrow(genome_to_gff_mapping), " genomes")
  
  # Update file paths in mapping to include directory prefixes
  genome_to_gff_mapping$gff_file_path <- file.path(gff_dir, genome_to_gff_mapping$gff_file_path)
  genome_to_gff_mapping$genome_file_path <- file.path(genome_dir, genome_to_gff_mapping$genome_file_path)
  
  # Verify that all mapped files exist
  missing_files <- c(
    genome_to_gff_mapping$genome_file_path[!file.exists(genome_to_gff_mapping$genome_file_path)],
    genome_to_gff_mapping$gff_file_path[!file.exists(genome_to_gff_mapping$gff_file_path)]
  )
  
  if (length(missing_files) > 0) {
    warning("Some files referenced in the index do not exist:", 
            paste("\n -", missing_files), 
            "\nContinuing with available files.")
  }
  
  # Create a lookup mapping from genome names to GFF file paths
  genome_name_to_gff <- setNames(
    genome_to_gff_mapping$gff_file_path,
    tools::file_path_sans_ext(basename(genome_to_gff_mapping$genome_file_path))
  )
  
  # Step 1: Run BLAST searches
  message("\n== STEP 1: Running BLAST searches ==")
  blast_results <- create_blastdb_and_run_tblastn(
    genome_dir = genome_dir, 
    query_dir = query_dir, 
    output_dir = output_dir, 
    evalue = evalue
  )
  
  # Step 2: Process BLAST hits to find corresponding genes
  message("\n== STEP 2: Finding genes corresponding to BLAST hits ==")
  
  # Create directory for gene hits
  gene_hit_dir_path <- file.path(output_dir, "gene_matches")
  if (!dir.exists(gene_hit_dir_path)) {
    dir.create(gene_hit_dir_path, recursive = TRUE)
  }
  
  # Initialize data frame to track gene identification results
  gene_identification_results <- data.frame(
    blast_hit_file_path = character(0),
    target_genome = character(0),
    gff_file_path = character(0),
    gene_match_file_path = character(0),
    stringsAsFactors = FALSE
  )
  
  # Get all BLAST hit files
  blast_hit_files <- list.files(
    path = file.path(output_dir, "blast_search_results"),
    pattern = "blast_hits\\.tsv$",
    full.names = TRUE
  )
  
  message("Processing ", length(blast_hit_files), " BLAST hit files")
  
  # Process each BLAST hit file
  for (blast_hit_path in blast_hit_files) {
    blast_filename <- basename(blast_hit_path)
    
    # Extract target genome name from the filename
    target_genome <- sub("^.*_vs_(.*)_blast_hits\\.tsv$", "\\1", blast_filename)
    message("  Processing hits in genome: ", target_genome)
    
    # Read BLAST hits
    blast_hit_table <- read.table(
      file = blast_hit_path,
      sep = "\t",
      col.names = c(
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"
      )
    ) |>
      mutate(sseqid = as.character(sseqid))
    
    # Find corresponding GFF file
    gff_file_path <- genome_name_to_gff[target_genome]
    
    if (!is.na(gff_file_path) && file.exists(gff_file_path)) {
      message("    Using GFF file: ", basename(gff_file_path))
      
      # Read GFF file
      gff_file_table <- read.table(
        file = gff_file_path,
        sep = "\t",
        col.names = c(
          "sseqid", "source", "type", "gstart", "gend", "score", 
          "strand", "phase", "gname"
        )
      ) |>
        mutate(sseqid = as.character(sseqid)) |>
        filter(type %in% c("mRNA", "transcript"))
      
      # Find genes from BLAST hits
      gene_matches <- find_gene_from_blast_hit(
        blast_hit = blast_hit_table,
        gff_file = gff_file_table
      )
      
      if (!is.null(gene_matches) && nrow(gene_matches) > 0) {
        # Extract query name from blast filename
        query_name <- sub("^(.*)_vs_.*_blast_hits\\.tsv$", "\\1", blast_filename)
        
        # Create descriptive gene match filename
        gene_match_file_path <- file.path(
          gene_hit_dir_path,
          paste0(query_name, "_matches_in_", target_genome, ".csv")
        )
        
        # Save gene match results
        write_csv(gene_matches, file = gene_match_file_path)
        message("    Found ", nrow(gene_matches), " gene matches, saved to: ", 
                basename(gene_match_file_path))
        
        # Update tracking data frame
        gene_identification_results <- rbind(
          gene_identification_results,
          data.frame(
            blast_hit_file_path = blast_hit_path,
            target_genome = target_genome,
            gff_file_path = gff_file_path,
            gene_match_file_path = gene_match_file_path,
            stringsAsFactors = FALSE
          )
        )
      } else {
        message("    No gene matches found for BLAST hits in: ", target_genome)
      }
    } else {
      warning("No matching GFF file found for genome: ", target_genome)
    }
  }
  
  # Save gene identification results
  gene_mapping_file <- file.path(output_dir, "gene_identification_mapping.csv")
  write_csv(gene_identification_results, file = gene_mapping_file)
  message("Gene identification results saved to: ", gene_mapping_file)
  
  # Step 3: Extract and translate gene sequences to proteins
  message("\n== STEP 3: Extracting and translating gene sequences ==")
  
  # Combine our results with the original genome mapping
  pipeline_results <- gene_identification_results |>
    left_join(
      genome_to_gff_mapping, 
      by = "gff_file_path"
    ) |>
    # remove target_genome, so it won't appear twice in the next join
    dplyr::select(-target_genome)
  
  # Also add query information
  pipeline_results <- left_join(pipeline_results,
                                blast_results,
                                by = "blast_hit_file_path")
  
  
  # Initialize data frame to track protein generation
  protein_generation_results <- data.frame(
    gene_match_file_path = character(0),
    transcript_id = character(0),
    protein_sequence_path = character(0),
    stringsAsFactors = FALSE
  )
  
  # Process each gene match file
  message("Processing ", nrow(pipeline_results), " gene match files for protein generation")
  
  for (i in seq_len(nrow(pipeline_results))) {
    gene_match_path <- pipeline_results$gene_match_file_path[i]
    gff_file_path <- pipeline_results$gff_file_path[i]
    genome_file_path <- pipeline_results$genome_file_path[i]
    
    message("  Processing gene matches from: ", basename(gene_match_path))
    
    # Read gene match data
    gene_match_data <- read_csv(gene_match_path, show_col_types = FALSE)
    
    # Extract transcript IDs
    transcript_id_raw <- gene_match_data$gname
    transcript_ids <- unique(stri_extract_first(transcript_id_raw, regex = "(?<=ID=)[^;]+"))
    
    message("    Found ", length(transcript_ids), " unique transcript(s) to process")
    
    
    
    # Process each transcript
    for (transcript_id in transcript_ids) {
      message("    Processing transcript: ", transcript_id)
      
      tryCatch({
        # Extract and translate sequence
        protein_sequence_path <- gene_seq_to_prot(
          genome_file = genome_file_path,
          gff_file = gff_file_path,
          transcript_id = transcript_id,
          output_dir = output_dir,
          gene_hit_path = gene_match_path
        )
        
        # Record protein generation result
        protein_generation_results <- rbind(
          protein_generation_results,
          data.frame(
            gene_match_file_path = gene_match_path,
            transcript_id = transcript_id,
            protein_sequence_path = protein_sequence_path,
            stringsAsFactors = FALSE
          )
        )
      }, error = function(e) {
        warning("Error processing transcript ", transcript_id, ": ", e$message)
      })
    }
  }
  
  # Save protein generation results
  protein_mapping_file <- file.path(output_dir, "protein_sequence_mapping.csv")
  write_csv(protein_generation_results, file = protein_mapping_file)
  message("Protein generation results saved to: ", protein_mapping_file)
  
  # Combine all results into a single comprehensive mapping
  final_results <- pipeline_results |>
    left_join(protein_generation_results, by = "gene_match_file_path") |>
    mutate(
      query_name = ifelse(is.na(query_name), 
                          sub("^(.*)_matches_in_.*\\.csv$", "\\1", basename(gene_match_file_path)),
                          query_name)
    )
  
  # Create a meaningful final mapping file
  final_mapping_file <- file.path(output_dir, "problaster_complete_pipeline_results.csv")
  write_csv(final_results, file = final_mapping_file)
  message("Complete pipeline results saved to: ", final_mapping_file)
  
  # Step 4: Generate sequence alignments (if requested)
  if (generate_alignments && nrow(protein_generation_results) > 0) {
    message("\n== STEP 4: Generating sequence alignments ==")
    
    alignment_results <- data.frame(
      protein_sequence_path = character(0),
      query_path = character(0),
      alignment_file_path = character(0),
      stringsAsFactors = FALSE
    )
    
    # Match each protein sequence with its original query
    for (i in seq_len(nrow(final_results))) {
      if (!is.na(final_results$protein_sequence_path[i]) && 
          !is.na(final_results$query_path[i])) {
        
        protein_file <- final_results$protein_sequence_path[i]
        query_file <- final_results$query_path[i]
        
        message("  Aligning: ", basename(protein_file), " with ", basename(query_file))
        
        # Generate the alignment
        alignment_file <- align_multiple_sequences(reference_file = reference_file, 
                                 comparison_file1 = protein_file,
                                 comparison_file2 = query_file,
                                 output_dir = output_dir)
        
        
        # Record alignment result
        if (!is.na(alignment_file)) {
          alignment_results <- rbind(
            alignment_results,
            data.frame(
              protein_sequence_path = protein_file,
              query_path = query_file,
              alignment_file_path = alignment_file,
              stringsAsFactors = FALSE
            )
          )
        }
      }
    }
    
    # Save alignment results
    if (nrow(alignment_results) > 0) {
      alignment_mapping_file <- file.path(output_dir, "sequence_alignment_mapping.csv")
      write_csv(alignment_results, file = alignment_mapping_file)
      message("Sequence alignment results saved to: ", alignment_mapping_file)
      
      # Add alignment information to final results
      final_results <- final_results |>
        left_join(
          alignment_results |> 
            dplyr::select(protein_sequence_path, alignment_file_path),
          by = "protein_sequence_path"
        )
      
      # Update the complete results file
      write_csv(final_results, file = final_mapping_file)
    } else {
      message("No alignments were generated.")
    }
  }
  
  # Generate a summary report
  pipeline_summary <- data.frame(
    stage = c("BLAST searches", "Gene matches", "Protein sequences", "Alignments"),
    count = c(
      nrow(blast_results),
      nrow(gene_identification_results),
      nrow(protein_generation_results),
      if(generate_alignments) nrow(alignment_results) else NA
    )
  )
  
  summary_file <- file.path(output_dir, "pipeline_summary.csv")
  write_csv(pipeline_summary, file = summary_file)
  
  # Calculate execution time
  end_time <- Sys.time()
  execution_time <- difftime(end_time, start_time, units = "mins")
  
  message("\n== Pipeline Execution Summary ==")
  message("Total execution time: ", round(as.numeric(execution_time), 2), " minutes")
  message("BLAST searches completed: ", pipeline_summary$count[1])
  message("Gene matches found: ", pipeline_summary$count[2])
  message("Protein sequences generated: ", pipeline_summary$count[3])
  if (generate_alignments) {
    message("Sequence alignments created: ", pipeline_summary$count[4])
  }
  message("All results saved to: ", output_dir)
  message("Complete pipeline mapping: ", final_mapping_file)
  
  # Return the final results dataframe
  return(final_results)
}


# Utility Functions -------------------------------------------------------
#' Filters gene hits based on position and score 
#'
#' This function filters gene hits to identify the most significant matches 
#' by including only those above a certain bitscore and keeps only one gene
#' of those which are close together.
#'
#' @param blast_hits DataFrame containing BLAST hits
#' @param max_distance Maximum distance between hit regions to be considered separate
#' @param min_bitscore Minimum bitscore threshold for filtering
#' @return Filtered DataFrame with the most significant hits
filter_gene_hits <- function(blast_hits, max_distance = 10000, min_bitscore = 30) {
  if (nrow(blast_hits) == 0) {
    return(blast_hits)
  }
  
  # Filter by bitscore
  filtered_hits <- blast_hits |>
    filter(bitscore >= min_bitscore)
  
  if (nrow(filtered_hits) == 0) {
    return(filtered_hits)
  }
  
  # Group by sequence ID and order by position
  ordered_hits <- filtered_hits |>
    group_by(sseqid) |>
    arrange(sseqid, pmin(sstart, send)) |>
    mutate(
      hit_start = pmin(sstart, send),
      hit_end = pmax(sstart, send)
    )
  
  # Identify clusters of hits that are close to each other
  clustered_hits <- ordered_hits |>
    mutate(
      cluster = cumsum(c(1, diff(hit_start) > max_distance))
    )
  
  # Get the best hit from each cluster
  best_hits <- clustered_hits |>
    group_by(sseqid, cluster) |>
    slice_max(order_by = bitscore, n = 1) |>
    ungroup() |>
    dplyr::select(-cluster)
  
  return(best_hits)
}


# Example call ------------------------------------------------------------
# Example of how to run the complete pipeline
if (T) {  # Set to TRUE to execute when sourcing this file
  results <- run_problaster_pipeline(
    genome_dir = "genomes/test_genomes",
    query_dir = "meiotic_genes_protein_fasta/test_meiotic_prot", 
    gff_dir = "gff_files",
    output_dir = "out_2025_02_28v8_SPO11_O_hap1_w_A_tha_aln",
    evalue = 1e-5,
    generate_alignments = TRUE,
    reference_file = "A_tha/A_tha_SPO11_1.fasta"
  )
}
