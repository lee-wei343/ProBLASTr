#!/usr/bin/env Rscript
#
# Script Name: ProBLASTr.R
# Description: A pipeline for protein BLAST analysis to identify homologous genes
#              across genomic datasets
#
# Copyright (C) [2025] 
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Version: 1.2
# Last Updated: 2025-03-05
#
# Author:
# Lee Weinand
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
  library(pwalign)
})

check_dependencies <- function() {
  required_packages <- c(
    "tidyverse", "GenomicRanges", "Biostrings", "GenomicFeatures",
    "rtracklayer", "Rsamtools", "stringi", "msa", "pwalign"
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
setup_logging <- function() {
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
    
    # Store the database path for indexing
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
  
  # Run tBLASTn for each query against each database
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
  
  # Only perform the join if matching GFF entries were found
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
                             gene_hit_path,
                             seqid = NULL) {
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
  if (!is.null(seqid)) {
    protein_header <- paste0(
      protein_header, " | seqid: ", seqid
    )
  }
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
  if (!is.null(seqid)) {
    transcript_header <- paste0(
      transcript_header, " | seqid: ", seqid
    )
  }
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

# Sequence Alignment Functions ---------------------------------------------
#' Pairwise align the query protein sequence against the subject protein sequence and 
#' score the percentage identity and the BLOSUM65 score
pairwise_align_protein_sequences_and_score <- function(query_file,
                                                       protein_files,
                                                       output_dir,
                                                       query_name,
                                                       gap_opening = -10,
                                                       gap_extension = -0.5) {
  # Create directory for alignments
  align_dir <- file.path(output_dir, "pairwise_alignments")
  if (!dir.exists(align_dir)) {
    dir.create(align_dir, recursive = TRUE)
  }
  
  # Generate descriptive output file name
  alignment_name <- paste0(query_name, "_pairwise_alignment.aln")
  alignment_file <- file.path(align_dir, alignment_name)
  
  message("Creating pairwise sequence alignment for query: ", query_name)
  message("  Including ",
          length(protein_files),
          " sequences from different genomes")
  
  # Read sequences
  tryCatch({
    # Read query sequence
    query_seq <- readAAStringSet(query_file)
    if (length(query_seq) > 1) {
      # If multiple sequences in the query file, just use the first one
      query_seq <- query_seq[1]
    }
    
    # Read protein sequence
    for (protein_file in protein_files) {
      if (file.exists(protein_file)) {
        protein_seq <- readAAStringSet(protein_file)
        # Add genome identifier from filename to make sequence names unique
        genome_id <- sub(".*_(.*?)_protein\\.fasta$",
                         "\\1",
                         basename(protein_file))
        names(protein_seq) <- paste0(names(protein_seq), " [", genome_id, "]")
        
        
        # Rename query sequence for clarity
        names(query_seq) <- paste0(names(query_seq), " [QUERY]")
        
        # load alignment data
        data("BLOSUM62")
        
        # Pairwise align query and subject protein sequence
        palign <- pairwiseAlignment(
          query_seq,
          protein_seq,
          substitutionMatrix = BLOSUM62,
          # Using BLOSUM62 for homologous proteins
          gapOpening = gap_opening,
          # Gap opening penalty
          gapExtension = gap_extension,
          # Gap extension penalty
          scoreOnly = FALSE,
          # Return full alignment, not just score
          type = "global"
          # Use global alignment for homologous proteins
        )
        
        # Score alignment
        score <- score(palign)
        pid <- pid(palign)
        
        # Write alignment as FASTA file
        writeXStringSet(as(palign, "AAStringSet"),
                        filepath = alignment_file,
                        format = "fasta")
      }
    }
    
    message("  Alignment completed: ", basename(alignment_file))
    message("    Pairwise alignment score: ", score)
    message("    Pairwise alignment identity: ", round(pid))
    return(list(pscore = score, pid = pid))
    
  }, error = function(e) {
    warning("Error performing pairwise sequence alignment: ", e$message)
    return(NA)
  })
}

#' Align multiple protein sequences from different genomes with the query sequence
#'
#' @param query_file Path to the original query protein file
#' @param protein_files Vector of paths to the translated protein sequences from different genomes
#' @param output_dir Directory for output files
#' @param query_name Name of the query for file naming
#' @return Path to the alignment file
align_multiple_sequences <- function(query_file, 
                                     protein_files, 
                                     output_dir,
                                     query_name) {
  
  # Create directory for alignments
  align_dir <- file.path(output_dir, "multiple_alignments")
  if (!dir.exists(align_dir)) {
    dir.create(align_dir, recursive = TRUE)
  }
  
  # Generate descriptive output file name
  alignment_name <- paste0(query_name, "_multiple_alignment.aln")
  alignment_file <- file.path(align_dir, alignment_name)
  
  message("Creating multiple sequence alignment for query: ", query_name)
  message("  Including ", length(protein_files), " sequences from different genomes")
  
  # Read sequences
  tryCatch({
    # Read query sequence
    query_seq <- readAAStringSet(query_file)
    if (length(query_seq) > 1) {
      # If multiple sequences in the query file, just use the first one
      query_seq <- query_seq[1]
    }
    
    # Read all protein sequences
    all_seqs <- query_seq
    for (protein_file in protein_files) {
      if (file.exists(protein_file)) {
        protein_seq <- readAAStringSet(protein_file)
        # Add genome identifier from filename to make sequence names unique
        genome_id <- sub(".*_(.*?)_protein\\.fasta$", "\\1", basename(protein_file))
        names(protein_seq) <- paste0(names(protein_seq), " [", genome_id, "]")
        all_seqs <- c(all_seqs, protein_seq)
        
        # Rename query sequence for clarity
        names(all_seqs)[1] <- paste0(names(all_seqs)[1], " [QUERY]")
        
        # Perform alignment using MUSCLE algorithm
        alignment <- msa(all_seqs, method = "Muscle")
        
        # score alignment
        
        
        # Write alignment as FASTA file
        writeXStringSet(
          as(alignment, "AAStringSet"),
          filepath = alignment_file,
          format = "fasta"
        )
        
        
        
        
      }
    }
    
    
    
    message("  Alignment completed: ", basename(alignment_file))
    return(alignment_file)
    
  }, error = function(e) {
    warning("Error performing multiple sequence alignment: ", e$message)
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
# Main Pipeline Function - modified to track qseqid values
run_problaster_pipeline <- function(genome_dir,
                                    query_dir,
                                    gff_dir,
                                    output_dir,
                                    evalue = 1e-5,
                                    generate_alignments = TRUE,
                                    min_score_threshold,
                                    min_pid_threshold,
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
  message("\n== STEP 1: Running tBLASTn searches ==")
  blast_results <- create_blastdb_and_run_tblastn(
    genome_dir = genome_dir, 
    query_dir = query_dir, 
    output_dir = output_dir, 
    evalue = evalue
  )
  
  # Step 2: Process BLAST hits to find corresponding genes
  message("\n== STEP 2: Finding annotated genes corresponding to tBLASTn hits ==")
  
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
      mutate(sseqid = as.character(sseqid),
             qseqid = as.character(qseqid))  # Ensure qseqid is character type
    
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
  
  # Add query information
  pipeline_results <- left_join(pipeline_results,
                                blast_results,
                                by = "blast_hit_file_path")
  
  # Initialize data frame to track protein generation
  protein_generation_results <- data.frame(
    gene_match_file_path = character(0),
    transcript_id = character(0),
    protein_sequence_path = character(0),
    qseqid = character(0),  # Query sequence ID
    seqid = character(0),   # Added column to store seqid (subject sequence ID)
    gstart = character(0),
    gend = character(0),
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
    
    # Extract transcript IDs and their corresponding qseqid's and seqid's from column gname
    transcript_info <- gene_match_data %>%
      group_by(gname) %>%
      summarize(
        transcript_id = stri_extract_first(gname[1], regex = "(?<=ID=)[^;]+"),
        qseqid = qseqid[1],  # Extract the first Qery Sequence ID
        seqid = sseqid[1],    # Extract the first Subject Sequence ID 
        gstart = gstart[1],  # added
        gend = gend[1]  # added
      ) %>%
      filter(!is.na(transcript_id))
    
    message("    Found ", nrow(transcript_info), " unique transcript(s) to process")
    
    # Process each transcript
    for (j in seq_len(nrow(transcript_info))) {
      transcript_id <- transcript_info$transcript_id[j]
      current_qseqid <- transcript_info$qseqid[j]
      current_seqid <- transcript_info$seqid[j]  # Get seqid
      current_gstart <- transcript_info$gstart[j]
      current_gend <- transcript_info$gend[j]
      
      message("    Processing transcript: ", transcript_id, 
              " (qseqid: ", current_qseqid, ", seqid: ", current_seqid, ")")
      
      tryCatch({
        # Extract and translate sequence
        protein_sequence_path <- gene_seq_to_prot(
          genome_file = genome_file_path,
          gff_file = gff_file_path,
          transcript_id = transcript_id,
          output_dir = output_dir,
          gene_hit_path = gene_match_path,
          seqid = current_seqid
        )
        
        # Record protein generation result with qseqid, seqid, gstart and gend
        protein_generation_results <- rbind(
          protein_generation_results,
          data.frame(
            gene_match_file_path = gene_match_path,
            transcript_id = transcript_id,
            protein_sequence_path = protein_sequence_path,
            qseqid = current_qseqid,
            seqid = current_seqid,
            gstart = current_gstart,
            gend = current_gend,
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
  
  # Step 4: Filter sequences based on pairwise alignment score and generate alignments
  alignment_results <- data.frame()
  filtered_protein_results <- data.frame()
  
  if (generate_alignments && nrow(protein_generation_results) > 0) {
    message("\n== STEP 4: Filtering sequences by alignment score and generating alignments ==")
    
    alignment_results <- data.frame(
      query_name = character(0),
      query_path = character(0),
      alignment_file_path = character(0),
      num_sequences = integer(0),
      stringsAsFactors = FALSE
    )
    
    filtered_protein_results <- data.frame(
      query_name = character(0),
      protein_sequence_path = character(0),
      pscore = numeric(0),
      pid = numeric(0),
      passed_filter = logical(0),
      stringsAsFactors = FALSE
    )
    
    # Group protein sequences by query
    protein_by_query <- pipeline_results %>%
      left_join(protein_generation_results, by = "gene_match_file_path") %>%
      filter(!is.na(protein_sequence_path)) %>%
      group_by(query_name, query_path) %>%
      summarise(
        protein_files = list(protein_sequence_path),
        num_sequences = n(),
        .groups = "drop"
      )
    
    # Assess and filter protein sequences for each query
    for (i in seq_len(nrow(protein_by_query))) {
      query_name <- protein_by_query$query_name[i]
      query_path <- protein_by_query$query_path[i]
      protein_files <- unlist(protein_by_query$protein_files[i])
      
      message("Evaluating and filtering sequences for query: ", query_name)
      message("  Analyzing ", length(protein_files), " sequences from different genomes")
      
      # Score each protein sequence against the query
      filtered_proteins <- c()
      scores_data <- data.frame(
        protein_sequence_path = character(0),
        pscore = numeric(0),
        pid = numeric(0),
        passed_filter = logical(0),
        stringsAsFactors = FALSE
      )
      
      for (protein_file in protein_files) {
        message("  Scoring: ", basename(protein_file))
        
        # Get alignment score and percent identity
        alignment_result <- pairwise_align_protein_sequences_and_score(
          query_file = query_path,
          protein_files = c(protein_file),
          output_dir = output_dir,
          query_name = paste0(query_name, "_", basename(protein_file))
        )
        
        # Check if alignment was successful
        if (!is.na(alignment_result$pscore)) {
          # Apply filtering criteria
          passed_filter <- (alignment_result$pscore >= min_score_threshold && 
                              alignment_result$pid >= min_pid_threshold)
          
          # Add to scores data
          scores_data <- rbind(
            scores_data,
            data.frame(
              protein_sequence_path = protein_file,
              pscore = alignment_result$pscore,
              pid = alignment_result$pid,
              passed_filter = passed_filter,
              stringsAsFactors = FALSE
            )
          )
          
          # Create score vs pid scatter plot
          if (nrow(scores_data) > 0) {
            plot <- ggplot(data = scores_data,
                           mapping = aes(x = pscore, y = pid)) +
              geom_point() +
              labs(title = paste("Alignment Scores for", query_name),
                   x = "BLOSUM Score",
                   y = "Percent Identity")
            
            # Create plot directory if it doesn't exist
            plot_dir <- file.path(output_dir, "alignment_plots")
            if (!dir.exists(plot_dir)) {
              dir.create(plot_dir, recursive = TRUE)
            }
            
            # Save plot with more descriptive filename
            ggsave(
              filename = paste0(query_name, "_alignment_scores.png"), 
              plot = plot, 
              path = plot_dir,
              width = 8,
              height = 6
            )
          }
          
          if (passed_filter) {
            filtered_proteins <- c(filtered_proteins, protein_file)
            message("    Passed filter: Score=", round(alignment_result$pscore, 1), 
                    ", PID=", round(alignment_result$pid, 1), "%")
          } else {
            message("    Failed filter: Score=", round(alignment_result$pscore, 1), 
                    ", PID=", round(alignment_result$pid, 1), "% (below thresholds)")
          }
        } else {
          message("    Alignment failed, excluding from further analysis")
        }
      }
      
      # Save filtering results for this query
      filtered_protein_results <- rbind(
        filtered_protein_results,
        data.frame(
          query_name = query_name,
          scores_data,
          stringsAsFactors = FALSE
        )
      )
      
      # Only proceed with multiple sequence alignment if enough sequences passed filtering
      if (length(filtered_proteins) >= 2) {
        message("  Proceeding with multiple sequence alignment using ", 
                length(filtered_proteins), " filtered sequences")
        
        # Generate multiple sequence alignment with filtered sequences
        alignment_file <- align_multiple_sequences(
          query_file = query_path,
          protein_files = filtered_proteins,
          output_dir = output_dir,
          query_name = paste0(query_name, "_filtered")
        )
        
        # Record alignment result
        if (!is.na(alignment_file)) {
          alignment_results <- rbind(
            alignment_results,
            data.frame(
              query_name = query_name,
              query_path = query_path,
              alignment_file_path = alignment_file,
              num_sequences = length(filtered_proteins),
              stringsAsFactors = FALSE
            )
          )
        }
      } else {
        message("  Not enough sequences passed filtering criteria for ", query_name, 
                ". Need at least 2 sequences for multiple alignment.")
      }
    }
    
    # Save filtering results
    filtering_results_file <- file.path(output_dir, "sequence_filtering_results.csv")
    write_csv(filtered_protein_results, file = filtering_results_file)
    message("Sequence filtering results saved to: ", filtering_results_file)
    
    # Save alignment results
    if (nrow(alignment_results) > 0) {
      alignment_mapping_file <- file.path(output_dir, "multiple_alignment_mapping.csv")
      write_csv(alignment_results, file = alignment_mapping_file)
      message("Multiple sequence alignment results saved to: ", alignment_mapping_file)
    } else {
      message("No multiple sequence alignments were generated after filtering.")
    }
  }
  
  # Combine all results into a single comprehensive mapping - MOVED AFTER STEP 4
  # Start with basic joining of results
  intermediate_results <- pipeline_results %>%
    left_join(protein_generation_results, by = "gene_match_file_path") %>%
    mutate(
      query_name = ifelse(is.na(query_name), 
                          sub("^(.*)_matches_in_.*\\.csv$", "\\1", basename(gene_match_file_path)),
                          query_name)
    )
  
  # Add filtering information if available
  if (nrow(filtered_protein_results) > 0) {
    intermediate_results <- intermediate_results %>%
      left_join(
        filtered_protein_results %>% 
          select(query_name, protein_sequence_path, pscore, pid, passed_filter),
        by = c("query_name", "protein_sequence_path")
      )
    
    # Filter to keep only hits that passed the threshold
    final_results <- intermediate_results %>%
      filter(is.na(passed_filter) | passed_filter == TRUE)
    
    message("Filtering final results: kept ", nrow(final_results), 
            " out of ", nrow(intermediate_results), " hits that passed threshold criteria")
  } else {
    # If no filtering was done, use all results
    final_results <- intermediate_results
    message("No alignment filtering was performed, including all hits in final results")
  }
  
  # Add alignment information if available
  if (nrow(alignment_results) > 0) {
    final_results <- final_results %>%
      left_join(alignment_results %>% select(query_name, alignment_file_path), by = "query_name")
  }
  
  # Create a meaningful final mapping file
  final_mapping_file <- file.path(output_dir, "problaster_complete_pipeline_results.csv")
  write_csv(final_results, file = final_mapping_file)
  message("Complete pipeline results saved to: ", final_mapping_file)
  
  # Generate a summary report
  pipeline_summary <- data.frame(
    stage = c("BLAST searches", "Gene matches", "Protein sequences", 
              "Sequences passing filters", "Final alignments"),
    count = c(
      nrow(blast_results),
      nrow(gene_identification_results),
      nrow(protein_generation_results),
      if(generate_alignments) sum(filtered_protein_results$passed_filter) else NA,
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
    message("Sequences passing alignment filters: ", pipeline_summary$count[4])
    message("Multiple sequence alignments created: ", pipeline_summary$count[5])
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
#' by including only those above a certain bitscore.
#'
#' @param blast_hits DataFrame containing BLAST hits
#' @param min_bitscore Minimum bitscore threshold for filtering
#' @return Filtered DataFrame with the most significant hits
filter_gene_hits <- function(blast_hits, min_bitscore = 30) {
  if (nrow(blast_hits) == 0) {
    return(blast_hits)
  }
  
  # Filter by bitscore
  filtered_hits <- blast_hits |>
    filter(bitscore >= min_bitscore)
  
  return(filtered_hits)
}


# Example call ------------------------------------------------------------
# Example of how to run the complete pipeline
if (T) {  # Set to TRUE to execute when sourcing this file
  results <- run_problaster_pipeline(
    genome_dir = "genomes/test_genomes",
    query_dir = "meiotic_genes_protein_fasta/test_meiotic_prot", 
    gff_dir = "gff_files",
    output_dir = "out_2025_03_05_SPO11_test",
    evalue = 1e-5,
    generate_alignments = T,
    min_score_threshold = 300, # Minimum acceptable BLOSUM score
    min_pid_threshold = 25, # Minimum percent identity (0-100)
    reference_file = "A_tha/A_tha_SPO11_1.fasta")}
