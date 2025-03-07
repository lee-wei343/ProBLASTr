#!/usr/bin/env Rscript
#
# Script Name: ProBLASTr.R
# Description: A pipeline for protein BLAST analysis to identify homologous genes
#              across genomic datasets, with integrated phylogenetic tree building
#
# Copyright (C) [2025]
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Version: 1.5
# Last Updated: 2025-03-07
#
# Author:
# Lee Weinand
#
# Dependencies:
# R version 4.4.1 or higher
# Required packages: tidyverse, GenomicRanges, Biostrings, GenomicFeatures,
#                    rtracklayer, Rsamtools, stringi, msa, pwalign
# For tree building: ape, phangorn, DECIPHER, ggtree
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
  library(ape)
  library(phangorn)
  library(DECIPHER)
  library(ggtree)
  library(ggplot2)
  library(seqinr)
})

check_dependencies <- function() {
  required_packages <- c(
    "tidyverse",
    "GenomicRanges",
    "Biostrings",
    "GenomicFeatures",
    "rtracklayer",
    "Rsamtools",
    "stringi",
    "msa",
    "pwalign",
    "ape",
    "phangorn",
    "DECIPHER",
    "ggtree",
    "ggplot2",
    "seqinr"
  )
  
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
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
    stop("No genome files found matching the pattern: ",
         gen_file_pattern)
  }
  
  message("Creating BLAST databases from ",
          length(genome_files),
          " genome files...")
  db_paths <- character(length(genome_files))
  
  # Create BLAST database for each genome file
  for (i in seq_along(genome_files)) {
    genome_file <- genome_files[i]
    genome_name <- tools::file_path_sans_ext(basename(genome_file))
    db_name <- file.path(blast_database_dir, paste0(genome_name, "_blast_db"))
    
    message("  Creating database for: ", genome_name)
    
    system2("makeblastdb",
            args = c("-in", genome_file, "-dbtype", dbtype, "-out", db_name))
    
    # Store the database path for indexing
    db_paths[i] <- db_name
  }
  
  # List and validate query files
  query_files <- list.files(query_dir, pattern = qur_file_pattern, full.names = TRUE)
  if (length(query_files) == 0) {
    stop("No query files found matching the pattern: ",
         qur_file_pattern)
  }
  
  message("Running BLAST searches with ",
          length(query_files),
          " query files...")
  
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
          "-query",
          query,
          "-db",
          db,
          "-out",
          output_file,
          "-evalue",
          evalue,
          "-outfmt",
          "\"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\""
        )
      )
      
      # Handle empty files
      if (file.exists(output_file) && file.size(output_file) == 0) {
        message("    No hits found for ", query_name, " in ", genome_name)
        file.remove(output_file)
      } else if (file.exists(output_file) &&
                 file.size(output_file) > 0) {
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
  
  message("BLAST searches completed. Found hits in ",
          nrow(blast_results),
          " searches.")
  
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
        dplyr::filter(gend >= blast_sstart & gstart <= blast_send)
    } else {
      # Reverse strand - find GFF entries that overlap the blast hit
      gff_file_in_range <- gff_file |>
        dplyr::filter(gend >= blast_send & gstart <= blast_sstart)
    }
    
    # Append the filtered GFF entries to the results
    if (nrow(gff_file_in_range) > 0) {
      df <- rbind(df, gff_file_in_range)
    }
  }
  
  # Only perform the join if matching GFF entries were found
  if (nrow(df) > 0) {
    result <- inner_join(blast_hit, df, by = "sseqid", relationship = "many-to-many")
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
#' @param seqid Subject sequence ID (optional)
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
  
  protein <- Biostrings::translate(
    full_transcript,
    no.init.codon = FALSE,                # Translate regardless of start codon
    if.fuzzy.codon = "solve"              
    # Fuzzy codons that can be translated non ambiguously to an amino acid or 
    # to * (stop codon) will be translated. Ambiguous fuzzy codons will be 
    # translated to X.
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
  
  # Make sure it's only a single prot 
  protein_file <- protein_files[1]
  
  # Generate descriptive output file name
  protein_basename <- basename(protein_file)
  alignment_name <- paste0(query_name, "_", protein_basename, "_pairwise.aln")
  alignment_file <- file.path(align_dir, alignment_name)
  
  message("Creating pairwise sequence alignment for query: ", query_name)
  message("  Comparing with protein: ", protein_basename)
  
  # Read sequences
  tryCatch({
    # Read query sequence
    query_seq <- readAAStringSet(query_file)
    if (length(query_seq) > 1) {
      # If multiple sequences in the query file, just use the first one
      query_seq <- query_seq[1]
    }
    
    # Read protein sequence
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
      
      message("  Alignment completed: ", basename(alignment_file))
      message("    Pairwise alignment score: ", round(score, 1))
      message("    Pairwise alignment identity: ", round(pid, 1), "%")
      return(list(pscore = score, pid = pid))
    } else {
      warning("Protein file not found: ", protein_file)
      return(list(pscore = NA, pid = NA))
    }
  }, error = function(e) {
    warning("Error performing pairwise sequence alignment: ", e$message)
    return(list(pscore = NA, pid = NA))
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
    
    # Read all protein sequences
    all_seqs <- query_seq
    for (protein_file in protein_files) {
      if (file.exists(protein_file)) {
        protein_seq <- readAAStringSet(protein_file)
        # Add genome identifier from filename to make sequence names unique
        genome_id <- sub(".*_(.*?)_protein\\.fasta$",
                         "\\1",
                         basename(protein_file))
        names(protein_seq) <- paste0(names(protein_seq), " [", genome_id, "]")
        all_seqs <- c(all_seqs, protein_seq)
      }
    }
    
    # Rename query sequence for clarity
    names(all_seqs)[1] <- paste0(names(all_seqs)[1], " [QUERY]")
    
    # Perform alignment using MUSCLE algorithm
    alignment <- msa(all_seqs, method = "Muscle")
    
    # Write alignment as FASTA file
    writeXStringSet(
      as(alignment, "AAStringSet"),
      filepath = alignment_file,
      format = "fasta"
    )
    
    message("  Alignment completed: ", basename(alignment_file))
    return(alignment_file)
    
  }, error = function(e) {
    warning("Error performing multiple sequence alignment: ", e$message)
    return(NA)
  })
}


# Phylogenetic Tree function ----------------------------------------------
create_phylogenetic_trees <- function(df, min_sequences = 2) {
  # Initialize list to store tree results
  tree_results <- list()
  
  # Filter for rows that passed the filter
  df_passed <- df |> dplyr::filter(passed_filter == TRUE)
  
  # Get distinct query names
  distinct_queries <- df_passed |> 
    dplyr::distinct(query_name) |> 
    pull(query_name)
  
  message("Found ", length(distinct_queries), " distinct queries")
  
  # Process each query separately
  for (query_name in distinct_queries) {
    message("Processing query: ", query_name)
    
    # Filter dataframe for current query
    query_df <- df_passed |> dplyr::filter(query_name == !!query_name)
    
    # Extract chromosomes from seqids and get unique values
    chromosomes <- query_df$seqid |>
      stringi::stri_extract(regex = "chr[0-9]{2}") |>
      unique() |>
      na.omit()
    
    message("  Found ", length(chromosomes), " chromosomes for query ", query_name)
    
    # For each chromosome, create a separate tree
    for (chrom in chromosomes) {
      message("  Processing chromosome: ", chrom)
      
      # Filter for sequences from this chromosome
      chrom_df <- query_df |>
        dplyr::filter(stringi::stri_detect(seqid, regex = paste0("^", chrom, ".*")))
      
      # Check if enough sequences for this query-chromosome combination
      if (nrow(chrom_df) < min_sequences) {
        message("    Not enough sequences for ", query_name, " on ", chrom, 
                " (found ", nrow(chrom_df), ", minimum required: ", min_sequences, ")")
        next
      }
      
      # Get protein files for this query-chromosome combination
      protein_files <- chrom_df$protein_sequence_path
      
      # Get original query file
      original_query_file <- unique(chrom_df$query_path)[1]
      
      # Determine output directory from alignment file path
      output_dir <- dirname(dirname(chrom_df$alignment_file_path[1]))
      
      # Create tree key for storing results
      tree_key <- paste0(query_name, "_", chrom)
      
      # Call the create_query_tree function with the chromosome in the name
      tree_result <- create_query_tree(
        query_name = tree_key,
        protein_files = protein_files,
        output_dir = output_dir,
        original_query_file = original_query_file,
        chromosome_name = chrom
      )
      
      if (!is.null(tree_result)) {
        tree_results[[tree_key]] <- tree_result
      }
    }
  }
  
  message("Created ", length(tree_results), " phylogenetic trees")
  return(tree_results)
}



create_query_tree <- function(query_name,
                              protein_files,
                              output_dir,
                              original_query_file = NULL,
                              chromosome_name = NULL) {
  # Create directory for trees
  tree_dir <- file.path(output_dir, "phylogenetic_trees")
  if (!dir.exists(tree_dir)) {
    dir.create(tree_dir, recursive = TRUE)
  }
  
  # Create directory for tree visualizations
  tree_plot_dir <- file.path(output_dir, "tree_plots")
  if (!dir.exists(tree_plot_dir)) {
    dir.create(tree_plot_dir, recursive = TRUE)
  }
  
  # Generate descriptive output file names
  tree_name <- paste0(query_name, "_phylogenetic_tree")
  tree_file <- file.path(tree_dir, paste0(tree_name, ".tre"))
  tree_plot_file <- file.path(tree_plot_dir, paste0(tree_name, ".png"))
  
  message("Creating phylogenetic tree for query: ", query_name)
  message("  Including ", length(protein_files), " sequences from different genomes")
  
  tryCatch({
    # Prepare sequences for alignment
    all_seqs <- Biostrings::AAStringSet()
    
    # Add original query sequence if provided
    if (!is.null(original_query_file) && file.exists(original_query_file)) {
      query_seq <- readAAStringSet(original_query_file)
      if (length(query_seq) > 1) {
        # If multiple sequences in the query file, just use the first one
        query_seq <- query_seq[1]
      }
      # Rename query sequence for clarity
      names(query_seq) <- "[QUERY]"
      all_seqs <- c(all_seqs, query_seq)
    }
    
    # Add all protein sequences
    for (protein_file in protein_files) {
      if (file.exists(protein_file)) {
        protein_seq <- readAAStringSet(protein_file)
        all_seqs <- c(all_seqs, protein_seq)
      }
    }
    
    # Proceed only if there're enough sequences
    if (length(all_seqs) < 2) {
      message("  Not enough sequences to build a tree (minimum 2 required)")
      return(NULL)
    }
    
    # Perform multiple sequence alignment
    message("  Aligning sequences using MUSCLE algorithm...")
    alignment <- msa(all_seqs, method = "Muscle")
    
    # Convert alignment to DNAbin format for tree building
    alignment_seqinr <- msaConvert(alignment, type = "seqinr::alignment")
    alignment_ape <- as.DNAbin(alignment_seqinr)
    
    # Compute distance matrix
    message("  Computing distance matrix...")
    dist_matrix <- dist.alignment(alignment_seqinr, "identity")
    
    # Build tree using Neighbor-Joining
    message("  Building Neighbor-Joining tree...")
    tree <- nj(dist_matrix)
    
    # Root tree at midpoint
    tree <- midpoint(tree)
    
    # Save tree to file
    write.tree(tree, file = tree_file)
    message("  Saved phylogenetic tree to: ", basename(tree_file))
    
    browser()
    # Create tree visualization
    message("  Creating tree visualization...")
    # clean tip label name
    tree$tip.label <- clean_tip_label(tree$tip.label)
    
    plot <- ggtree(tree, layout = "circular", open.angle = 340) + 
      geom_tiplab(aes(angle = angle), color='blue') +
      labs(title = paste0("Phylogenetic Tree for Query ", query_name)) 
      
    
    # Save plot to file
    ggsave(tree_plot_file, plot, width = 14, height = 12, dpi = 300)
    message("  Saved tree visualization to: ", basename(tree_plot_file))
    
    # Return tree results
    return(list(
      tree = tree,
      tree_file = tree_file,
      plot_file = tree_plot_file,
      num_sequences = length(all_seqs)
    ))
    
  }, error = function(e) {
    warning("Error creating phylogenetic tree: ", e$message)
    return(NULL)
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
#' @param generate_trees Whether to generate phylogenetic trees (default: TRUE)
#' @param min_score_threshold Minimum alignment score threshold for filtering
#' @param min_pid_threshold Minimum percent identity threshold for filtering
#' @param min_tree_sequences Minimum number of sequences required for tree building
#' @return Data frame with complete pipeline results
run_problaster_pipeline <- function(genome_dir,
                                    query_dir,
                                    gff_dir,
                                    index_path,
                                    output_dir,
                                    evalue = 1e-5,
                                    generate_alignments = TRUE,
                                    generate_trees = TRUE,
                                    min_score_threshold = 300, # Added default value
                                    min_pid_threshold = 25,    # Added default value
                                    min_tree_sequences = 2) {
  # Initialize error logging
  setup_logging()
  
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
      missing(output_dir) ||
      missing(index_path)) {
    stop("Required parameter missing! Please provide all required directories.")
  }
  
  # Create main output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  }
  
  # Read and validate index file mapping genomes to GFF files
  index_file_path <- index_path
  if (!file.exists(index_file_path)) {
    stop(
      "Required index file 'index.csv' does not exist. This file should map genome files to GFF files."
    )
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
  missing_files <- c(genome_to_gff_mapping$genome_file_path[!file.exists(genome_to_gff_mapping$genome_file_path)],
                     genome_to_gff_mapping$gff_file_path[!file.exists(genome_to_gff_mapping$gff_file_path)])
  
  if (length(missing_files) > 0) {
    warning(
      "Some files referenced in the index do not exist:",
      paste("\n -", missing_files),
      "\nContinuing with available files."
    )
  }
  
  # Create a lookup mapping from genome names to GFF file paths
  genome_name_to_gff <- setNames(genome_to_gff_mapping$gff_file_path,
                                 tools::file_path_sans_ext(basename(
                                   genome_to_gff_mapping$genome_file_path
                                 )))
  
  
  
  ## Step 1: Run BLAST searches 
  message("\n== STEP 1: Running tBLASTn searches ==")
  blast_results <- create_blastdb_and_run_tblastn(
    genome_dir = genome_dir,
    query_dir = query_dir,
    output_dir = output_dir,
    evalue = evalue
  )
  
  
  
  ## Step 2: Process BLAST hits to find corresponding genes 
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
    ) |>
      mutate(sseqid = as.character(sseqid), qseqid = as.character(qseqid))  # Ensure qseqid is character type
    
    # Find corresponding GFF file
    gff_file_path <- genome_name_to_gff[target_genome]
    
    if (!is.na(gff_file_path) && file.exists(gff_file_path)) {
      message("    Using GFF file: ", basename(gff_file_path))
      
      # Read GFF file
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
        dplyr::filter(type %in% c("mRNA", "transcript"))
      
      # Find genes from BLAST hits
      gene_matches <- find_gene_from_blast_hit(blast_hit = blast_hit_table, gff_file = gff_file_table)
      
      if (!is.null(gene_matches) && nrow(gene_matches) > 0) {
        # Extract query name from blast filename
        query_name <- sub("^(.*)_vs_.*_blast_hits\\.tsv$",
                          "\\1",
                          blast_filename)
        
        # Create descriptive gene match filename
        gene_match_file_path <- file.path(
          gene_hit_dir_path,
          paste0(query_name, "_matches_in_", target_genome, ".csv")
        )
        
        # Save gene match results
        write_csv(gene_matches, file = gene_match_file_path)
        message(
          "    Found ",
          nrow(gene_matches),
          " gene matches, saved to: ",
          basename(gene_match_file_path)
        )
        
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
  
  
  
  ## Step 3: Extract and translate gene sequences to proteins 
  message("\n== STEP 3: Extracting and translating gene sequences ==")
  
  # Combine our results with the original genome mapping
  pipeline_results <- gene_identification_results |>
    left_join(genome_to_gff_mapping, by = "gff_file_path") |>
    # remove target_genome, so it won't appear twice in the next join
    dplyr::select(-target_genome)
  
  # Add query information
  pipeline_results <- left_join(pipeline_results, blast_results, by = "blast_hit_file_path")
  
  # Initialize data frame to track protein generation
  protein_generation_results <- data.frame(
    gene_match_file_path = character(0),
    transcript_id = character(0),
    protein_sequence_path = character(0),
    qseqid = character(0),
    # Query sequence ID
    seqid = character(0),
    # Added column to store seqid (subject sequence ID)
    gstart = numeric(0),  # Changed from character to numeric
    gend = numeric(0),   # Changed from character to numeric
    stringsAsFactors = FALSE
  )
  
  # Process each gene match file
  message("\nProcessing ",
          nrow(pipeline_results),
          " gene match files for protein generation")
  
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
        qseqid = qseqid[1],
        # Extract the first Query Sequence ID
        seqid = sseqid[1],
        # Extract the first Subject Sequence ID
        gstart = as.numeric(gstart[1]),  # Ensure numeric
        gend = as.numeric(gend[1]),     # Ensure numeric
        .groups = "drop"
      ) %>%
      dplyr::filter(!is.na(transcript_id))
    
    message("    Found ",
            nrow(transcript_info),
            " unique transcript(s) to process")
    
    # Process each transcript
    for (j in seq_len(nrow(transcript_info))) {
      transcript_id <- transcript_info$transcript_id[j]
      current_qseqid <- transcript_info$qseqid[j]
      current_seqid <- transcript_info$seqid[j] 
      current_gstart <- transcript_info$gstart[j]
      current_gend <- transcript_info$gend[j]
      
      message(
        "    Processing transcript: ",
        transcript_id,
        " (qseqid: ",
        current_qseqid,
        ", seqid: ",
        current_seqid,
        ")"
      )
      
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
        warning("Error processing transcript ",
                transcript_id,
                ": ",
                e$message)
      })
    }
  }
  
  # Save protein generation results
  protein_mapping_file <- file.path(output_dir, "protein_sequence_mapping.csv")
  write_csv(protein_generation_results, file = protein_mapping_file)
  message("Protein generation results saved to: ", protein_mapping_file)
  
  
  
  
  
  
  ## Step 4: Filter sequences by pair alignment and gen. aln's 
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
  
  if (generate_alignments && nrow(protein_generation_results) > 0) {
    message("\n== STEP 4: Filtering sequences by alignment score and generating alignments ==")
    
    # Group protein sequences by query
    # First, just group by query_name without worrying about query_path
    protein_by_query <- pipeline_results %>%
      left_join(protein_generation_results, by = "gene_match_file_path") %>%
      dplyr::filter(!is.na(protein_sequence_path)) %>%
      # Make sure query_name is defined
      mutate(
        query_name = ifelse(
          is.na(query_name),
          sub("^(.*)_matches_in_.*\\.csv$", "\\1", basename(gene_match_file_path)),
          query_name
        )
      ) %>%
      group_by(query_name) %>%
      summarise(
        protein_files = list(protein_sequence_path),
        num_sequences = n(),
        .groups = "drop"
      )
    
    # Assess and filter protein sequences for each query
    for (i in seq_len(nrow(protein_by_query))) {
      query_name <- protein_by_query$query_name[i]
      protein_files <- unlist(protein_by_query$protein_files[i])
      
      # Look up the query path from blast results
      query_path <- find_query_path_by_name(query_name, blast_results)
      
      # Skip if can't find query path
      if (is.na(query_path) || !file.exists(query_path)) {
        message("Warning: Cannot find query file for ", query_name, ". Skipping this query.")
        next
      }
      
      message("Evaluating and filtering sequences for query: ", query_name)
      message("  Analyzing ", length(protein_files), " sequences from different genomes")
      
      # Initialize scores_data for this query with correct columns
      scores_data <- data.frame(
        protein_sequence_path = character(0),
        pscore = numeric(0),
        pid = numeric(0),
        passed_filter = logical(0),
        stringsAsFactors = FALSE
      )
      
      ### Each protein seq: pairwise align and score 
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
          passed_filter <- (
            alignment_result$pscore >= min_score_threshold &&
              alignment_result$pid >= min_pid_threshold
          )
          
          # Add to scores_data with explicit columns
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
          
          if (passed_filter) {
            message(
              "    Passed filter: Score=",
              round(alignment_result$pscore, 1),
              ", PID=",
              round(alignment_result$pid, 1),
              "%"
            )
          } else {
            message(
              "    Failed filter: Score=",
              round(alignment_result$pscore, 1),
              ", PID=",
              round(alignment_result$pid, 1),
              "% (below thresholds)"
            )
          }
        } else {
          message("    Alignment failed, excluding from further analysis")
        }
      }
      
      # Add scores to the filtered_protein_results
      if (nrow(scores_data) > 0) {
        # Add a label for each point based on the filename
        scores_data$label <- basename(scores_data$protein_sequence_path)
        
        # Add query_name and append to filtered_protein_results
        filtered_protein_results <- rbind(
          filtered_protein_results,
          cbind(
            data.frame(query_name = query_name, stringsAsFactors = FALSE),
            scores_data
          )
        )
        
        ### Create score vs pid scatter plot 
        plot <- ggplot(data = scores_data, 
                       mapping = aes(x = pscore, y = pid, color = passed_filter)) +
          geom_point(size = 3) +
          labs(
            title = paste("Alignment Scores for", query_name),
            x = "BLOSUM Score",
            y = "Percent Identity"
          ) +
          scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "red"),
                             name = "Passed Filter") +
          geom_hline(yintercept = min_pid_threshold, linetype = "dashed") +
          geom_vline(xintercept = min_score_threshold, linetype = "dashed")
        
        # Create plot directory 
        plot_dir <- file.path(output_dir, "alignment_plots")
        if (!dir.exists(plot_dir)) {
          dir.create(plot_dir, recursive = TRUE)
        }
        
        # Save one plot per query
        ggsave(
          filename = paste0(query_name, "_alignment_scores.png"),
          plot = plot,
          path = plot_dir,
          width = 8,
          height = 6
        )
      }
      
      # Get filtered protein files
      filtered_proteins <- scores_data %>%
        dplyr::filter(passed_filter == TRUE) %>%
        pull(protein_sequence_path)
      
      # Only proceed with multiple sequence alignment if enough sequences passed filtering
      if (length(filtered_proteins) >= 2) {
        message(
          "  Proceeding with multiple sequence alignment using ",
          length(filtered_proteins),
          " filtered sequences"
        )
        
        ## Generate multiple sequence alignment with filtered sequences 
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
        message(
          "  Not enough sequences passed filtering criteria for ",
          query_name,
          ". Need at least 2 sequences for multiple alignment."
        )
      }
    }
    
    # Save filtering results
    filtering_results_file <- file.path(output_dir, "sequence_filtering_results.csv")
    write_csv(filtered_protein_results, file = filtering_results_file)
    message("Sequence filtering results saved to: ",
            filtering_results_file)
    
    # Save alignment results
    if (nrow(alignment_results) > 0) {
      alignment_mapping_file <- file.path(output_dir, "multiple_alignment_mapping.csv")
      write_csv(alignment_results, file = alignment_mapping_file)
      message("Multiple sequence alignment results saved to: ",
              alignment_mapping_file)
    } else {
      message("No multiple sequence alignments were generated after filtering.")
    }
  }
  
  
  
  
  ## Step 5: Combine pipeline results 
  message("\n== STEP 5: Combining pipeline results ==")
  
  # Join the protein generation results with pipeline results
  intermediate_results <- pipeline_results %>%
    left_join(protein_generation_results, by = "gene_match_file_path") %>%
    mutate(query_name = ifelse(
      is.na(query_name),
      sub(
        "^(.*)_matches_in_.*\\.csv$",
        "\\1",
        basename(gene_match_file_path)
      ),
      query_name
    ))
  
  # Add filtering information if available
  if (nrow(filtered_protein_results) > 0) {
    intermediate_results <- intermediate_results %>%
      left_join(
        filtered_protein_results %>%
          dplyr::select(
            query_name,
            protein_sequence_path,
            pscore,
            pid,
            passed_filter
          ),
        by = c("query_name", "protein_sequence_path")
      )
    
    message("Added alignment scores and filtering results to pipeline results")
  }
  
  # Add alignment information if available
  if (nrow(alignment_results) > 0) {
    intermediate_results <- intermediate_results %>%
      left_join(alignment_results %>% 
                  dplyr::select(query_name, alignment_file_path),
                by = "query_name")
    
    message("Added alignment file information to pipeline results")
  }
  
  # Final results
  final_results <- intermediate_results
  
  # Create a meaningful final mapping file
  final_mapping_file <- file.path(output_dir, "problaster_complete_pipeline_results.csv")
  write_csv(final_results, file = final_mapping_file)
  message("Complete pipeline results saved to: ", final_mapping_file)
  
  
  
  
  
  
  
  
  ## Step 6: Create phylogenetic trees if requested 
  tree_results <- NULL
  if (generate_trees && nrow(final_results) > 0) {
    message("\n== STEP 6: Creating phylogenetic trees ==")
    
    # Create trees using filtered sequences
    tree_results <- create_phylogenetic_trees(df = final_results, 
                                              min_sequences = min_tree_sequences)
  }
  
  
  
  
  
  
  
  
  
  # Generate a summary report
  pipeline_summary <- data.frame(
    stage = c(
      "BLAST searches",
      "Gene matches",
      "Protein sequences",
      "Sequences passing filters",
      "Multiple sequence alignments",
      "Phylogenetic trees"
    ),
    count = c(
      nrow(blast_results),
      nrow(gene_identification_results),
      nrow(protein_generation_results),
      if (generate_alignments)
        sum(filtered_protein_results$passed_filter, na.rm = TRUE)
      else
        NA,
      if (generate_alignments)
        nrow(alignment_results)
      else
        NA,
      if (generate_trees &&
          !is.null(tree_results))
        length(tree_results)
      else
        NA
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
    message("Sequences passing alignment filters: ",
            pipeline_summary$count[4])
    message("Multiple sequence alignments created: ",
            pipeline_summary$count[5])
  }
  if (generate_trees) {
    message("Phylogenetic trees created: ", pipeline_summary$count[6])
  }
  message("All results saved to: ", output_dir)
  message("Complete pipeline mapping: ", final_mapping_file)
  
  # Return the final results dataframe
  return(final_results)
}


# Utility Functions -------------------------------------------------------
#' Find query path by query name in the blast results
#'
#' @param query_name The name of the query to look up
#' @param blast_results The blast results data frame
#' @return The query path if found, NA otherwise
find_query_path_by_name <- function(query_name, blast_results) {
  query_info <- blast_results %>%
    dplyr::filter(query_name == !!query_name) %>%
    dplyr::select(query_path) %>%
    distinct()
  
  if (nrow(query_info) > 0) {
    return(query_info$query_path[1])
  } else {
    return(NA)
  }
}

clean_tip_label <- function(tip_labels) {
  # Initialize output vector
  new_labels <- character(length(tip_labels))
  
  for (i in 1:length(tip_labels)) {
    tip_label <- tip_labels[i]
    
    # Special case for query label
    if (tip_label == "[QUERY]") {
      new_labels[i] <- tip_label
      next
    }
    
    # Try to extract gene ID using various patterns
    gene_id <- NULL
    # Pattern for standard _g12345.t1 format
    gene_match <- regexpr("_g[0-9]+\\.t[0-9]+", tip_label)
    if (gene_match > 0) {
      gene_id <- substr(tip_label, gene_match + 1, 
                        gene_match + attr(gene_match, "match.length") - 1)
    } else {
      # Alternative pattern for other gene ID formats
      gene_match <- regexpr("g[0-9]+\\.t[0-9]+", tip_label)
      if (gene_match > 0) {
        gene_id <- substr(tip_label, gene_match, 
                          gene_match + attr(gene_match, "match.length") - 1)
      }
    }
    
    # Try to extract chromosome info
    chrom_info <- NULL
    chrom_patterns <- c("chr[0-9]+_hap[0-9]+", "chr[0-9]+", "scaffold[0-9]+")
    
    for (pattern in chrom_patterns) {
      chrom_match <- regexpr(pattern, tip_label)
      if (chrom_match > 0) {
        chrom_info <- substr(tip_label, chrom_match, 
                             chrom_match + attr(chrom_match, "match.length") - 1)
        break  # Exit loop after first match
      }
    }
    
    # Create new label based on what was found
    if (!is.null(gene_id) && !is.null(chrom_info)) {
      # Gene ID and chromosome found
      new_labels[i] <- paste0(gene_id, " (", chrom_info, ")")
    } else if (!is.null(gene_id)) {
      # Gene ID found
      new_labels[i] <- gene_id
    } else if (!is.null(chrom_info)) {
      # Only Chromosome info found
      new_labels[i] <- paste0("Unknown gene (", chrom_info, ")")
    } else {
      # Fallback: Extract everything before "| Source:" if present
      source_pos <- regexpr(" \\| Source:", tip_label)
      if (source_pos > 1) {
        new_labels[i] <- substr(tip_label, 1, source_pos - 1)
      } else {
        # If all else fails, keep original label
        new_labels[i] <- tip_label
      }
    }
  }
  
  return(new_labels)
}



# Example call ------------------------------------------------------------
# Example of how to run the complete pipeline
# don't end path with "/" or "\"
  results <- run_problaster_pipeline(
    genome_dir = "genomes/test_genomes",
    query_dir = "meiotic_genes_protein_fasta/test_meiotic_prot",
    gff_dir = "gff_files",
    index_path = "index.csv",
    # index of genom to gff mapping
    output_dir = "out_2025_03_07_tree_test",
    evalue = 1e-5,
    # tBLASTn evalue
    generate_alignments = TRUE,
    generate_trees = TRUE,
    min_score_threshold = 300,
    # Minimum acceptable BLOSUM62 score
    min_pid_threshold = 25,
    # Minimum percent identity (0-100)
    min_tree_sequences = 2
    # Minimum sequences to build a tree
  )
