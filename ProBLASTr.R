# Script Name: ProBLASTr.R
# Description: An R script for protein BLAST analysis and homology verification
#
# Copyright (C) [2025]
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
# Version: 1.0
# Last Updated: [2025.02.26]
#
# Project Repository: [URL to your repository]
#
# Dependencies:
# some missing
# blastn: 2.16.0+
#   Package: blast 2.16.0, build Aug  7 2024 01:45:09
#
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
# ProBLASTr functions -----------------------------------------------------


# Load required libraries
library(tidyverse)
library(Biostrings)
library(GenomicRanges)
library(rtracklayer)
library(stringr)
library(msa)

# -- BLAST Database Functions ----------------------------------------------------------

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
  
  # Create directories
  blast_database_dir <- file.path(output_dir, "blast_database")
  blast_hits_dir <- file.path(output_dir, "blast_hits")
  for (dir in c(blast_database_dir, blast_hits_dir)) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
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
        db_name
      )
    )
    # keep track of the database paths
    db_paths[i] <- db_name
  }
  
  # List and validate query files
  query_files <- list.files(query_dir, pattern = qur_file_pattern, full.names = TRUE)
  if (length(query_files) == 0) {
    stop("No query files found!")
  }
  
  message("Running BLAST searches...")
  
  # Initialize index for successful BLAST runs
  index_blast <- data.frame(
    query_path = character(0),
    blast_hit_file_path = character(0),
    genome_file = character(0)
  )
  
  # Run BLAST for each query against each database
  for (query in query_files) {
    query_name <- tools::file_path_sans_ext(basename(query))
    for (i in seq_along(db_paths)) {
      db <- db_paths[i]
      genome_file <- genome_files[i]
      message(sprintf("Query: '%s'\nDb: '%s'\n", query, db))
      
      db_name <- basename(db)
      output_file <- file.path(blast_hits_dir,
                               paste0(query_name, "_in_", db_name, ".tsv"))
      
      # Run tBLASTn with properly escaped quotes
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
        file.remove(output_file)
      }
      
      # Record successful BLAST runs
      if (file.exists(output_file) && file.size(output_file) > 0) {
        index_blast <- rbind(
          index_blast,
          data.frame(
            query_path = query,
            blast_hit_file_path = output_file,
            genome_file = genome_file
          )
        )
      }
    }
  }
  
  return(index_blast)
}

# -- Sequence Translation Functions ---------------------------------------------------

#' Extract and translate DNA sequence for every BLAST hit
#'
#' @param blast_hit_file Path to the BLAST hit TSV file
#' @param genome_file Path to the genome FASTA file
#' @param output_dir Directory for output files
#' @param extend_bp Number of base pairs to extend the hit region (to capture complete coding frames)
#' @return DataFrame with paths to translated protein sequences
translate_blast_hits <- function(blast_hit_file, genome_file, output_dir, extend_bp = 300) {
  # Create output directory
  protein_dir <- file.path(output_dir, "direct_translations")
  if (!dir.exists(protein_dir)) {
    dir.create(protein_dir, recursive = TRUE)
  }
  
  # Read BLAST hits
  tryCatch({
    hits <- read.table(
      file = blast_hit_file,
      sep = "\t",
      col.names = c(
        "qseqid", "sseqid", "pident", "length", "mismatch", 
        "gapopen", "qstart", "qend", "sstart", "send", 
        "evalue", "bitscore"
      )
    )
  }, error = function(e) {
    warning(paste("Error reading BLAST hit file:", blast_hit_file, "-", e$message))
    return(data.frame())
  })
  
  if (nrow(hits) == 0) {
    return(data.frame())
  }
  
  # Load genome
  tryCatch({
    genome <- readDNAStringSet(genome_file)
    names(genome) <- sub(" .*", "", names(genome))
  }, error = function(e) {
    warning(paste("Error reading genome file:", genome_file, "-", e$message))
    return(data.frame())
  })
  
  results <- data.frame(
    hit_id = character(0),
    protein_file = character(0),
    frame = integer(0),
    aa_length = integer(0),
    query_id = character(0),
    hit_chromosome = character(0),
    hit_start = integer(0),
    hit_end = integer(0),
    e_value = numeric(0),
    percent_identity = numeric(0)
  )
  
  # Process each hit
  for (i in 1:nrow(hits)) {
    hit <- hits[i, ]
    
    # Determine coordinates, accounting for reverse strand hits
    is_reverse <- hit$sstart > hit$send
    start_pos <- min(hit$sstart, hit$send) - extend_bp
    end_pos <- max(hit$sstart, hit$send) + extend_bp
    
    # Ensure coordinates are valid
    start_pos <- max(1, start_pos)
    
    # Only proceed if the hit chromosome is in the genome
    if (hit$sseqid %in% names(genome)) {
      # Extract DNA sequence
      seq_length <- length(genome[[hit$sseqid]])
      end_pos <- min(end_pos, seq_length)
      dna_seq <- subseq(genome[[hit$sseqid]], start = start_pos, end = end_pos)
      
      # Translate in all 6 frames
      frames <- list()
      # Forward frames
      frames[[1]] <- translate(dna_seq, if.fuzzy.codon = "solve")  # Frame 0
      frames[[2]] <- translate(subseq(dna_seq, start=2), if.fuzzy.codon = "solve")  # Frame 1
      frames[[3]] <- translate(subseq(dna_seq, start=3), if.fuzzy.codon = "solve")  # Frame 2
      # Reverse frames
      rev_seq <- reverseComplement(dna_seq)
      frames[[4]] <- translate(rev_seq, if.fuzzy.codon = "solve")  # Frame 0, reverse
      frames[[5]] <- translate(subseq(rev_seq, start=2), if.fuzzy.codon = "solve")  # Frame 1, reverse
      frames[[6]] <- translate(subseq(rev_seq, start=3), if.fuzzy.codon = "solve")  # Frame 2, reverse
      
      # Find the frame that gives the longest ORF
      orfs <- lapply(frames, function(f) {
        # Split by stop codons
        orfs <- strsplit(as.character(f), "\\*")[[1]]
        # Return the longest ORF
        return(orfs[which.max(nchar(orfs))])
      })
      
      best_frame <- which.max(sapply(orfs, nchar))
      best_orf <- orfs[[best_frame]]
      
      # Generate unique hit ID
      hit_id <- paste0(
        hit$qseqid, "_", 
        hit$sseqid, "_", 
        start_pos, "-", end_pos, 
        "_frame", best_frame
      )
      
      # Clean ID for file naming
      clean_id <- gsub("[^a-zA-Z0-9._-]", "_", hit_id)
      
      # Save protein sequence to file
      protein_file <- file.path(protein_dir, paste0(clean_id, ".fasta"))
      writeLines(
        c(
          paste0(">", hit_id, " | ", 
                 "Query:", hit$qseqid, 
                 " Genome:", basename(genome_file),
                 " Identity:", hit$pident, "%",
                 " E-value:", hit$evalue),
          best_orf
        ),
        protein_file
      )
      
      # Record result
      results <- rbind(
        results,
        data.frame(
          hit_id = hit_id,
          protein_file = protein_file,
          frame = best_frame,
          aa_length = nchar(best_orf),
          query_id = hit$qseqid,
          hit_chromosome = hit$sseqid,
          hit_start = start_pos,
          hit_end = end_pos,
          e_value = hit$evalue,
          percent_identity = hit$pident
        )
      )
    }
  }
  
  return(results)
}

# -- Sequence Analysis Functions -----------------------------------------------------
#' Align translated sequences with original query sequences
#'
#' @param protein_file Path to the translated protein sequence file
#' @param query_file Path to the original query protein file
#' @param output_dir Directory for output files
#' @return Path to the alignment file
align_sequences <- function(protein_file, query_file, output_dir) {
  
  align_dir <- file.path(output_dir, "alignments")
  
  if (!(dir.exists(align_dir))) {
    dir.create(align_dir, recursive = TRUE)
  }
  
  # Generate output file name
  alignment_name <- paste0(
    tools::file_path_sans_ext(basename(protein_file)), 
    "_vs_", 
    tools::file_path_sans_ext(basename(query_file)),
    ".aln"
  )
  alignment_file <- file.path(align_dir, alignment_name)
  
  # Read sequences
  tryCatch({
    protein_seq <- readAAStringSet(protein_file)
    query_seq <- readAAStringSet(query_file)
    
    # If there are multiple sequences in the query file, try to find the right one
    if (length(query_seq) > 1) {
      # Extract query ID from protein file header
      protein_header <- names(protein_seq)[1]
      query_id_match <- str_extract(protein_header, "Query:[^\\s|]+")
      if (!is.na(query_id_match)) {
        query_id <- gsub("Query:", "", query_id_match)
        # Find the matching query sequence
        matching_indices <- grep(query_id, names(query_seq))
        if (length(matching_indices) > 0) {
          query_seq <- query_seq[matching_indices[1]]
        } else {
          # If no match found, just use the first query sequence
          query_seq <- query_seq[1]
        }
      } else {
        # If no query ID found in header, use the first query sequence
        query_seq <- query_seq[1]
      }
    }
    
    # Add a description to indicate which is query and which is translated hit
    names(protein_seq) <- paste0(names(protein_seq), " [TRANSLATED]")
    names(query_seq) <- paste0(names(query_seq), " [QUERY]")
    
    # Combine sequences
    combined_seqs <- c(query_seq, protein_seq)
    
    # Perform alignment using MUSCLE algorithm from msa package
    alignment <- msa(combined_seqs, method = "Muscle")
    
    # Write alignment as FASTA file instead of CLUSTAL format
    # This is the key change - using "fasta" format instead of "clustal"
    writeXStringSet(
      as(alignment, "AAStringSet"),
      filepath = alignment_file,
      format = "fasta"  # Changed from "clustal" to "fasta"
    )
    
    # Check if alignment file was created successfully
    if (!file.exists(alignment_file) || file.size(alignment_file) == 0) {
      warning("Alignment failed for ", basename(protein_file))
      return(NA)
    }
    
    return(alignment_file)
    
  }, error = function(e) {
    warning("Error performing alignment: ", e$message)
    return(NA)
  })
}

#' Run HMMER search to identify domains in protein sequences
#'
#' @param protein_file Path to the protein sequence file
#' @param hmm_db Path to the HMMER database (Pfam-A.hmm)
#' @param output_dir Directory for output files
#' @param evalue E-value threshold for domain hits
#' @return Path to the HMMER output file
identify_domains <- function(protein_file, hmm_db, output_dir, evalue = 1e-3) {
  # Create domains directory
  domains_dir <- file.path(output_dir, "domains")
  if (!dir.exists(domains_dir)) {
    dir.create(domains_dir, recursive = TRUE)
  }
  
  # Generate output file name
  output_name <- paste0(tools::file_path_sans_ext(basename(protein_file)), "_domains.txt")
  output_file <- file.path(domains_dir, output_name)
  
  # Run hmmscan
  tryCatch({
    system2(
      "hmmscan",
      args = c(
        "--domtblout", output_file,
        "-E", evalue,
        "--cpu", 2,
        hmm_db,
        protein_file
      )
    )
    
    # Check if output file was created successfully
    if (!file.exists(output_file) || file.size(output_file) == 0) {
      warning("No domains found in ", basename(protein_file))
      # Create empty file with header to indicate no domains found
      writeLines("# No domains found", output_file)
    }
    
  }, error = function(e) {
    warning("Error running HMMER: ", e$message)
    return(NA)
  })
  
  return(output_file)
}

#' Parse HMMER domain output into a data frame
#'
#' @param domain_file Path to the HMMER domain output file
#' @return Data frame containing domain information
parse_domain_file <- function(domain_file) {
  # Check if file exists and is not empty
  if (!file.exists(domain_file) || file.size(domain_file) == 0) {
    return(data.frame())
  }
  
  # Read domain file
  lines <- readLines(domain_file)
  
  # Filter out comment lines
  lines <- lines[!grepl("^#", lines)]
  
  # If no content lines, return empty data frame
  if (length(lines) == 0) {
    return(data.frame())
  }
  
  # Parse each line
  domains <- lapply(lines, function(line) {
    # Split by whitespace
    fields <- strsplit(line, "\\s+")[[1]]
    fields <- fields[fields != ""]
    
    # HMMER domtblout format has specific columns
    if (length(fields) >= 22) {
      return(data.frame(
        target_name = fields[1],
        target_accession = fields[2],
        query_name = fields[4],
        query_accession = fields[5],
        e_value = as.numeric(fields[6]),
        score = as.numeric(fields[7]),
        domain_start = as.numeric(fields[17]),
        domain_end = as.numeric(fields[18]),
        stringsAsFactors = FALSE
      ))
    } else {
      return(NULL)
    }
  })
  
  # Combine into a single data frame
  domains <- do.call(rbind, domains[!sapply(domains, is.null)])
  
  return(domains)
}

# -- Main Function -------------------------------------------------------------------

#' Main function to run the complete analysis pipeline
#'
#' @param genome_dir Directory containing genome FASTA files
#' @param query_dir Directory containing query protein FASTA files
#' @param output_dir Directory for output files
#' @param evalue E-value threshold for BLAST searches (default: 1e-5)
#' @param extend_bp Number of base pairs to extend around each hit (default: 300)
#' @param hmm_db Path to HMMER database for domain identification (default: NULL)
#' @param run_alignment Whether to run sequence alignments (default: TRUE)
#' @param run_domain_search Whether to run domain searches (default: TRUE if hmm_db is provided)
#' @return A list containing results and summary information
main <- function(genome_dir,
                               query_dir,
                               output_dir,
                               evalue = 1e-5,
                               extend_bp = 300,
                               hmm_db = NULL,
                               run_alignment = TRUE,
                               run_domain_search = !is.null(hmm_db)) {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Step 1: Create BLAST databases and run searches
  message("Step 1: Running BLAST searches...")
  blast_index <- create_blastdb_and_run_tblastn(
    genome_dir = genome_dir,
    query_dir = query_dir,
    output_dir = output_dir,
    evalue = evalue
  )
  
  # Save the BLAST index
  write_csv(blast_index, file.path(output_dir, "blast_index.csv"))
  
  # Step 2: Translate all BLAST hits
  message("Step 2: Translating BLAST hits to protein sequences...")
  all_translations <- data.frame()
  
  for (i in 1:nrow(blast_index)) {
    blast_hit_file <- blast_index$blast_hit_file_path[i]
    genome_file <- blast_index$genome_file[i]
    query_file <- blast_index$query_path[i]
    
    message(paste0("Processing hit file ", i, " of ", nrow(blast_index), ": ", 
                   basename(blast_hit_file)))
    
    translations <- translate_blast_hits(
      blast_hit_file = blast_hit_file,
      genome_file = genome_file,
      output_dir = output_dir,
      extend_bp = extend_bp
    )
    
    if (nrow(translations) > 0) {
      # Add source file information
      translations$blast_hit_file <- blast_hit_file
      translations$genome_file <- genome_file
      translations$query_file <- query_file
      
      all_translations <- rbind(all_translations, translations)
    }
  }
  
  # Save all translation results
  write_csv(all_translations, file.path(output_dir, "all_translations_index.csv"))
  
  # Step 3: Sequence Analysis - Alignment
  if (run_alignment && nrow(all_translations) > 0) {
    message("Step 3a: Aligning translated proteins with query sequences...")
    
    # Initialize columns for alignment results
    all_translations$alignment_file <- NA
    
    # Process each translation
    for (i in 1:nrow(all_translations)) {
      protein_file <- all_translations$protein_file[i]
      query_file <- all_translations$query_file[i]
      
      if (i %% 10 == 0 || i == 1 || i == nrow(all_translations)) {
        message(paste0("Aligning sequence ", i, " of ", nrow(all_translations)))
      }
      
      # Run alignment
      alignment_file <- align_sequences(
        protein_file = protein_file,
        query_file = query_file,
        output_dir = output_dir
      )
      
      # Update data frame
      all_translations$alignment_file[i] <- alignment_file
    }
    
    # Save updated translations index
    write_csv(all_translations, file.path(output_dir, "all_translations_index.csv"))
  }
  
  # Step 4: Sequence Analysis - Domain Identification
  if (run_domain_search && nrow(all_translations) > 0) {
    if (is.null(hmm_db) || !file.exists(hmm_db)) {
      warning("HMMER database not found. Skipping domain identification.")
    } else {
      message("Step 3b: Identifying protein domains...")
      
      # Initialize columns for domain results
      all_translations$domain_file <- NA
      all_domains <- data.frame()
      
      # Process each translation
      for (i in 1:nrow(all_translations)) {
        protein_file <- all_translations$protein_file[i]
        
        if (i %% 10 == 0 || i == 1 || i == nrow(all_translations)) {
          message(paste0("Searching domains in sequence ", i, " of ", nrow(all_translations)))
        }
        
        # Run domain search
        domain_file <- identify_domains(
          protein_file = protein_file,
          hmm_db = hmm_db,
          output_dir = output_dir,
          evalue = 1e-3
        )
        
        # Update data frame
        all_translations$domain_file[i] <- domain_file
        
        # Parse domain results
        domains <- parse_domain_file(domain_file)
        if (nrow(domains) > 0) {
          domains$protein_file <- protein_file
          domains$hit_id <- all_translations$hit_id[i]
          domains$query_id <- all_translations$query_id[i]
          all_domains <- rbind(all_domains, domains)
        }
      }
      
      # Save domain results
      if (nrow(all_domains) > 0) {
        write_csv(all_domains, file.path(output_dir, "domain_results.csv"))
        
        # Generate domain summary
        domain_summary <- all_domains %>%
          group_by(query_id, target_name) %>%
          summarize(
            count = n(),
            avg_score = mean(score),
            min_e_value = min(e_value)
          )
        
        write_csv(domain_summary, file.path(output_dir, "domain_summary.csv"))
      }
      
      # Save updated translations index
      write_csv(all_translations, file.path(output_dir, "all_translations_index.csv"))
    }
  }
  
  # Generate summary statistics
  if (nrow(all_translations) > 0) {
    summary_stats <- all_translations %>%
      group_by(query_id) %>%
      summarize(
        total_hits = n(),
        avg_length = mean(aa_length),
        min_length = min(aa_length),
        max_length = max(aa_length),
        avg_identity = mean(percent_identity)
      )
    
    # Save summary statistics
    write_csv(summary_stats, file.path(output_dir, "translation_summary_stats.csv"))
    
    # Print summary
    message("\nSummary of results:")
    message(paste0("Total translated protein sequences: ", nrow(all_translations)))
    message(paste0("Number of queries with hits: ", length(unique(all_translations$query_id))))
    message(paste0("Number of genomes with hits: ", length(unique(all_translations$genome_file))))
    
    if (exists("all_domains") && nrow(all_domains) > 0) {
      message(paste0("Number of domains identified: ", nrow(all_domains)))
      message(paste0("Number of unique domain types: ", length(unique(all_domains$target_name))))
    }
  } else {
    message("\nNo BLAST hits were found or translated.")
  }
  
  # Return all results
  results <- list(
    blast_index = blast_index,
    translations = all_translations
  )
  
  if (exists("all_domains") && nrow(all_domains) > 0) {
    results$domains <- all_domains
  }
  
  return(results)
}

# -- Example Usage ------------------------------------------------------------------

# Example command to run the analysis:
results <- main(
  genome_dir = "genomes/test_genomes",
  query_dir = "meiotic_genes_protein_fasta/test_meiotic_prot",
  output_dir = "out_2025_02_26v3",
  evalue = 1e-5,
  extend_bp = 300,
  hmm_db = "pfam_db/Pfam-A.hmm"
)
