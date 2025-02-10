suppressMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(Biostrings)
  library(GenomicFeatures)
  library(rtracklayer)
  library(Rsamtools)
})

create_blastdb_and_run_blast <- function(genome_dir, query_dir, output_dir, evalue = 1e-5) {
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  
  # Create BLAST databases
  genome_files <- list.files(genome_dir, pattern = "\\.fasta$", full.names = TRUE)
  message("Creating BLAST databases...")
  db_paths <- lapply(genome_files, function(genome_file) {
    db_name <- file.path(output_dir, paste0(tools::file_path_sans_ext(basename(genome_file)), "_db"))
    system2("makeblastdb",
            args = c("-in", genome_file,
                     "-dbtype", "nucl",
                     "-title", "potato_genome_db",
                     "-out", db_name,
                     "-parse_seqids"))
    return(db_name)
  })
  
  # Run BLAST searches
  message("Running BLAST searches...")
  query_files <- list.files(query_dir, pattern = "\\.fasta$", full.names = TRUE)
  
  for (query in query_files) {
    query_name <- tools::file_path_sans_ext(basename(query))
    for (db in unlist(db_paths)) {
      db_name <- basename(db)
      output_file <- file.path(output_dir, paste0(query_name, "_in_", db_name, ".tsv"))
      
      system2("tblastn",
              args = c("-query", query,
                       "-db", db,
                       "-out", output_file,
                       "-evalue", evalue,
                       "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"))
    }
  }
}

#find var names fix function
find_gene_from_blast_hit <- function(blast_hit#s, gff_file#s, min_seq_len) {
  
  #min_seq_len = 50000 # minimum sequence length
  
  # get the file names from the blast hit
  blast_hits <- list.files(path = "output_dir", pattern = "*.tsv")
  
  # get the file names from the annotated gff3 files
  gff_files <- list.files(path = "gff_files", pattern = "*.gff3")
  
  # Create output directory if it doesn't exist
  dir.create("output_dir", showWarnings = FALSE)
  
  #
  for (blast_hit in blast_hits) {
    for (gff_file in gff_files) {
      print(blast_hit)
      print(gff_file)
      
      
      # read database file
      blast_hits <- read.table(
        file = file.path("output_dir", blast_hit),
        sep = "\t",
        col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                      "qstart", "qend", "sstart", "send", "evalue", "bitscore")
      ) |>
        mutate(sseqid = as.character(sseqid)) # make sure that ssequid of both df will be of same type
      
      # filter the blast hits out which overlap or are more than 50 kb appart or have a
      # bitscore below 50
      filtered_blast_hit <- blast_hits |>
        mutate(low_start = pmin(sstart,send),
               high_end = pmax(sstart, send)) |>
        group_by(sseqid) |>
        arrange(low_start) |>
        mutate(
          prev_end = lag(high_end, default = head(high_end, 1)), # the first one is always a overlap?
          distance_to_prev = low_start - prev_end, # calculate the length of a sequence
          is_overlap = low_start <= prev_end # filter out sequence starts in a previous sequence
        ) |> 
        filter((distance_to_prev <= min_seq_len | is.na(distance_to_prev)) & !is_overlap) |>
        filter(bitscore >= 50) |> # only keep high bitscore
        ungroup() |>
        dplyr::select(-prev_end, -distance_to_prev, -is_overlap, low_start, high_end)
      
      # read gff file
      gff <- read.table(
        file = file.path("gff_files", gff_file),
        sep = "\t",
        col.names = c("sseqid", "source", "type", "gstart", "gend", 
                      "score", "strand", "phase", "gname")
      ) |>
        mutate(sseqid = as.character(sseqid)) |> # make sure that ssequid of both df will be of same type
        filter(type %in% c("mRNA", "transcript")) # type should be "mRNA" or "transcript"
      
      # merge the blast hits with the anotated genome based on the sequence id
      result <- inner_join(filtered_blast_hit, gff, by = "sseqid", relationship = "many-to-many") |>
        select(-contains("idk")) |>
        mutate(seq_in_gene = ifelse(sstart >= gstart & send <= gend, T, F)) |>
        filter(seq_in_gene == T) |>
        select(-seq_in_gene)
      
      # only save the csv files with results in them
      if (nrow(result) > 0) {
        write_csv(result,
                  file = file.path(
                    "output_dir",
                    paste0(
                      tools::file_path_sans_ext(blast_hit), 
                      "_with_",
                      tools::file_path_sans_ext(gff_file), 
                      "_blast_matches.csv")))
        
      } else {
        next
      }
    }
    
  }
}


# gene seq to prot --------------------------------------------------------


# Set working directory
setwd("~/Documents/potato_DNA")

# Files
genome_file <- "genomes/O_hap1_genome.fasta"
gff_file <- "gff_files/O_hap1_liftoff_low_high_confidence_valid_ORFs_only_no_Ns.gff3"

# Read gff file
gff <- import(gff_file)

# Filter for the gene and for CDS's
gff <- gff[gff$type == "CDS"]
cds <- gff[as.character(gff$Parent) == "anno1.g72247.t1"]

# makes sence?
# Sort cds based on strand
if (unique(strand(cds)) == "-") {
  cds <- sort(cds, decreasing = TRUE)
} else {
  cds <- sort(cds)
}

# Select genome file to extract the sequence further ahead
fa <- FaFile(genome_file)

# Extract sequences for each CDS
exon_seqs <- sapply(seq_along(cds), function(i) {
  chr <- as.character(seqnames(cds[i]))
  start <- start(cds[i])
  end <- end(cds[i])
  
  # Extract sequence from the genome
  seq <- scanFa(fa, param = GRanges(seqnames = chr, 
                                    ranges = IRanges(start = start, end = end)))
  
  # Reverse complement if on minus strand
  if (as.character(strand(cds[i])) == "-") {
    seq <- reverseComplement(seq)
  }
  
  return(as.character(seq))
})

# Join CDS sequences
full_transcript <- DNAString(paste(exon_seqs, collapse=""))

# Function to find and translate the longest Open Reading Frame
find_and_translate_longest_ORF <- function(sequence) {
  # Find all ATG positions
  atg_positions <- start(matchPattern("ATG", sequence))
  longest_protein <- ""
  
  for (start_pos in atg_positions) {
    # Get subsequence starting from ATG
    current_seq <- subseq(sequence, start_pos)
    
    # Ensure sequence length is a multiple of 3 before translating
    remainder <- length(current_seq) %% 3
    if (remainder != 0) {
      print("not multiple of 3") # testing
    }
    
    # Only go ahead if sequence length is longer than 3
    if (length(current_seq) >= 3) {
      
      # Translate the sequence to the protein sequence
      protein <- translate(current_seq)
      
      # Find the position of the first stop codon
      stop_pos <- min(c(
        start(matchPattern("*", protein)),
        length(protein) + 1
      ))
      
      # Get the protein sequence up to the stop codon
      current_protein <- substr(as.character(protein), 1, stop_pos - 1)
      
      # Update if this is the longest protein found
      if (nchar(current_protein) > nchar(longest_protein)) {
        longest_protein <- current_protein
      }
    }
  }
  
  return(longest_protein)
}

# Get the protein sequence
protein <- find_and_translate_longest_ORF(full_transcript)

# Get the gene name from the cds object
gene_name <- unique(as.character(mcols(cds)$Parent))

# Write protein output with header
protein_header <- paste0(">", gene_name, " - Longest Open Reading Frame")
writeLines(c(protein_header, protein), "protein.fasta")

# Write transcript output with header
transcript_header <- paste0(">", gene_name, " - Full Transcript Sequence")
writeLines(c(transcript_header, as.character(full_transcript)), "transcript.fasta")









# main --------------------------------------------------------------------

main <- function(variables) {
  # empty
}
