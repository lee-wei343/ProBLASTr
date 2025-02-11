library(tidyverse)
library(GenomicRanges)
library(Biostrings)
library(GenomicFeatures)
library(rtracklayer)
library(Rsamtools)
library(stringi)

create_blastdb_and_run_tblastn <- function(genome_dir, query_dir, output_dir, evalue) {
  
  # Create the folder "blast_database" if not exist already
  blast_database_dir <- file.path(output_dir, "blast_database")
  if (!dir.exists(blast_database_dir)) {
    dir.create(blast_database_dir, recursive = TRUE)
  }
  
  # Create BLAST databases
  genome_files <- list.files(genome_dir, pattern = "\\.fasta$", full.names = TRUE)
  message("Creating BLAST databases...")
  db_paths <- lapply(genome_files, function(genome_file) {
    db_name <- file.path(blast_database_dir, paste0(tools::file_path_sans_ext(basename(genome_file)), "_db"))
    system2("makeblastdb",
            args = c("-in", genome_file,
                     "-dbtype", "nucl",
                     "-out", db_name,
                     "-parse_seqids"))
    return(db_name)
  })
  
  # Run BLAST searches
  message("Running BLAST searches...")
  query_files <- list.files(query_dir, pattern = "\\.fasta$", full.names = TRUE)
  
  # Create the folder "blast_hits" if not exist already
  blast_hits_dir <- file.path(output_dir, "blast_hits")
  if (!dir.exists(blast_hits_dir)) {
    dir.create(blast_hits_dir, recursive = TRUE)
  }
  
  for (query in query_files) {
    query_name <- tools::file_path_sans_ext(basename(query))
    for (db in unlist(db_paths)) {
      db_name <- basename(db)
      output_file <- file.path(output_dir, "blast_hits", 
                               paste0(query_name, "_in_", db_name, ".tsv"))
      
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
find_gene_from_blast_hit <- function(blast_hit, gff_file, min_seq_len) {
  
  # filter the blast hits out which overlap or are more than 50 kb appart or have a
  # bitscore below 50
  filtered_blast_hit <- blast_hit |>
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
  
  
  # merge the blast hits with the anotated genome based on the sequence id
  result <- inner_join(filtered_blast_hit, gff_file, by = "sseqid", 
                       relationship = "many-to-many") |>
    dplyr::select(-contains("idk")) |>
    mutate(seq_in_gene = ifelse(sstart >= gstart & send <= gend, T, F)) |>
    filter(seq_in_gene == T) |>
    dplyr::select(-seq_in_gene)
  
  # only return not empty df's
  if (nrow(result) > 0) {
    return(result)
  } else {
    return(NULL)
  }
}


# gene seq to prot --------------------------------------------------------
gene_seq_to_prot <- function(genome_file, gff_file, transcript_id) {
  # parameters are corresponding genome file, gff file and the transcript_id found in
  # the blast_hit output
  
# genome_file <- "genomes/O_hap1_genome.fasta"
# gff_file <- "gff_files/O_hap1_liftoff_low_high_confidence_valid_ORFs_only_no_Ns.gff3"

# Read gff file
gff <- import(gff_file)

# Filter for the gene and for CDS's
gff <- gff[gff$type == "CDS"]
cds <- gff[as.character(gff$Parent) == transcript_id] # maybe parent or transcript_id column

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
cds_seqs <- sapply(seq_along(cds), function(i) {
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
}


# main --------------------------------------------------------------------
main <- function(genome_dir, 
                 query_dir, 
                 gff_dir, 
                 output_dir, 
                 evalue = 1e-5, 
                 min_seq_len = 50000) {
  
  # create outputdir if not exist already
  if (!(dir.exists(output_dir))) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  # # create a blast datavase and run tblastn with the query fastas against it
  # create_blastdb_and_run_tblastn(genome_dir, query_dir, output_dir, evalue)
  
  # get the file names from the blast hit
  blast_hits <- list.files(path = file.path(output_dir, "blast_hits"),
                           pattern = "\\.tsv$")
  
  # get the file names from the annotated gff3 files
  gff_files <- list.files(path = gff_file, pattern = "\\.gff3$")
  
  # create the gene_hits directory if not exist already
  gene_hit_path <- file.path(output_dir, "gene_hits")
  if (!dir.exists(gene_hit_path)) {
    dir.create(gene_hit_path, recursive = TRUE)
  }
  
  # create empty dataframe for the file paths, to save the iterated path to an index 
  # dataframe
  index <- data.frame(
    blast_hit_file_path <- character(0),
    gff_file_path <- character(0),
    stringsAsFactors = FALSE
  )
  
  # loop over the blast results and read them in
  for (blast_hit in blast_hits) {
    
    # get pathname
    blast_hit_file_path <- file.path(output_dir, "blast_hits", blast_hit)
    
    # read database file
    blast_hit_table <- read.table(
      file = blast_hit_file_path,
      sep = "\t",
      col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                    "qstart", "qend", "sstart", "send", "evalue", "bitscore")
    ) |>
      mutate(sseqid = as.character(sseqid)) # make sure that ssequid of both df will be of same type
    
    # get the unique sseqid 
    blast_hit_unique_sseqid <- paste(unique(blast_hit_table$sseqid), collapse = "|")
    
    
    
    # find gff files that match the unique sseqid
    matching_gff_file <- system2("egrep",
                                 args = c("-l", shQuote(blast_hit_unique_sseqid), 
                                          file.path(gff_dir, "*")),
                                 stdout = TRUE) 
    
    # loop over the matching gff files 
    for (gff_file_hit in basename(matching_gff_file)) {
      
      # get pathname
      gff_file_path <- file.path(gff_dir, gff_file_hit)
      
      # read gff file
      gff_file_table <- read.table(
        file = gff_file_path,
        sep = "\t",
        col.names = c("sseqid", "source", "type", "gstart", "gend", 
                      "score", "strand", "phase", "gname")
      ) |>
        mutate(sseqid = as.character(sseqid)) |> # make sure that ssequid of both df will be of same type
        filter(type %in% c("mRNA", "transcript"))
      
      
      result <- find_gene_from_blast_hit(blast_hit = blast_hit_table,
                               gff_file = gff_file_table,
                               min_seq_len = min_seq_len)
      if (!(is.null(result)) && is.data.frame(result) && nrow(result) > 0) {
       write_csv(result,
                  file = file.path(
                    gene_hit_path,
                    paste0(
                      tools::file_path_sans_ext(blast_hit), 
                      "_with_",
                      tools::file_path_sans_ext(gff_file_hit), 
                      "_gene_hit.csv"))) 
        
       # add pathnames to index data frame
       index <- rbind(index, data.frame(
         blast_hit_file_path = blast_hit_file_path,
         gff_file_path = gff_file_path
       )) 
      }
      next
    }
  }
  
  # read in an index file which should contain gff_file_path with matching genome_file_path
  index_from_file <- read_csv(file = "index.csv")
  
  # join the index's
  index <- inner_join(index, index_from_file, by = "gff_file_path") |>
    na.omit()
  
  for (genome_file_path in index$genome_file_path) {
    for (gff_file_path in index$gff_file_path) {
      for (blast_hit_file_path in index$blast_hit_file_path) {
        
        blast_hit_file_path_df <- read_csv(file = blast_hit_file_path)
        
        transcript_id_raw <- blast_hit_file_path_df$gname # ID=TARGET;
        transcript_id <- stri_extract(transcript_id_raw,
                                               regex = "ID=.*?;")
        
        # call based on index
        gene_seq_to_prot(genome_file, gff_file, transcript_id = transcript_id)
      }
    }
  }
}

main(genome_dir = "genomes",
     query_dir = "meiotic_genes_protein_fasta",
     gff_dir = "gff_files",
     output_dir = "tco_out")


create_blastdb_and_run_tblastn(genome_dir = "genomes",
     query_dir = "meiotic_genes_protein_fasta",
     output_dir = "tco_out",
     evalue = 1e-5)


