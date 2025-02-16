library(tidyverse)
library(GenomicRanges)
library(Biostrings)
library(GenomicFeatures)
library(rtracklayer)
library(Rsamtools)
library(stringi)


##
#  options(error = NULL)
# 
options(error = function(){
  beepr::beep(9)
  Sys.sleep(2)
})

##

# create BLAST database and run tBLASTn -----------------------------------
create_blastdb_and_run_tblastn <- function(genome_dir, query_dir, output_dir,
                                           evalue = 1e-5, 
                                           dbtype = "nucl", 
                                           gen_file_pattern = "\\.fasta$",
                                           qur_file_pattern = "\\.fasta$") {
  # Check if genome_dir, query_dir exist
  if (!(dir.exists(genome_dir))) {
    stop(paste0("Genome directory '", genome_dir, "' not found!"))
  }
  if (!(dir.exists(query_dir))) {
    stop(paste0("Query directory '", query_dir, "' not found!"))
  }
  
  # Create the folder "blast_database" if not exist already
  blast_database_dir <- file.path(output_dir, "blast_database")
  if (!dir.exists(blast_database_dir)) {
    dir.create(blast_database_dir, recursive = TRUE)
  }
  
  # Create BLAST databases
  # List genome files
  genome_files <- list.files(genome_dir, pattern = gen_file_pattern, full.names = TRUE)
  
  # Check if there are any files
  if (length(genome_files) == 0) {
    stop(paste0("No genome file found!"))
  }
  
  message("Creating BLAST databases...")
  
  # initialise empty vector for the database path's
  db_paths <- character(length(genome_files)) 
  
  for (i in seq_along(genome_files)) {
    genome_file <- genome_files[i]
    db_name <- file.path(blast_database_dir,
                         paste0(tools::file_path_sans_ext(basename(genome_file)), "_db"))
    
    # Run Bash makeblastdb
    system2("makeblastdb",
            args = c("-in", genome_file,
                     "-dbtype", dbtype,
                     "-out", db_name,
                     "-parse_seqids"))
    
    # add databasename to database path vector
    db_paths[i] <- db_name
  }
  
  
  # Run BLAST searches
  # List query files
  query_files <- list.files(query_dir, pattern = qur_file_pattern, full.names = TRUE)
  
  # Check if there are any files
  if (length(query_files) == 0) {
    stop(paste0("No query files found!"))
  }
  
  message("Running BLAST searches...")
  
  # Create the folder "blast_hits" if not exist already
  blast_hits_dir <- file.path(output_dir, "blast_hits")
  if (!dir.exists(blast_hits_dir)) {
    dir.create(blast_hits_dir, recursive = TRUE)
  }
  
  for (query in query_files) {
    query_name <- tools::file_path_sans_ext(basename(query))
    for (db in unlist(db_paths)) {
      
      # keep track of tblastn run
      cat(paste0("Query: '", query, "'\n", "Db: '", db, "'\n"))
      
      db_name <- basename(db)
      output_file <- file.path(output_dir, "blast_hits", 
                               paste0(query_name, "_in_", db_name, ".tsv"))
      
      # Run Bash tblastn
      system2("tblastn",
              args = c("-query", query,
                       "-db", db,
                       "-out", output_file,
                       "-evalue", evalue,
                       "-outfmt", "\"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\""))
      
      # ensure no empty files are left behind if the evalue threshhold filters out all hits
      if (file.exists(output_file)) {
        if (file.size(output_file) == 0) {
          file.remove(output_file)
        }
      } 
    }
  }
}


# find gene from blast hit ------------------------------------------------
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
gene_seq_to_prot <- function(genome_file, gff_file, transcript_id, output_dir,
                             gene_hit_path) {
  
  ##
  # Debug prints at start of function
  message("Inside gene_seq_to_prot(), received parameters:")
  message("genome_file = ", genome_file)
  message("gff_file = ", gff_file)
  message("transcript_id = ", transcript_id)
  message("output_dir = ", output_dir)
  ##
  
  # parameters are corresponding genome file, gff file and the transcript_id found in
  # the blast_hit output
  
  # genome_file <- "genomes/O_hap1_genome.fasta"
  # gff_file <- "gff_files/O_hap1_liftoff_low_high_confidence_valid_ORFs_only_no_Ns.gff3"
  
  # Read gff file
  gff <- import(gff_file)
  
  # Filter for the gene and for CDS's
  gff <- gff[gff$type == "CDS"]
  cds <- gff[as.character(gff$Parent) == transcript_id] # maybe parent or transcript_id column
  
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
  full_transcript <- DNAString(paste(cds_seqs, collapse=""))
  
  # Function to find and translate the longest Open Reading Frame
  find_and_translate_longest_ORF <- function(sequence) {
    # Find all ATG positions
    atg_positions <- start(matchPattern("ATG", sequence))
    longest_protein <- ""
    
    for (start_pos in atg_positions) {
      # Get subsequence starting from ATG
      current_seq <- subseq(sequence, start_pos)
      
      # # Ensure sequence length is a multiple of 3 before translating
      # remainder <- length(current_seq) %% 3
      # if (remainder != 0) {
      #   print("not multiple of 3") # testing
      # }
      # => sequence gets cut at the end later
      
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
  
  # Create folder for protein sequences
  protein_output_dir_path <- file.path(output_dir, "protein")
  if (!(dir.exists(protein_output_dir_path))) {
    dir.create(protein_output_dir_path)
  }
  
  # Create folder for transcript sequences
  transcript_output_dir_path <- file.path(output_dir, "transcript")
  if (!(dir.exists(transcript_output_dir_path))) {
    dir.create(transcript_output_dir_path)
  }
  
  # Write protein output with header
  protein_header <- paste0(">", gene_name, 
                           " Gene hit: ", basename(gene_hit_path),
                           " Genome: ", basename(genome_file),
                           " gff_file: ", basename(gff_file),
                           " - Longest Open Reading Frame")
  writeLines(c(protein_header, protein),
             file.path(protein_output_dir_path, 
                       paste0(gene_name, "_protein.fasta")))
  
  # Write transcript output with header
  transcript_header <- paste0(">", gene_name,
                              " Gene hit: ", basename(gene_hit_path),
                              " Genome: ", basename(genome_file),
                              " gff_file: ", basename(gff_file),
                              " - Full Transcript Sequence")
  writeLines(c(transcript_header, as.character(full_transcript)), 
             file.path(transcript_output_dir_path,
                       paste0(gene_name, "_transcript.fasta")))
  
}

# main --------------------------------------------------------------------
main <- function(genome_dir, 
                 query_dir, 
                 gff_dir, 
                 output_dir, 
                 evalue = 1e-5, 
                 min_seq_len = 50000) {
  
  # read in an index file which should contain gff_file_path with matching genome_file_path
  if (!(file.exists("index.csv"))) {
    stop("The file index.csv does not exist.")
  }
  index_from_file <- read_csv(file = "index.csv",
                              show_col_types = FALSE) |>
    na.omit()

  if (!("gff_file_path" %in% colnames(index_from_file))) {
    stop("index.csv must contain a column named 'gff_file_path'")
  }

  # create output dirextory if not exist already
  if (!(dir.exists(output_dir))) {
    dir.create(output_dir, recursive = TRUE)
  }

  # create empty dataframe for the file paths, to save the iterated path to an index
  # dataframe
  index <- data.frame(
    blast_hit_file_path <- character(0),
    gff_file_path = character(0),
    gene_hit_path = character(0),
    stringsAsFactors = FALSE
  )

  # create a blast datavase and run tblastn with the query fasta's against it
  create_blastdb_and_run_tblastn(genome_dir, query_dir, output_dir, evalue)

  # get the file names from the blast hit
  blast_hits <- list.files(path = file.path(output_dir, "blast_hits"),
                           pattern = "\\.tsv$")

  # create the gene_hits directory if not exist already
  gene_hit_dir_path <- file.path(output_dir, "gene_hits")
  if (!dir.exists(gene_hit_dir_path)) {
    dir.create(gene_hit_dir_path, recursive = TRUE)
  }

  # loop over the blast results and read them in
  for (blast_hit in blast_hits) {
    message("Searching gff3 files for annotated genes for blast hit: '", blast_hit, "'...")

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
        gene_hit_file_path <- file.path(gene_hit_dir_path,
                                        paste0(tools::file_path_sans_ext(blast_hit),"_with_",tools::file_path_sans_ext(gff_file_hit),"_gene_hit.csv"))

        write_csv(result,
                  file = gene_hit_file_path)

        # add pathnames to index data frame
        index <- rbind(index, data.frame(
          blast_hit_file_path = blast_hit_file_path,
          gff_file_path = file.path(gff_file_path),
          gene_hit_path = file.path(gene_hit_file_path)
        ))
      } else {
        next
      }

    }
    }

    # save the current indexing in a csv as backup incase the inner_join fails
    write_csv(index, file = "index.csv.bak")
    
    # fullfil path
    index_from_file$gff_file_path <- file.path(gff_dir, index_from_file$gff_file_path)
    index_from_file$genome_file_path <- file.path(genome_dir, index_from_file$genome_file_path)


    # join the index's
    index_de <- left_join(index, index_from_file, by = "gff_file_path") |>
      na.omit() |>
      mutate(query_path = str_replace(gff_file_path, "_in_.*", ".fasta"))
    
    # 
    # index_de$gene_hit_path <- file.path(output_dir, index_de$gene_hit_path)
  
    # Translate and transcripe the annotated blast hits to protein sequences
    # gene_hit are the best blast results with the gene name from the gff file
    for (i in seq_along(index_de$gene_hit_path)) {
      cat("Iteration: ", i)
      gene_hit_path <- index_de$gene_hit_path[i]
      gff_file_path <- index_de$gff_file_path[i]
      genome_file_path <- index_de$genome_file_path[i]
      
      # remove output_dir from paths
      gff_file_path <- sub(pattern = paste0("^", output_dir, "[/\\]"),
                           replacement = "", gff_file_path)
      genome_file_path <- sub(pattern = paste0("^", output_dir, "[/\\]"),
                              replacement = "", genome_file_path)
      
      # read in the blast hit
      gene_hit_file_path_df <- read_csv(file = gene_hit_path)
      
      # extract
      transcript_id_raw <- gene_hit_file_path_df$gname # ID=TARGET;
      print(head(transcript_id_raw))
      transcript_ids <- unique(stri_extract_first(transcript_id_raw,
                                            regex = "(?<=ID=)[^;]+")) # ID=TARGET; stops at first ";"
      print(head(transcript_ids))
      
      for (transcript_id in transcript_ids) {
        
        ##
        # Debug: Print values right before the call
        message("\nAbout to call gene_seq_to_prot with these values:")
        message("genome_file = '", genome_file_path, "'")
        message("gff_file = '", gff_file_path, "'")
        message("transcript_id = '", transcript_id, "'")
        message("output_dir = '", output_dir, "'\n")
        ##
        
        # call based on index
        gene_seq_to_prot(genome_file = genome_file_path,
                         gff_file = gff_file_path,
                         transcript_id = transcript_id,
                         output_dir = output_dir,
                         gene_hit_path = gene_hit_path)
        
      }
    }
  # save index
  write_csv(index_de, file = "final_index.csv")
}

#####

# tomato_blast <- function(genome_dir, gen_file_pattern, query_dir, tomato_out_dir) {
#   if (!(dir.exists(tomato_out_dir))) {
#     dir.create(tomato_out_dir)
#   }
#   
#   for (proteine_file in query_dir) {
#     # look if resulting protein sequences can be found in the tomato genome
#     create_blastdb_and_run_tblastn(genome_dir = genome_dir,
#                                    query_dir = query_dir,
#                                    output_dir = tomato_out_dir)
#   }
# }



# calls ---------------------------------------------------------------
main(genome_dir = "genomes",
     query_dir = "meiotic_genes_protein_fasta",
     gff_dir = "gff_files",
     output_dir = "out_2024_02_16",
     evalue = 1e-5)




# system2("makeblastdb",
#         args = c("-in", "tomato/UP000004994_4081.fasta", 
#                  "-dbtype", "prot",
#                  "-out", "tco_out/tomato_blast_hit/tomato_db"))

###
# for (protein_fasta in list.files(path = "tco_out/protein", full.names = T)) {
#   tomato_output_file_name <- paste0("tco_out/tomato_blast_hit_", basename(protein_fasta), ".tsv")
#   # Run Bash blastp! protein sequence to proteom
#   system2("blastp",
#           args = c("-query", protein_fasta,
#                    "-db", "tco_out/tomato_blast_hit/tomato_db",
#                    "-out", tomato_output_file_name,
#                    "-evalue", 1e-5,
#                    "-outfmt", "\"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\""))
# }
# 
# 
# # Run Bash blast**p**! protein sequence to proteom
# system2("blastp",
#         args = c("-query", "tco_out/protein/Soltub.DH_A07.04_2G002120.1_protein.fasta",
#                  "-db", "tco_out/tomato_blast_hit/tomato_db",
#                  "-out", "tco_out/tomato_hit.tsv",
#                  "-evalue", 1e-5,
#                  "-outfmt", "\"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\""))
# 
