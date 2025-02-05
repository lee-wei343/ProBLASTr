library(tidyverse)

db_files <- list.files(path = "output_dir", pattern = "*.tbv")
gff_files <- list.files(path = "gff_files", pattern = "*.gff3")

find_matching_genes <- function(db_file, gff_file) {
  # read database file
  db <- read.table(
    file = db_file,
    sep = "\t",
    col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                  "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  ) |>
    mutate(sseqid = as.character(sseqid)) # make sure that ssequid of both df will be of same type
  
  # read gff file
  gff <- read.table(
    file = gff_file,
    sep = "\t",
    col.names = c("sseqid", "source", "type", "gstart", "gend", 
                  "score", "strand", "phase", "gname")
  ) |>
    mutate(sseqid = as.character(sseqid)) |> # make sure that ssequid of both df will be of same type
    filter(type %in% c("mRNA", "transcript"))
  
  # join the database with the gff file where the sequenceid is the same
  # find a match from the sequence in in the gene sequence 
  result <- inner_join(db, gff, by = "sseqid", relationship = "many-to-many") |>
    mutate(seq_in_gene = between(sstart, gstart, gend) & between(send, gstart, gend)) %>%
    filter(seq_in_gene)
  
  return(result)
}

# apply the function to all files from the lists
results <- mapply(
  find_matching_genes, 
  file.path("output_dir", db_files), 
  file.path("gff_files", gff_files), 
  SIMPLIFY = FALSE
)
