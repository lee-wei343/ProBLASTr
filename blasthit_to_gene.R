# TCO
library(tidyverse)

min_seq_len = 50000 # minimum sequence length

# get the file names from the blast hit
db_files <- list.files(path = "output_dir", pattern = "*.tsv")

# get the file names from the annotated gff3 files
gff_files <- list.files(path = "gff_files", pattern = "*.gff3")

# Create output directory if it doesn't exist
dir.create("output_dir", showWarnings = FALSE)

#
for (db_file in db_files) {
  for (gff_file in gff_files) {
    print(db_file)
    print(gff_file)
    
    
    # read database file
    blast_hits <- read.table(
      file = file.path("output_dir", db_file),
      sep = "\t",
      col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                    "qstart", "qend", "sstart", "send", "evalue", "bitscore")
    ) |>
      mutate(sseqid = as.character(sseqid)) # make sure that ssequid of both df will be of same type
    
    # Filter overlaps only on the same chromosome
    filtered_blast <- blast_hits |>
      group_by(sseqid) |>
      arrange(sstart) |>
      mutate(
        prev_end = lag(send, default = first(send)), # find the previous sequence end
        distance_to_prev = sstart - prev_end, # calculate the length of a sequence
        is_overlap = sstart <= prev_end # filter out sequence starts in a previous sequence
      ) |> 
      filter((distance_to_prev <= min_seq_len | is.na(distance_to_prev)) & !is_overlap) |>
      filter(bitscore >= 50) |> # only keep high bitscore
      ungroup() |>
      select(-prev_end, -distance_to_prev, -is_overlap)
    
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
    result <- inner_join(filtered_blast, gff, by = "sseqid", relationship = "many-to-many") |>
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
                    tools::file_path_sans_ext(db_file), 
                    "_with_",
                    tools::file_path_sans_ext(gff_file), 
                    "_blast_matches.csv")))
      
    } else {
      next
    }
  }
  
}
