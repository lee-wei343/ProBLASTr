#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

library(GenomicRanges)
library(Biostrings)
library(GenomicFeatures)
library(rtracklayer)
library(Rsamtools)

# Set working directory
setwd("~/Documents/potato_DNA/protein_extraction_test")

# Files
genome_file <- "DH_A07_v1.asm.sm.fasta"
gff_file <- "gene_exons.gff3"

# Read gff file
gff <- import(gff_file)

# Filter for the gene and for exons
gff[gff$type == "exon"]
exons <- gff[as.character(gff$Parent) == "Soltub.DH_A07.04_2G002120.1"]

# Sort exons based on strand
if (unique(strand(exons)) == "-") {
  exons <- sort(exons, decreasing = TRUE)
} else {
  exons <- sort(exons)
}

# Select genome file to extract the sequence further ahead
fa <- FaFile(genome_file)

# Extract sequences for each exon
exon_seqs <- sapply(seq_along(exons), function(i) {
  chr <- as.character(seqnames(exons[i]))
  start <- start(exons[i])
  end <- end(exons[i])
  
  # Extract sequence from the genome
  seq <- scanFa(fa, param = GRanges(seqnames = chr, 
                                    ranges = IRanges(start = start, end = end)))
  
  # Reverse complement if on minus strand
  if (as.character(strand(exons[i])) == "-") {
    seq <- reverseComplement(seq)
  }
  
  return(as.character(seq))
})

# Join exon sequences
full_transcript <- DNAString(paste(exon_seqs, collapse=""))

# Function to find and translate the longest Open Reading Frame
findLongestORF <- function(sequence) {
  # Find all ATG positions
  atg_positions <- start(matchPattern("ATG", sequence))
  longest_protein <- ""
  
  for (start_pos in atg_positions) {
    # Get subsequence starting from ATG
    current_seq <- subseq(sequence, start_pos)
    
    # Ensure sequence length is a multiple of 3 before translating
    remainder <- length(current_seq) %% 3
    if (remainder != 0) {
      current_seq <- subseq(current_seq, 1, length(current_seq) - remainder)
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
protein <- findLongestORF(full_transcript)

# Get the gene name from the exons object
gene_name <- unique(as.character(mcols(exons)$Parent))

# Write protein output with header
protein_header <- paste0(">", gene_name, " - Longest Open Reading Frame")
writeLines(c(protein_header, protein), "protein.fasta")

# Write transcript output with header
transcript_header <- paste0(">", gene_name, " - Full Transcript Sequence")
writeLines(c(transcript_header, as.character(full_transcript)), "transcript.fasta")