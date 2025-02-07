#!/bin/bash

set -euo pipefail
# set -x # for debugging

QUERY_DIR=meiotic_genes_protein_fasta
DB_DIR=genomes
OUTPUT_DIR=output_dir

mkdir -p $OUTPUT_DIR

# the file ASY3_mafft.fasta is in the mafft format and results in the error: "BLAST query error: CFastaReader: Near line 2, there's a line that doesn't look like plausible data, but it's not marked as defline or comment."
# clean mafft files by removing the "-"
# for mafft_file in $QUERY_DIR/*mafft*; do
#     output_file="$QUERY_DIR/$(basename "$mafft_file" .fasta)_clean.fasta"
#     sed '/^>/!s/-//g' "$mafft_file" > "$output_file"
# done



# make a blast nucleotide database for each genome fasta file
for haplotype_fasta_file in $DB_DIR/*\.fasta; do
    echo "Processing ${haplotype_fasta_file}"
    makeblastdb \
        -in "$haplotype_fasta_file" \
        -dbtype nucl \
        -title "potato_genome_db" \
        -out "$OUTPUT_DIR/$(basename "$haplotype_fasta_file" .fasta)_db" \
        -parse_seqids
done

# loop through each protein target_seq and tblastn it against all the databases
for query in $QUERY_DIR/*\.fasta; do
    query_name=$(basename $query .fasta)

    for db in $OUTPUT_DIR/*_db.nin; do
        db_name=$(basename "$db" .nin)
        echo "Processing $query_name in database $db_name"

        tblastn \
            -query "$query" \
            -db "$OUTPUT_DIR/$db_name" \
            -out "$OUTPUT_DIR/${query_name}_in_${db_name}.tsv" \
            -evalue 1e-5 \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    done
done

