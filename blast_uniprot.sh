#!/bin/bash

CHR=$1
START=$2
END=$3
OUTPUT_DIR=$4
DB_PATH=$5

REF="/mnt/scratch/projects/biol-specgen-2018/yacine/Polymnia/Data/Gene_color_genes"
FASTA_FILE="GCA_965197365.1_ilMecPoly1.hap1.1_genomic.fna"

mkdir -p "$OUTPUT_DIR"

# Extract sequence region
(echo ">$CHR"; 
 cat $REF/$FASTA_FILE | awk '/^>/ { if(NR>1) print ""; printf("%s\n",$0); next; } { printf("%s",$0);} END {printf("\n"); }' \
 | grep -A1 "$CHR" | grep -v "$CHR" | cut -c "$START"-"$END") > "$OUTPUT_DIR/${CHR}_${START}_${END}.fasta"

echo "FASTA FILE READY"

# Run Augustus
module load AUGUSTUS/3.5.0-foss-2023a
augustus --progress=true --strand=both --species=heliconius_melpomene1 \
  "$OUTPUT_DIR/${CHR}_${START}_${END}.fasta" > "$OUTPUT_DIR/${CHR}_${START}_${END}.gff"

echo "AUGUSTUS FINISHED RUNNING"

# Extract protein fasta
getAnnoFasta.pl "$OUTPUT_DIR/${CHR}_${START}_${END}.gff"
echo "PROTEIN FASTA READY"

# Run BLAST
module purge
module load BLAST+/2.14.1-gompi-2023a

echo "STARTING BLASTP"
blastp -query "$OUTPUT_DIR/${CHR}_${START}_${END}.aa" \
  -db "$DB_PATH" \
  -outfmt 6 -max_target_seqs 5 | awk '$12 > 80' > "$OUTPUT_DIR/Polymnia_${CHR}_${START}_${END}_blastp_results.txt"

echo "BLASTP DONE"
