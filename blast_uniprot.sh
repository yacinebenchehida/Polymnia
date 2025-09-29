REF="/mnt/scratch/projects/biol-specgen-2018/yacine/Polymnia/Data/Gene_color_genes"
FASTA_FILE="GCA_965197365.1_ilMecPoly1.hap1.1_genomic.fna"
RESULTS="Polymnia_$1_$2_$3"

mkdir -p $RESULTS

(echo ">$1"; cat $REF/$FASTA_FILE|awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}'|grep -A1 $1|grep -v $1| cut -c "$2"-"$3") > $RESULTS/"$1"_"$2"_"$3".fasta

echo FASTA FILE READY

module load AUGUSTUS/3.5.0-foss-2023a
augustus --progress=true --strand=both --species=heliconius_melpomene1 $RESULTS/"$1"_"$2"_"$3".fasta > $RESULTS/"$1"_"$2"_"$3".gff

echo AUGUSTUS FINISHED RUNNING 

getAnnoFasta.pl  $RESULTS/"$1"_"$2"_"$3".gff
echo PROTEIN FASTA READY

module purge 
module load BLAST+/2.14.1-gompi-2023a

echo STARTING BLASTP
blastp -query $RESULTS/"$1"_"$2"_"$3".aa -db database_insect -outfmt 6 -max_target_seqs 5 | awk '$12 > 80' > $RESULTS/"Polymnia_$1_$2_$3"_blastp_gff_genes_against_uniprot.txt
echo BLASTP DONE




# USAGE = ./blast_uniprot.sh  ctg002890_1 1030001 1260001
