#! /bin/bash

#SBATCH --mem=20GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --job-name=GWAS_GEMMA
#SBATCH --time=0-01:00:00

module load R/4.2.1-foss-2022a
module load BCFtools/1.19-GCC-13.2.0

####################
# Set useful paths #
####################
RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Polymnia/Results/female_only/$1/$2"
PLINK="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/plink_linux_x86_64_20231018/plink"
GEMMA="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/gemma/gemma-0.98.5-linux-static-AMD64"
VCF="/mnt/scratch/projects/biol-specgen-2018/yacine/Polymnia/Data/females_polymnia.filtered_chr_maf5.snps.vcf.gz"
PHENOTYPES="${RESULTS}/GEMMA_encoding_phenotype_${2}.txt"
CSV="/mnt/scratch/projects/biol-specgen-2018/yacine/Polymnia/Data/$3"
PC="$4"

#  Set working directory
mkdir -p $RESULTS

# Create phenotype file
python3 ./Prepare_phenotype_files.py --file1 <(bcftools query -l $VCF) --file2 $CSV --col $PC > $PHENOTYPES

# Generate bed files with plink 
 $PLINK --vcf $VCF --double-id --allow-extra-chr --allow-no-sex --set-missing-var-ids @:# --pheno $PHENOTYPES --make-bed --pca --out $RESULTS/$2 --threads 8
 cat $RESULTS/"$2".fam|perl -pe 's/2$/0/g'|perl -pe 's/2.5$/0.5/g'|perl -pe 's/3$/1/g' > $RESULTS/tmp.fam
 mv $RESULTS/tmp.fam $RESULTS/"$2".fam

# Get pairwise relatedness between individuals
 $GEMMA -bfile $RESULTS/$2 -gk 1 -o "$2".relatedness -miss 0.2 -maf 0.05  -outdir $RESULTS > $RESULTS/gemma_"$2"_relatedness2.out

# GWAS 
 $GEMMA -bfile $RESULTS/$2 -lmm 4 -o "$2".assoc.gemma -miss 0.2 -maf 0.05  -outdir $RESULTS -k $RESULTS/"$2".relatedness.cXX.txt > $RESULTS/"$2".gemma.out

# Create output for R Manhattan plots
 cat  $RESULTS/*.assoc.gemma.assoc.txt| awk '$14 < 0.05 {print $1"\t"$2"\t"$3"\t"$14}'|grep -v "##" > $RESULTS/plot_"$2".txt

# Create file with scaffold sorted by ascending size
bcftools index -s $VCF|awk '{print $1"\t"$2}'|sort -V -rk 2 > $RESULTS/scaffold_order_$2.txt

# Get the number of SNPs used for the analyses (used to define the bonferroni threshold)
THRESHOLD=$(cat $RESULTS/*assoc.gemma.log.txt|grep "analyzed SNPs/var" |awk '{print $7}')

# Make a png plot of the gwas (low quality)
Rscript ./Plot.R $RESULTS/"$2".assoc.gemma.assoc.txt $RESULTS/scaffold_order_$2.txt $THRESHOLD $RESULTS $PHENOTYPES $VCF
mv $RESULTS/GWAS.png  $RESULTS/GWAS_"$2".png

# sbatch ./master_script.sh bF bF_S1_PC1 PCs_bF_S1.csv 3
