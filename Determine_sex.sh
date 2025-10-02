CHR1="OZ240203.1" # Chr Z
CHR2="OZ240191.1" # Chr 12
VCF="/mnt/scratch/projects/biol-specgen-2018/yacine/Polymnia/Data/polymnia.filtered_chr_maf10.snps.vcf.gz"

# Compute average DP per sample for both chromosomes and ratio
paste \
  <(bcftools query -l $VCF) \
  <(bcftools query -f '[%DP\t]\n' -r $CHR1 $VCF | awk '
    {for(i=1;i<=NF;i++){sum[i]+=$i; count[i]+=($i!=".")}} 
    END {for(i=1;i<=NF;i++) print sum[i]/count[i]}') \
  <(bcftools query -f '[%DP\t]\n' -r $CHR2 $VCF | awk '
    {for(i=1;i<=NF;i++){sum[i]+=$i; count[i]+=($i!=".")}} 
    END {for(i=1;i<=NF;i++) print sum[i]/count[i]}') \
  | awk '{printf "%s\t%.2f\t%.2f\t%.2f\n", $1, $2, $3, $2/$3}' \
> coverage_ratio_per_individual.txt

# Assign individuals with a ratio < 0.9 as females otherwise male
cat coverage_ratio_per_individual.txt |awk '{ if ($4 < 0.9) $4="F"; else $4="M"; print }' > ../Data/ind_with_sex.txt
