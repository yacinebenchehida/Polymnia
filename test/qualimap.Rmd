```{r setup, include=FALSE}
remotes::install_github('yihui/knitr',force=TRUE)
library(knitr)
library(rmarkdown)
opts_chunk$set(echo = TRUE)
```

# I) Qualimap
## Running qualimap
```{bash}
#!/bin/bash
# Author: Yacine Ben Chehida

#SBATCH --time=25:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=15G 
#SBATCH --account=BIOL-SPECGEN-2018 
#SBATCH --job-name=map_sort_RG_dup_merge

##################################
# 1 - Load necessar softwares #
##################################
module load bio/Qualimap/2.2.1-foss-2019b-R-3.6.2

###########################################                              
# 2 - check mapping quality with qualimap #
###########################################
mkdir -p $path_results/2_Qualimap/"$1"

qualimap bamqc --java-mem-size=20G \
-bam $path_results/1_sorted_dedup_bam/$1/"$1"_sorted_dedup.bam \
-c \
-gd HUMAN \
-nt 8 \
-outdir $path_results/2_Qualimap/"$1" \
-outformat HTML  
```

The reference genome used for each sample can be found below:
```{r}
setwd("/Users/yacinebenchehida/Desktop/porpoise_stuff/2025-2026")
dat = read.table("untitled folder/species.txt", comment.char = "")
colnames(dat) = c("ID","Species","Couleur")
dat
```
