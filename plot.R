##############################
# Install required libraries #
##############################
library(ggfortify)
library(qqman)
library(tidyverse)
library(ggtext)
library(normentR)
library(grid)
library(gridExtra)
library(dplyr)
library(sjlabelled)
library(ggplot2)
library(cowplot)
library(quotidieR)
library(ggthemes)
library(ggridges)

###############
# Import data #
###############
dat = commandArgs(trailingOnly=TRUE)
print(-log10(0.05/as.numeric(dat[3])))
gwas = read.table(dat[1],skip=1)
bonf_threshold = -log10(0.05/as.numeric(dat[3]))
path_results=dat[4]
pheno_file = dat[5]
vcf_file = dat[6]
print(path_results)

################          
# Prepare data #
################
gwas = gwas[,c(1,2,3,14)]
colnames(gwas) = c("CHR","SNP","BP","P")
order_contigs = read.table(dat[2])
order_contigs = as.data.frame(order_contigs)
order_contigs = as.factor(order_contigs$V1)

gwas$CHR=factor(gwas$CHR, levels=order_contigs)

sig_data <- gwas %>% 
  subset(P < 1)
notsig_data <- gwas %>%
  subset(P >= 1) %>% 
  group_by(CHR) %>% 
  sample_frac(1)
gwas_data <- bind_rows(sig_data, notsig_data)

data_cum <- gwas_data %>% 
  group_by(CHR) %>% 
  summarise(max_bp = max(BP)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(CHR, bp_add)

gwas_data <- gwas_data %>% 
  inner_join(data_cum, by = "CHR") %>% 
  mutate(bp_cum = BP + bp_add)

axis_set <- gwas_data %>% 
  group_by(CHR) %>% 
  summarize(center = mean(bp_cum))

ylim <- gwas_data %>% 
  filter(P == min(P)) %>% 
  mutate(ylim = abs(floor(log10(P))) + 2) %>% 
  pull(ylim)


########
# Plot #
########
png(file=paste(path_results,"/GWAS.png",sep=""),width=800,height=500,type="cairo")
ggplot(gwas_data, aes(x = bp_cum, y = -log10(P), 
                                  color =CHR)) +
  geom_point(alpha = 0.75) +
  scale_x_continuous() +
  scale_color_manual(values = rep(c("blue4", "orange3"), unique(length(axis_set$CHR)))) +
  labs(x = NULL, 
       y = "-log<sub>10</sub>(p)") + 
  theme_bw() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(size=10),
    plot.title = element_text(hjust = 0.5,size=10)) +
    geom_hline(yintercept=bonf_threshold, linetype="dashed", 
                color = "red", size=0.8)
dev.off()

###################
# Function: Zoomed Manhattan around top SNPs #
###################
plot_top_snps_zoom <- function(gwas_df, threshold, zoom_kb = 50, path_results) {
  
  plots_list <- list()
  
  # Find top SNP per scaffold above threshold
  top_snps <- gwas_df %>%
    filter(-log10(P) >= threshold) %>%
    group_by(CHR) %>%
    filter(P == min(P)) %>% 
    slice_sample(n = 1) %>%  # pick one if ties
    ungroup()
  
  for (i in 1:nrow(top_snps)) {
    snp <- top_snps[i, ]
    chr_data <- gwas_df %>% filter(CHR == snp$CHR)
    max_bp <- max(chr_data$BP)
    
    # Define zoom window
    if (snp$BP < zoom_kb*1000) {
      start_bp <- 0
      end_bp <- min(2*zoom_kb*1000, max_bp)
    } else if (snp$BP > (max_bp - zoom_kb*1000)) {
      start_bp <- max(0, max_bp - 2*zoom_kb*1000)
      end_bp <- max_bp
    } else {
      start_bp <- snp$BP - zoom_kb*1000
      end_bp <- snp$BP + zoom_kb*1000
    }
    
    chr_window <- chr_data %>% filter(BP >= start_bp & BP <= end_bp)
    
    p <- ggplot(chr_window, aes(x = BP, y = -log10(P))) +
      geom_point(alpha = 0.75, color = "blue4") +
      geom_hline(yintercept = threshold, linetype = "dashed", color = "red", size = 0.8) +
      labs(title = paste("Zoom:", snp$CHR, "Top SNP:", snp$SNP),
           x = "Position (bp)", y = "-log10(P)") +
      theme_bw() +
      theme(plot.title = element_text(size = 10, hjust = 0.5))
    
    plots_list[[i]] <- p
  }
  
  # Combine all zoomed plots in a grid
  n_col <- ceiling(sqrt(length(plots_list)))
  n_row <- ceiling(length(plots_list)/n_col)
  
  png(file = paste0(path_results, "/GWAS_zoom.png"), width = 1200, height = 800, type = "cairo")
  grid.arrange(grobs = plots_list, ncol = n_col, nrow = n_row)
  dev.off()
  
  return(top_snps)
}

###################
# Function: Genotype-Phenotype plots for top SNPs #
###################
plot_genotype_phenotype <- function(top_snps, pheno_file, vcf_file, path_results) {

  pdf(file = paste0(path_results, "/TopSNPs_genotype_phenotype.pdf"))
  
  for (i in 1:nrow(top_snps)) {
    snp <- top_snps[i, ]
    region <- paste0(snp$CHR, ":", snp$BP)
    
   tmp_gt <- file.path(tempdir(), paste0("gt_", snp$SNP, ".txt"))

    # Extract genotypes into a temporary file
    cmd_gt <- paste0(
      "bcftools view --regions ", region, " ", vcf_file,
      " | bcftools query -f '%POS [ %GT]\\n' | perl -pe 's/ /\\n/g' | awk 'NR>2' > ", tmp_gt
    )
    system(cmd_gt)
    
    # Read phenotype file (columns FID IID PHENO)
    pheno <- read.table(pheno_file, stringsAsFactors = FALSE)
    gt <- read.table(tmp_gt, stringsAsFactors = FALSE)
    print(gt)
    if(nrow(gt) != nrow(pheno)){
      warning("Number of genotypes does not match number of phenotypes for SNP: ", snp$SNP)
    }
    
    df <- data.frame(FID = pheno$V1,
                     IID = pheno$V2,
                     Genotype = gt$V1,
                     Phenotype_value = pheno$V3,
                     stringsAsFactors = FALSE)
    
    df$Genotype <- gsub("\\|","/", df$Genotype)
    df <- df %>% filter(Genotype %in% c("0/0","0/1","1/1"))
    df$Genotype_cat <- paste("Top SNP:", snp$SNP)
    
    # Ridge plot
    p1 <- ggplot(df, aes(x = Phenotype_value, y = Genotype, fill = Genotype)) +
      geom_density_ridges(alpha = 0.7, scale = 20) +
      labs(title = paste("Top SNP:", snp$SNP, "in", snp$CHR),
           x = "Phenotype", y = "Genotype") +
      theme_bw()
    print(p1)
    
    # Boxplot
    p2 <- ggplot(df, aes(x = Genotype, y = Phenotype_value, fill = Genotype)) +
      geom_boxplot(outlier.shape = NA) +
      geom_point(position = position_jitter(width = 0.2), size = 2) +
      labs(title = paste("Top SNP:", snp$SNP, "in", snp$CHR)) +
      theme_bw()
    print(p2)
    
    # Clean up
    unlink(tmp_gt)
  }
  
  dev.off()
}

###################
# Run functions
###################
top_snps <- plot_top_snps_zoom(gwas, bonf_threshold, zoom_kb = 100, path_results)
plot_genotype_phenotype(top_snps, pheno_file, vcf_file, path_results)
