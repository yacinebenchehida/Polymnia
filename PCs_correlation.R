library(corrplot)
library(csv)
library(tidyverse)

gca_files <- list.files("/Users/yacinebenchehida/Downloads/PCA\ data/", pattern = "*.csv$", full.names = TRUE)
dfs <- list()
id_lists <- list()
Number_of_PCs_to_keep = 20
for (i in seq_along(gca_files)) {
  data <- read.csv(gca_files[i])
  data <-data[,c(1:Number_of_PCs_to_keep)]
  
  print(dim(data))
  # standardize ID column name
  colnames(data)[2] <- "IID"
  
  # rename columns except IID
  comp <- gsub(".*PCs_([^.]*)\\.csv", "\\1", gca_files[i])
  colnames(data)[colnames(data) != "IID"] <- 
    paste0(colnames(data)[colnames(data) != "IID"], "_", comp)
  
  dfs[[i]] <- data
  id_lists[[i]] <- data$IID
}


# Find IDs present in all files
common_ids <- Reduce(intersect, id_lists)

# Keep only rows with common IDs in each df
dfs_common <- lapply(dfs, function(df) df[df$IID %in% common_ids, ])

# Remove IID column from all except the first before cbind
dfs_common <- lapply(dfs_common, function(df) df[, -c(1,2), drop = FALSE])

# Combine
final <- do.call(cbind, dfs_common)
final <- as.data.frame(lapply(final, as.numeric))
corr_matrix <- cor(final, use = "everything")
#corrplot(corr_matrix,method = 'color', tl.cex = 0.2)
#corrplot(corr_matrix,method = 'color',diag = FALSE, type = 'upper', tl.cex = 0.2)
#method = 'color', diag = FALSE, type = 'upper'
#corrplot(corr_matrix, order="hclust",method = 'color')



block_size <- dim(data)[2]-2
n <- ncol(corr_matrix)

# plot without re-ordering
pdf("/Users/yacinebenchehida/Desktop/Convergent_evolution/Polymnia/corr_PCs.pdf",50,50)
corrplot(corr_matrix, method = 'color', tl.cex = 0.6, tl.srt = 45, tl.offset=0.6)


nblocks <- n / block_size  # number of CSVs

# loop over all pairs of CSV blocks
for (a in 0:(nblocks-1)) {
  for (b in 0:(nblocks-1)) {
    
    i <- a * block_size + 1
    j <- i + block_size - 1
    
    k <- b * block_size + 1
    l <- k + block_size - 1
    
    # x spans columns (b block)
    x1 <- k - 0.5
    x2 <- l + 0.5
    
    # y spans rows (a block, inverted for corrplot coords)
    y1 <- n - j + 0.5
    y2 <- n - i + 1.5
    
    rect(x1, y1, x2, y2, border = "black", lwd = 2)
  }
}


# draw diagonal blocks of block_size
for (i in seq(1, n, by = block_size)) {
  j <- min(i + block_size - 1, n)
  
  # x-coords (left->right)
  x1 <- i - 0.5
  x2 <- j + 0.5
  
  # y-coords must be inverted for corrplot plotting orientation
  y1 <- n - j + 0.5
  y2 <- n - i + 1.5
  
  rect(x1, y1, x2, y2, border = "red", lwd = 2)
}
dev.off()
