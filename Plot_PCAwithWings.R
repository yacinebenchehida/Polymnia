# Define path to PCs csv
data <- read.table("/Users/yacinebenchehida/Desktop/Convergent_evolution/Polymnia/PCA\ data/PCs_rF_S1.csv",sep=",",header = T)

# Load images
png_files <- list.files("/Users/yacinebenchehida/Desktop/Convergent_evolution/Polymnia/photo/dorsal_forewing_cartoons_resized/", pattern = "*.png$", full.names = TRUE) # forewing
#png_files <- list.files("/Users/yacinebenchehida/Desktop/Convergent_evolution/Polymnia/photo/ventral_hindwing_cartoons_resized/", pattern = "*.png$", full.names = TRUE) # hindwing

# Subset individuals to plot
male_samples <- read.table("/Users/yacinebenchehida/Desktop/Convergent_evolution/Polymnia/male_list.txt")
male_samples <- unlist(c(male_samples))
data <- data[data$IID %in% male_samples, ]

# Extract prefix before third underscore
data$prefix <- sub("^(([^_]+_[^_]+)_).*", "\\1", data$IID)

# For each prefix, find the matching path and assign each png image to the right ID
data$url = "tobeadded"
for (i in 1:nrow(data)) {
  match <- png_files[grepl(data$prefix[i], png_files)]
  if (length(match) > 0) {
    data$url[i] <- match[1]   # take the first match
  } else {
    data$url[i] <- NA         # or leave blank if no match
  }
}

# Plot
library(ggplot2)
library(ggimage)

pdf("/Users/yacinebenchehida/Desktop/Convergent_evolution/Polymnia/photo/test_rapid.pdf",10,10)
ggplot(data, aes(PC1, PC1, image = url)) +
  theme_bw() +   geom_image(size = 0.06)
dev.off()
