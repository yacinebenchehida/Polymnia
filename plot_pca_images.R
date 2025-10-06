plot_pca_images <- function(pca_file,
                            image_dir,
                            x_pc = "PC1",
                            y_pc = "PC2",
                            subset_samples = NULL,
                            output_pdf = "pca_plot.pdf") {
  # Load PCA data
  data <- read.table(pca_file, sep = ",", header = TRUE, stringsAsFactors = FALSE)
  
  # Optional subset of samples
  if (!is.null(subset_samples)) {
    if (file.exists(subset_samples)) {
      subset_vec <- unlist(read.table(subset_samples, stringsAsFactors = FALSE))
    } else {
      subset_vec <- subset_samples
    }
    data <- data[data$IID %in% subset_vec, , drop = FALSE]
  }
  
  # List PNG files
  png_files <- list.files(image_dir, pattern = "*.png$", full.names = TRUE)
  
  # Initialize url column
  data$url <- NA_character_
  
  # Match images using 'my_names' column
  for (i in 1:nrow(data)) {
    match <- png_files[grepl(data$my_names[i], png_files)]
    if (length(match) > 0) {
      data$url[i] <- match[1]  # take the first match
    } else {
      data$url[i] <- NA         # leave NA if no match
    }
  }
  
  # Plot
  library(ggplot2)
  library(ggimage)
  
  pdf(output_pdf, 10, 10)
  print(
    ggplot(data, aes_string(x = x_pc, y = y_pc, image = "url")) +
      theme_bw() +
      geom_image(size = 0.04)
  )
  dev.off()
  
  invisible(data)
}

# Usage example:
plot_pca_images(pca_file="/Users/yacinebenchehida/Desktop/Convergent_evolution/Polymnia/PCA\ data/PCs_yH_S1.csv",
                image_dir="/Users/yacinebenchehida/Desktop/Convergent_evolution/Polymnia/photo/ventral_hindwing_cartoons_resized",
                x_pc = "PC2",
                y_pc = "PC3",
                subset_samples = NULL,
                output_pdf = "pca_plot.pdf")
