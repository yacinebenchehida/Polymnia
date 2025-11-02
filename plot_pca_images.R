plot_pca_images <- function(
    pca_file,
    image_dir,
    x_pc = "PC1",
    y_pc = "PC2",
    subset_samples = NULL,
    output_pdf = "pca_plot.pdf",
    single_PC = FALSE,
    subsetting = FALSE,
    remove_bg = TRUE,
    bg_color = "white",
    fuzz = 10,
    trim = TRUE,
    temp_dir = file.path(tempdir(), "pca_imgs")) {
  
  # Load required libraries
  suppressPackageStartupMessages({
    library(ggplot2)
    library(ggimage)
  })
  if (remove_bg) {
    if (!requireNamespace("magick", quietly = TRUE))
      stop("Please install the 'magick' package: install.packages('magick')")
  }
  
  # --- Clean up any previous runs to avoid caching issues ---
  if (dir.exists(temp_dir)) unlink(temp_dir, recursive = TRUE, force = TRUE)
  dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Clear ggimage’s internal cache (to force reload of all PNGs)
  if (exists(".image_cache", envir = asNamespace("ggimage"))) {
    cache_env <- get(".image_cache", envir = asNamespace("ggimage"))
    rm(list = ls(envir = cache_env), envir = cache_env)
  }
  gc()  # clear any cached magick objects
  
  # --- Load PCA data ---
  data <- read.table(pca_file, sep = ",", header = TRUE, stringsAsFactors = FALSE)
  
  # Optional subsetting
  if (!is.null(subset_samples)) {
    if (file.exists(subset_samples)) {
      subset_vec <- unlist(read.table(subset_samples, stringsAsFactors = FALSE))
    } else {
      subset_vec <- subset_samples
    }
    data <- data[data$IID %in% subset_vec, , drop = FALSE]
  }
  
  # --- Match PNG images ---
  png_files <- list.files(image_dir, pattern = "\\.png$", full.names = TRUE, ignore.case = TRUE)
  data$url <- NA_character_
  
  for (i in seq_len(nrow(data))) {
    matches <- png_files[grepl(data$my_names[i], png_files, fixed = TRUE)]
    if (length(matches) > 0) {
      matched_path <- matches[1]
      
      # --- Optional background removal ---
      if (remove_bg) {
        img <- magick::image_read(matched_path)
        if (trim) img <- magick::image_trim(img, fuzz = fuzz)
        img <- magick::image_transparent(img, color = bg_color, fuzz = fuzz)
        img <- magick::image_background(img, color = "none")
        
        outfile <- file.path(temp_dir, paste0(i, "_", basename(matched_path)))
        magick::image_write(img, path = outfile, format = "png")
        data$url[i] <- outfile
      } else {
        data$url[i] <- matched_path
      }
    } else {
      data$url[i] <- NA_character_
    }
  }
  
  reordered_data <- data[order(data[[x_pc]]), ]
  reordered_data$order <- seq_len(nrow(reordered_data))
  
  if(subsetting==TRUE){
    reordered_data <- reordered_data[seq(1, nrow(reordered_data), by = 4), ]
  }
  
  
  
  # --- Plot using cairo_pdf for proper transparency ---
  if(single_PC==TRUE){
  if (capabilities("cairo")) {
    grDevices::cairo_pdf(output_pdf, width = 10, height = 10)
  } else {
    pdf(output_pdf, width = 10, height = 10)
  }
  size = 1.8 / nrow(reordered_data)
  p_one_pc <- ggplot(reordered_data, aes(x = order, y = order, image = url)) +
    theme_bw() + ylab(paste(x_pc)) +
    geom_image(size = size, by = "width", na.rm = TRUE)
  print(p_one_pc)
  
  dev.off()
  }else{

  # --- Plot using cairo_pdf for proper transparency ---
  if (capabilities("cairo")) {
    grDevices::cairo_pdf(output_pdf, width = 10, height = 10)
  } else {
    pdf(output_pdf, width = 10, height = 10)
  }
  
  p <- ggplot(data, aes_string(x = x_pc, y = y_pc, image = "url")) +
    theme_bw() +
    geom_image(size = 0.05, by = "width", na.rm = TRUE)
  
  print(p)
  dev.off()
  } 
  # Example debug output (adjust as needed)
  #print(data[data$PC1 > 0 & data$PC1 < 12 & data$PC2 < -12 & data$PC2 > -25, ])
  
  invisible(data)
}

 # Usage example:   plot_pca_images(pca_file = "PCA data/PCs_yH_S1.csv", image_dir="photo/ventral_hindwing_cartoons_resized/", x_pc = "PC1", y_pc = "PC2", subset_samples = NULL, output_pdf = "/Users/yacinebenchehida/Desktop/Convergent_evolution/Polymnia/Results/No_Pop_structure/yH/yH_S1_PC1/yH_S1_PC1_single.pdf", remove_bg = TRUE, bg_color = "white", fuzz = 10, trim = TRUE, single_PC=TRUE, subsetting=TRUE)
