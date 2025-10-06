plot_pca_images <- function(pca_file,
                            image_dir,
                            x_pc = "PC1",
                            y_pc = "PC2",
                            subset_samples = NULL,
                            output_pdf = "pca_plot.pdf",
                            ## new args:
                            remove_bg = TRUE,      # set FALSE to keep originals
                            bg_color = "white",    # color to make transparent
                            fuzz = 10,             # 0..100 percent tolerance (start ~10)
                            trim = TRUE,           # trim border before/after making transparent
                            temp_dir = file.path(tempdir(), "pca_imgs")) {
  
  # read PCA
  data <- read.table(pca_file, sep = ",", header = TRUE, stringsAsFactors = FALSE)
  
  # optional subset (keeps your existing behavior)
  if (!is.null(subset_samples)) {
    if (file.exists(subset_samples)) {
      subset_vec <- unlist(read.table(subset_samples, stringsAsFactors = FALSE))
    } else {
      subset_vec <- subset_samples
    }
    data <- data[data$IID %in% subset_vec, , drop = FALSE]
  }
  
  # list PNGs (use regex that actually matches .png)
  png_files <- list.files(image_dir, pattern = "\\.png$", full.names = TRUE, ignore.case = TRUE)
  
  # prepare output column
  data$url <- NA_character_
  
  if (remove_bg) {
    if (!requireNamespace("magick", quietly = TRUE)) {
      stop("Please install the 'magick' package (install.packages('magick')) to remove backgrounds.")
    }
    dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  # match and (optionally) process images
  for (i in seq_len(nrow(data))) {
    # match using my_names (same approach as your original). Use fixed = TRUE to avoid regex surprises:
    matches <- png_files[grepl(data$my_names[i], png_files, fixed = TRUE)]
    
    if (length(matches) > 0) {
      matched_path <- matches[1]
      
      if (remove_bg) {
        img <- magick::image_read(matched_path)
        
        # optional trim (removes edges equal to the background color)
        if (trim) {
          img <- magick::image_trim(img, fuzz = fuzz)
        }
        
        # set pixels matching bg_color to transparent; fuzz is 0..100
        img <- magick::image_transparent(img, color = bg_color, fuzz = fuzz)
        
        # ensure background is 'none' so alpha is preserved when we write
        img <- magick::image_background(img, color = "none")
        
        outfile <- file.path(temp_dir, paste0(sprintf("%04d", i), "_", basename(matched_path)))
        magick::image_write(img, path = outfile, format = "png")
        data$url[i] <- outfile
      } else {
        data$url[i] <- matched_path
      }
    } else {
      data$url[i] <- NA_character_
    }
  }
  
  # plotting: use cairo_pdf if available so transparency is preserved in PDF outputs
  if (grepl("\\.pdf$", output_pdf, ignore.case = TRUE)) {
    if (capabilities("cairo")) {
      grDevices::cairo_pdf(output_pdf, width = 10, height = 10)
    } else {
      warning("cairo not available on this R build — falling back to pdf(). Transparency may not be preserved.")
      pdf(output_pdf, width = 10, height = 10)
    }
  } else {
    # fallback: a normal pdf() or other device depending on extension could be used;
    # for simplicity, use pdf() if the filename doesn't end with .pdf
    pdf(output_pdf, width = 10, height = 10)
  }
  
  # plotting with ggimage
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("ggimage", quietly = TRUE)) {
    dev.off()
    stop("Please install 'ggplot2' and 'ggimage' (install.packages(c('ggplot2','ggimage'))).")
  }
  
  library(ggplot2)
  library(ggimage)
  
  p <- ggplot(data, aes_string(x = x_pc, y = y_pc, image = "url")) +
    theme_bw() +
    geom_image(size = 0.04, by = "width", na.rm = TRUE)
  
  print(p)
  dev.off()
  
  invisible(data)
}



plot_pca_images(
  pca_file="/Users/yacinebenchehida/Desktop/Convergent_evolution/Polymnia/PCA\ data/PCs_yF_S1.csv",
  image_dir="/Users/yacinebenchehida/Desktop/Convergent_evolution/Polymnia/photo/dorsal_forewing_cartoons_resized//",
  x_pc = "PC1",
  y_pc = "PC2",
  remove_bg = TRUE,    # enable background removal
  bg_color = "white",  # color to remove
  fuzz = 12,           
  output_pdf = "/Users/yacinebenchehida/Desktop/Convergent_evolution/Polymnia/photo/forewing_yellow.pdf"
)

