file_list <- list.files(pattern = "*readcount")

# Loop through each file
for (file in file_list) {
  # Check if file is empty or non-existent before reading
  if (!file.exists(file) || file.size(file) == 0) {
    cat("Skipping empty or missing file:", file, "\n")
    next
  }
  
  # Read the ribo data with error handling
  df <- tryCatch(
    read.delim(file, header = FALSE),
    error = function(e) {
      cat("Error reading file:", file, "-", e$message, "\n")
      return(NULL)
    }
  )
  
  if (is.null(df) || nrow(df) == 0) {
    cat("No data in file:", file, "\n")
    next
  }
  
  # Extract the sample name from the file name
  sample_name <- sub(".readcount", "", file)

  A <- sapply(as.character(df$V6), function(x) {
    split_string <- strsplit(x, ":")[[1]]
    if (length(split_string) >= 2) {
      split_string[2]
    } else {
      NA
    }
  })
  df$A <- gsub('"', '', A)
  df$A <- as.numeric(df$A)
  
  C <- sapply(as.character(df$V7), function(x) {
    split_string <- strsplit(x, ":")[[1]]
    if (length(split_string) >= 2) {
      split_string[2]
    } else {
      NA
    }
  })
  df$C <- gsub('"', '', C)
  df$C <- as.numeric(df$C)
  
  G <- sapply(as.character(df$V8), function(x) {
    split_string <- strsplit(x, ":")[[1]]
    if (length(split_string) >= 2) {
      split_string[2]
    } else {
      NA
    }
  })
  df$G <- gsub('"', '', G)
  df$G <- as.numeric(df$G)
  
  T <- sapply(as.character(df$V9), function(x) {
    split_string <- strsplit(x, ":")[[1]]
    if (length(split_string) >= 2) {
      split_string[2]
    } else {
      NA
    }
  })
  df$T <- gsub('"', '', T)
  df$T <- as.numeric(df$T)

  # Assuming your data frame is named df
  df <- subset(df, select = -c(V5, V6, V7, V8, V9, V10))

  colnames(df)[colnames(df) == "V1"] <- "chr"
  colnames(df)[colnames(df) == "V2"] <- "pos"
  colnames(df)[colnames(df) == "V3"] <- "ref_base"
  colnames(df)[colnames(df) == "V4"] <- "Total_count"
 
  output_file_name <- paste0(sample_name, ".txt")
  write.table(df, file = output_file_name, sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat("Processed:", file, "->", output_file_name, "\n")
}
