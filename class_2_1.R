
# function to classify DEGs
classify_gene <- function(logFC, padj) {
  if (is.na(padj)) padj <- 1  # handle missing p-value
  if (logFC > 1 & padj < 0.05) {
    return("Upregulated")
  } else if (logFC < -1 & padj < 0.05) {
    return("Downregulated")
  } else {
    return("Not_Significant")
  }
}

getwd()
print("hello world")

# input/output directories
input_dir <- "raw_data"
output_dir <- "results"

# create output folder if not exist
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}

# list files
files_to_process <- c("DEGs_Data_1.csv", "DEGs_Data_2.csv")

# prepare empty list to store results
deg_results_list <- list()

# loop through each file
for (file_name in files_to_process) {
  cat("\nProcessing:", file_name, "\n")
  
  # build full file path
  input_file_path <- file.path(input_dir, file_name)
  
  # load csv
  data_2r <- read.csv(input_file_path, header = TRUE)
  cat("File imported. Checking for missing values...\n")
  
  
  # replace missing padj values with 1
  missing_count <- sum(is.na( data_2r$padj))
  cat("Missing padj values:", missing_count, "\n")
  data_2r$padj[is.na( data_2r$padj)] <- 1  # safe operation
  
  # classify gene status
  data_2r$status <- mapply(classify_gene,  data_2r$logFC,  data_2r$padj)
  
  # save to list
  deg_results_list[[file_name]] <-  data_2r
  
  # save output to file
  output_file_path <- file.path(output_dir, paste0("DEG_results_", file_name))
  write.csv( data_2r, output_file_path, row.names = FALSE)
  
  cat("Results saved to:", output_file_path, "\n")
  cat("Gene classification summary:\n")
  print(table( data_2r$status))
}

# save the entire R workspace

save.image(file = "AmnaMAnis_Class_2_Assignment.RData")
