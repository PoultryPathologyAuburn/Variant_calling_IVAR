###########Number 2 - Finding the total SNPs of P1, P10 and Shared SNPs excluding Indels 
# Set the working directory to where the files are located
setwd("/Users/deepachaudhary/Documents/aminoacidanalysis/bam/iVaroutput/parameterchange/0.03_tsv")

# Define the list of isolate numbers
isolates <- c("iso1", "iso5", "iso6", "iso7", "iso8", "iso9", "iso10")

# Initialize a list to store the counts of intersections and totals for each isolate
intersection_counts <- list()
total_snps_p1 <- list()
total_snps_p10 <- list()

# Loop through each isolate
for (isolate in isolates) {
  # Construct the file paths for P1 and P10 for the current isolate
  file_p1 <- paste0(isolate, "p1C.tsv")
  file_p10 <- paste0(isolate, "p10C.tsv")
  
  # Check if files exist before attempting to read them
  if (file.exists(file_p1) && file.exists(file_p10)) {
    # Load P1 and P10 variant files
    p1 <- read.table(file_p1, header = TRUE, sep = "\t")
    p10 <- read.table(file_p10, header = TRUE, sep = "\t")
    
    # Filter out indels by ensuring the REF and ALT columns have only single nucleotide characters
    p1_snps <- p1[nchar(p1$REF) == 1 & nchar(p1$ALT) == 1, c("POS", "REF", "ALT")]
    p10_snps <- p10[nchar(p10$REF) == 1 & nchar(p10$ALT) == 1, c("POS", "REF", "ALT")]
    
    # Count total SNPs for P1 and P10
    total_snps_p1[[isolate]] <- nrow(p1_snps)
    total_snps_p10[[isolate]] <- nrow(p10_snps)
    
    # Find intersecting SNPs
    intersection <- merge(p1_snps, p10_snps, by = c("POS", "REF", "ALT"))
    
    # Count the number of intersected positions
    num_intersections <- nrow(intersection)
    
    # Save the intersection to a new file only if there are intersecting changes
    if (num_intersections > 0) {
      output_file <- paste0(isolate, "_P1_P10_intersected_SNPs.tsv")
      write.table(intersection, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
      cat("Isolate:", isolate, "- Number of intersected SNPs: ", num_intersections, "\n")
      cat("Intersected SNPs saved to", output_file, "\n")
    } else {
      cat("Isolate:", isolate, "- No intersected SNPs found.\n")
    }
    
    # Store the count of intersections in the list
    intersection_counts[[isolate]] <- num_intersections
  } else {
    cat("Files missing for isolate: ", isolate, " - Check files: ", file_p1, " and ", file_p10, "\n")
  }
}

# Print summary of intersections and total SNPs for all isolates
cat("\nSummary of intersected SNPs for all isolates:\n")
print(intersection_counts)
cat("\nTotal P1 SNPs for all isolates:\n")
print(total_snps_p1)
cat("\nTotal P10 SNPs for all isolates:\n")
print(total_snps_p10)


########### Finding the number of unique SNPs P1 and P10
# Set the working directory to where the files are located
setwd("/Users/deepachaudhary/Documents/aminoacidanalysis/bam/iVaroutput/parameterchange/0.03_tsv")

# Define the list of isolate numbers
isolates <- c("iso1", "iso5", "iso6", "iso7", "iso8", "iso9", "iso10")

# Initialize lists to store the counts of intersections, totals, and unique SNPs for each isolate
intersection_counts <- list()
total_snps_p1 <- list()
total_snps_p10 <- list()
unique_snps_p1 <- list()
unique_snps_p10 <- list()

# Loop through each isolate
for (isolate in isolates) {
  # Construct the file paths for P1 and P10 for the current isolate
  file_p1 <- paste0(isolate, "p1C.tsv")
  file_p10 <- paste0(isolate, "p10C.tsv")
  
  # Check if files exist before attempting to read them
  if (file.exists(file_p1) && file.exists(file_p10)) {
    # Load P1 and P10 variant files
    p1 <- read.table(file_p1, header = TRUE, sep = "\t")
    p10 <- read.table(file_p10, header = TRUE, sep = "\t")
    
    # Filter out indels by ensuring the REF and ALT columns have only single nucleotide characters
    p1_snps <- p1[nchar(p1$REF) == 1 & nchar(p1$ALT) == 1, c("POS", "REF", "ALT")]
    p10_snps <- p10[nchar(p10$REF) == 1 & nchar(p10$ALT) == 1, c("POS", "REF", "ALT")]
    
    # Count total SNPs for P1 and P10
    total_snps_p1[[isolate]] <- nrow(p1_snps)
    total_snps_p10[[isolate]] <- nrow(p10_snps)
    
    # Find intersecting SNPs
    intersection <- merge(p1_snps, p10_snps, by = c("POS", "REF", "ALT"))
    
    # Find unique SNPs in P1
    unique_p1 <- setdiff(p1_snps, p10_snps)
    unique_snps_p1[[isolate]] <- nrow(unique_p1)
    
    # Find unique SNPs in P10
    unique_p10 <- setdiff(p10_snps, p1_snps)
    unique_snps_p10[[isolate]] <- nrow(unique_p10)
    
    # Count the number of intersected positions
    num_intersections <- nrow(intersection)
    
    # Save the intersection and unique SNPs to new files only if there are intersecting/unique changes
    if (num_intersections > 0) {
      output_file_inter <- paste0(isolate, "_P1_P10_intersected_SNPs.tsv")
      write.table(intersection, file = output_file_inter, sep = "\t", row.names = FALSE, quote = FALSE)
      cat("Isolate:", isolate, "- Number of intersected SNPs: ", num_intersections, "\n")
      cat("Intersected SNPs saved to", output_file_inter, "\n")
    }
    if (nrow(unique_p1) > 0) {
      output_file_p1 <- paste0(isolate, "_P1_unique_SNPs.tsv")
      write.table(unique_p1, file = output_file_p1, sep = "\t", row.names = FALSE, quote = FALSE)
      cat("Unique P1 SNPs for", isolate, ": ", nrow(unique_p1), "\n")
      cat("Unique P1 SNPs saved to", output_file_p1, "\n")
    }
    if (nrow(unique_p10) > 0) {
      output_file_p10 <- paste0(isolate, "_P10_unique_SNPs.tsv")
      write.table(unique_p10, file = output_file_p10, sep = "\t", row.names = FALSE, quote = FALSE)
      cat("Unique P10 SNPs for", isolate, ": ", nrow(unique_p10), "\n")
      cat("Unique P10 SNPs saved to", output_file_p10, "\n")
    }
  } else {
    cat("Files missing for isolate: ", isolate, " - Check files: ", file_p1, " and ", file_p10, "\n")
  }
}

# Print summary of intersections and total SNPs for all isolates
cat("\nSummary of intersected SNPs for all isolates:\n")
print(intersection_counts)
cat("\nTotal P1 SNPs for all isolates:\n")
print(total_snps_p1)
cat("\nTotal P10 SNPs for all isolates:\n")
print(total_snps_p10)
cat("\nUnique P1 SNPs for all isolates:\n")
print(unique_snps_p1)
cat("\nUnique P10 SNPs for all isolates:\n")
print(unique_snps_p10)

