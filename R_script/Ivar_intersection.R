# Define the list of isolate numbers
isolates <- c("iso1", "iso5", "iso6", "iso7", "iso8", "iso9", "iso10")

# Initialize a list to store the counts of intersections for each isolate
intersection_counts <- list()

# Loop through each isolate
for (isolate in isolates) {
  # Construct the file paths for P1 and P10 for the current isolate
  file_p1 <- paste0(isolate, "p1.tsv")
  file_p10 <- paste0(isolate, "p10.tsv")
  
  # Load P1 and P10 variant files
  p1 <- read.table(file_p1, header = TRUE, sep = "\t")
  p10 <- read.table(file_p10, header = TRUE, sep = "\t")
  
  # Select relevant columns: POS, REF_AA, ALT_AA
  p1_aa <- p1[, c("POS", "REF_AA", "ALT_AA")]
  p10_aa <- p10[, c("POS", "REF_AA", "ALT_AA")]
  
  # Remove rows with NA values
  p1_aa <- na.omit(p1_aa)
  p10_aa <- na.omit(p10_aa)
  
  # Remove rows where REF_AA is the same as ALT_AA
  p1_aa <- p1_aa[p1_aa$REF_AA != p1_aa$ALT_AA, ]
  p10_aa <- p10_aa[p10_aa$REF_AA != p10_aa$ALT_AA, ]
  
  # Find intersecting amino acid changes
  intersection <- merge(p1_aa, p10_aa, by = c("POS", "REF_AA", "ALT_AA"))
  
  # Count the number of intersected positions
  num_intersections <- nrow(intersection)
  
  # Save the intersection to a new file only if there are intersecting changes
  if (num_intersections > 0) {
    output_file <- paste0(isolate, "_P1_P10_intersected_AA_changes.tsv")
    write.table(intersection, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
    cat("Isolate:", isolate, "- Number of intersected amino acid changes: ", num_intersections, "\n")
    cat("Intersected amino acid changes saved to", output_file, "\n")
  } else {
    cat("Isolate:", isolate, "- No intersected amino acid changes found.\n")
  }
  
  # Store the count of intersections in the list
  intersection_counts[[isolate]] <- num_intersections
}

# Print summary of intersections for all isolates
cat("\nSummary of intersected amino acid changes for all isolates:\n")
print(intersection_counts)
