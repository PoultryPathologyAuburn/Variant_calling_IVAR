
##################
# Set the working directory to where the files are located
setwd("/Users/deepachaudhary/Documents/aminoacidanalysis/bam/iVaroutput/parameterchange/0.03_tsv")

# Define the list of isolate numbers
isolates <- c("iso1", "iso5", "iso6", "iso7", "iso8", "iso9", "iso10")

# Initialize lists to store counts of intersections, totals, and unique indels for each isolate
intersection_indels <- list()
total_indels_p1 <- list()
total_indels_p10 <- list()
unique_indels_p1 <- list()
unique_indels_p10 <- list()

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
    
    # Filter for indels (where length of REF or ALT > 1)
    p1_indels <- p1[nchar(p1$REF) != 1 | nchar(p1$ALT) != 1, c("POS", "REF", "ALT")]
    p10_indels <- p10[nchar(p10$REF) != 1 | nchar(p10$ALT) != 1, c("POS", "REF", "ALT")]
    
    # Count total indels for P1 and P10
    total_indels_p1[[isolate]] <- nrow(p1_indels)
    total_indels_p10[[isolate]] <- nrow(p10_indels)
    
    # Find intersecting indels
    intersection <- merge(p1_indels, p10_indels, by = c("POS", "REF", "ALT"))
    intersection_indels[[isolate]] <- nrow(intersection)
    
    # Find unique indels in P1 and P10
    unique_p1 <- setdiff(p1_indels, p10_indels)
    unique_indels_p1[[isolate]] <- nrow(unique_p1)
    
    unique_p10 <- setdiff(p10_indels, p1_indels)
    unique_indels_p10[[isolate]] <- nrow(unique_p10)
    
    # Save intersection indels to a file
    write.table(intersection, file = paste0(isolate, "_P1_P10_intersected_indels.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    
    # Save unique P1 indels to a file
    write.table(unique_p1, file = paste0(isolate, "_P1_unique_indels.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    
    # Save unique P10 indels to a file
    write.table(unique_p10, file = paste0(isolate, "_P10_unique_indels.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    
    # Save total counts as summaries in separate text files
    write(paste0("Total P1 indels for ", isolate, ": ", nrow(p1_indels)), file = paste0(isolate, "_total_indels_summary.txt"))
    write(paste0("Total P10 indels for ", isolate, ": ", nrow(p10_indels)), file = paste0(isolate, "_total_indels_summary.txt"), append = TRUE)
    write(paste0("Intersected indels for ", isolate, ": ", nrow(intersection)), file = paste0(isolate, "_total_indels_summary.txt"), append = TRUE)
    write(paste0("Unique P1 indels for ", isolate, ": ", nrow(unique_p1)), file = paste0(isolate, "_total_indels_summary.txt"), append = TRUE)
    write(paste0("Unique P10 indels for ", isolate, ": ", nrow(unique_p10)), file = paste0(isolate, "_total_indels_summary.txt"), append = TRUE)
  } else {
    cat("Files missing for isolate: ", isolate, " - Check files: ", file_p1, " and ", file_p10, "\n")
  }
}

# Print summary of indels for all isolates
cat("\nSummary of intersected indels for all isolates:\n")
print(intersection_indels)
cat("\nTotal P1 indels for all isolates:\n")
print(total_indels_p1)
cat("\nTotal P10 indels for all isolates:\n")
print(total_indels_p10)
cat("\nUnique P1 indels for all isolates:\n")
print(unique_indels_p1)
cat("\nUnique P10 indels for all isolates:\n")
print(unique_indels_p10)

################# Annotation
# Set the working directory to where the files are located
setwd("/Users/deepachaudhary/Documents/aminoacidanalysis/bam/iVaroutput/parameterchange/0.03_tsv")

# Define genomic regions (including HN and other regions)
genomic_regions <- data.frame(
  Region = c("NP", "P", "M", "F", "HN", "L"),
  Start = c(56, 1804, 3256, 4498, 6321, 8370),
  End = c(1792, 3244, 4487, 6279, 8312, 15073)
)

# Define the list of isolate numbers
isolates <- c("iso1", "iso5", "iso6", "iso7", "iso8", "iso9", "iso10")

# Function to annotate indels based on genomic regions
annotate_by_position <- function(indel_data) {
  indel_data$Region <- NA  # Create a column for region annotation
  
  # Loop through each genomic region
  for (i in 1:nrow(genomic_regions)) {
    region_name <- genomic_regions$Region[i]
    region_start <- genomic_regions$Start[i]
    region_end <- genomic_regions$End[i]
    
    # Assign the region name to indels falling within the region's start and end
    indel_data$Region[indel_data$POS >= region_start & indel_data$POS <= region_end] <- region_name
  }
  
  return(indel_data)
}

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
    
    # Filter for indels (where length of REF or ALT > 1)
    p1_indels <- p1[nchar(p1$REF) != 1 | nchar(p1$ALT) != 1, ]
    p10_indels <- p10[nchar(p10$REF) != 1 | nchar(p10$ALT) != 1, ]
    
    # Annotate indels by position using the genomic regions
    p1_annotated <- annotate_by_position(p1_indels)
    p10_annotated <- annotate_by_position(p10_indels)
    
    # Save annotated P1 indels to a file
    write.table(p1_annotated, file = paste0(isolate, "_P1_annotated_indels.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    
    # Save annotated P10 indels to a file
    write.table(p10_annotated, file = paste0(isolate, "_P10_annotated_indels.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    
  } else {
    cat("Files missing for isolate: ", isolate, " - Check files: ", file_p1, " and ", file_p10, "\n")
  }
}

cat("\nIndel extraction and annotation completed for P1 and P10 for all isolates.\n")

#######classification
# Set the working directory to where the files are located
setwd("/Users/deepachaudhary/Documents/aminoacidanalysis/bam/iVaroutput/parameterchange/0.03_tsv")

# Define genomic regions (including HN and other regions)
genomic_regions <- data.frame(
  Region = c("NP", "P", "M", "F", "HN", "L"),
  Start = c(56, 1804, 3256, 4498, 6321, 8370),
  End = c(1792, 3244, 4487, 6279, 8312, 15073)
)

# Define the list of isolate numbers
isolates <- c("iso1", "iso5", "iso6", "iso7", "iso8", "iso9", "iso10")

# Function to annotate indels based on genomic regions and classify as in-frame or frameshift
annotate_and_classify_indels <- function(indel_data) {
  indel_data$Region <- NA  # Create a column for region annotation
  indel_data$Classification <- NA  # Create a column for in-frame or frameshift classification
  
  # Annotate the indel based on genomic regions
  for (i in 1:nrow(genomic_regions)) {
    region_name <- genomic_regions$Region[i]
    region_start <- genomic_regions$Start[i]
    region_end <- genomic_regions$End[i]
    
    indel_data$Region[indel_data$POS >= region_start & indel_data$POS <= region_end] <- region_name
  }
  
  # Classify indels as in-frame or frameshift
  for (i in 1:nrow(indel_data)) {
    ref_length <- nchar(indel_data$REF[i])
    alt_length <- nchar(indel_data$ALT[i])
    indel_length_diff <- abs(ref_length - alt_length)
    
    if (indel_length_diff %% 3 == 0) {
      indel_data$Classification[i] <- "In-frame"
    } else {
      indel_data$Classification[i] <- "Frameshift"
    }
  }
  
  return(indel_data)
}

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
    
    # Filter for indels (where length of REF or ALT > 1)
    p1_indels <- p1[nchar(p1$REF) != 1 | nchar(p1$ALT) != 1, ]
    p10_indels <- p10[nchar(p10$REF) != 1 | nchar(p10$ALT) != 1, ]
    
    # Annotate and classify indels for P1 and P10
    p1_annotated <- annotate_and_classify_indels(p1_indels)
    p10_annotated <- annotate_and_classify_indels(p10_indels)
    
    # Save annotated and classified P1 indels to a file
    write.table(p1_annotated, file = paste0(isolate, "_P1_annotated_indels.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    
    # Save annotated and classified P10 indels to a file
    write.table(p10_annotated, file = paste0(isolate, "_P10_annotated_indels.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    
  } else {
    cat("Files missing for isolate: ", isolate, " - Check files: ", file_p1, " and ", file_p10, "\n")
  }
}

cat("\nIndel extraction, annotation, and classification completed for P1 and P10 for all isolates.\n")

