##To find the total number of nSNPs of P1 and P10, shared and unique AA

##########
# Set the working directory
setwd("/Users/deepachaudhary/Documents/aminoacidanalysis/bam/iVaroutput/parameterchange/0.03_tsv")

# Define the list of isolate numbers
isolates <- c("iso1", "iso5", "iso6", "iso7", "iso8", "iso9", "iso10")

# Loop through each isolate
for (isolate in isolates) {
  # Construct file paths for P1 and P10 for the current isolate
  file_p1 <- paste0(isolate, "p1C.tsv")
  file_p10 <- paste0(isolate, "p10C.tsv")
  
  # Check if files exist before attempting to read them
  if (file.exists(file_p1) && file.exists(file_p10)) {
    # Load P1 and P10 variant files
    p1 <- read.table(file_p1, header = TRUE, sep = "\t")
    p10 <- read.table(file_p10, header = TRUE, sep = "\t")
    
    # Remove rows where any of the key fields are NA
    p1 <- na.omit(p1[, c("POS", "REF_AA", "ALT_AA")])
    p10 <- na.omit(p10[, c("POS", "REF_AA", "ALT_AA")])
    
    # Filter for non-synonymous SNPs by checking amino acid changes and ensuring single nucleotide changes
    p1_nsnp <- p1[p1$REF_AA != p1$ALT_AA & nchar(p1$REF) == 1 & nchar(p1$ALT) == 1, ]
    p10_nsnp <- p10[p10$REF_AA != p10$ALT_AA & nchar(p10$REF) == 1 & nchar(p10$ALT) == 1, ]
    
    # Identify unique rows by making unique identifiers
    p1_nsnp$unique_id <- paste(p1_nsnp$POS, p1_nsnp$REF_AA, p1_nsnp$ALT_AA, sep="_")
    p10_nsnp$unique_id <- paste(p10_nsnp$POS, p10_nsnp$REF_AA, p10_nsnp$ALT_AA, sep="_")
    
    # Find shared and unique non-synonymous SNPs
    shared_nsnp <- p1_nsnp[p1_nsnp$unique_id %in% p10_nsnp$unique_id, ]
    unique_p1_nsnp <- p1_nsnp[!p1_nsnp$unique_id %in% p10_nsnp$unique_id, ]
    unique_p10_nsnp <- p10_nsnp[!p10_nsnp$unique_id %in% p1_nsnp$unique_id, ]
    
    # Save the filtered non-synonymous SNP data and shared/unique data to files
    write.table(p1_nsnp, paste0(isolate, "_P1_non_synonymous_SNPs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(p10_nsnp, paste0(isolate, "_P10_non_synonymous_SNPs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(shared_nsnp, paste0(isolate, "_Shared_non_synonymous_SNPs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(unique_p1_nsnp, paste0(isolate, "_Unique_P1_non_synonymous_SNPs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(unique_p10_nsnp, paste0(isolate, "_Unique_P10_non_synonymous_SNPs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Print summary to console
    cat(isolate, "- Total non-synonymous SNPs in P1: ", nrow(p1_nsnp), ", in P10: ", nrow(p10_nsnp), "\n")
    cat(isolate, "- Shared non-synonymous SNPs: ", nrow(shared_nsnp), "\n")
    cat(isolate, "- Unique non-synonymous SNPs in P1: ", nrow(unique_p1_nsnp), ", in P10: ", nrow(unique_p10_nsnp), "\n")
    
    # Clean up added columns
    p1_nsnp$unique_id <- NULL
    p10_nsnp$unique_id <- NULL
  } else {
    cat("Files missing for isolate: ", isolate, " - Check files: ", file_p1, " and ", file_p10, "\n")
  }
}

###########Adding annotation 
# Set the working directory
setwd("/Users/deepachaudhary/Documents/aminoacidanalysis/bam/iVaroutput/parameterchange/0.03_tsv")

# Define the list of isolate numbers
isolates <- c("iso1", "iso5", "iso6", "iso7", "iso8", "iso9", "iso10")

# Define genomic regions
genomic_regions <- data.frame(
  Region = c("NP", "P", "M", "F", "HN", "L"),
  Start = c(56, 1804, 3256, 4498, 6321, 8370),
  End = c(1792, 3244, 4487, 6279, 8312, 15073)
)

# Function to annotate genomic region based on position
annotate_genomic_region <- function(position) {
  for (i in 1:nrow(genomic_regions)) {
    if (position >= genomic_regions$Start[i] && position <= genomic_regions$End[i]) {
      return(genomic_regions$Region[i])
    }
  }
  return(NA)  # Return NA if no region is found
}

# Loop through each isolate
for (isolate in isolates) {
  # Construct file paths for P1 and P10 for the current isolate
  file_p1 <- paste0(isolate, "p1C.tsv")
  file_p10 <- paste0(isolate, "p10C.tsv")
  
  # Check if files exist before attempting to read them
  if (file.exists(file_p1) && file.exists(file_p10)) {
    # Load P1 and P10 variant files
    p1 <- read.table(file_p1, header = TRUE, sep = "\t")
    p10 <- read.table(file_p10, header = TRUE, sep = "\t")
    
    # Remove rows where any of the key fields are NA
    p1 <- na.omit(p1[, c("POS", "REF_AA", "ALT_AA")])
    p10 <- na.omit(p10[, c("POS", "REF_AA", "ALT_AA")])
    
    # Filter for non-synonymous SNPs by checking amino acid changes and ensuring single nucleotide changes
    p1_nsnp <- p1[p1$REF_AA != p1$ALT_AA & nchar(p1$REF) == 1 & nchar(p1$ALT) == 1, ]
    p10_nsnp <- p10[p10$REF_AA != p10$ALT_AA & nchar(p10$REF) == 1 & nchar(p10$ALT) == 1, ]
    
    # Annotate genomic regions
    p1_nsnp$Region <- sapply(p1_nsnp$POS, annotate_genomic_region)
    p10_nsnp$Region <- sapply(p10_nsnp$POS, annotate_genomic_region)
    
    # Identify unique rows by making unique identifiers
    p1_nsnp$unique_id <- paste(p1_nsnp$POS, p1_nsnp$REF_AA, p1_nsnp$ALT_AA, sep="_")
    p10_nsnp$unique_id <- paste(p10_nsnp$POS, p10_nsnp$REF_AA, p10_nsnp$ALT_AA, sep="_")
    
    # Find shared and unique non-synonymous SNPs
    shared_nsnp <- p1_nsnp[p1_nsnp$unique_id %in% p10_nsnp$unique_id, ]
    unique_p1_nsnp <- p1_nsnp[!p1_nsnp$unique_id %in% p10_nsnp$unique_id, ]
    unique_p10_nsnp <- p10_nsnp[!p10_nsnp$unique_id %in% p1_nsnp$unique_id, ]
    
    # Save the filtered non-synonymous SNP data and shared/unique data to files
    write.table(p1_nsnp, paste0(isolate, "_P1_non_synonymous_SNPs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(p10_nsnp, paste0(isolate, "_P10_non_synonymous_SNPs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(shared_nsnp, paste0(isolate, "_Shared_non_synonymous_SNPs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(unique_p1_nsnp, paste0(isolate, "_Unique_P1_non_synonymous_SNPs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(unique_p10_nsnp, paste0(isolate, "_Unique_P10_non_synonymous_SNPs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Print summary to console
    cat(isolate, "- Total non-synonymous SNPs in P1: ", nrow(p1_nsnp), ", in P10: ", nrow(p10_nsnp), "\n")
    cat(isolate, "- Shared non-synonymous SNPs: ", nrow(shared_nsnp), "\n")
    cat(isolate, "- Unique non-synonymous SNPs in P1: ", nrow(unique_p1_nsnp), ", in P10: ", nrow(unique_p10_nsnp), "\n")
    
    # Clean up added columns
    p1_nsnp$unique_id <- NULL
    p10_nsnp$unique_id <- NULL
  } else {
    cat("Files missing for isolate: ", isolate, " - Check files: ", file_p1, " and ", file_p10, "\n")
  }
}

