 #######Extraction of nSNPs directly from IVar output and proceed the change in frequency

library(readr)
library(dplyr)
library(ggplot2)

# Set the working directory where the original IVAR output files are located
setwd("/Users/deepachaudhary/Documents/aminoacidanalysis/bam/iVaroutput/parameterchange/0.03_tsv")

# List of isolate identifiers
isolates <- c("Iso1", "Iso5", "Iso6", "Iso7", "Iso8", "Iso9", "Iso10")

# Function to analyze non-synonymous SNP frequency changes for one isolate
analyze_isolate <- function(isolate) {
  # Construct file paths for P1 and P10 TSV files
  p1_file <- paste0(isolate, "p1C.tsv")
  p10_file <- paste0(isolate, "p10C.tsv")
  
  # Ensure files exist before proceeding
  if (file.exists(p1_file) && file.exists(p10_file)) {
    # Load SNP data for P1 and P10 and filter for non-synonymous SNPs
    data_p1 <- read_tsv(p1_file) %>%
      filter(REF_AA != ALT_AA) %>%
      mutate(Isolates = isolate,  # Add isolate column manually
             POS = as.numeric(POS),
             SNP_Frequency_P1 = ifelse(TOTAL_DP > 0, ALT_DP / TOTAL_DP, 0))
    
    data_p10 <- read_tsv(p10_file) %>%
      filter(REF_AA != ALT_AA) %>%
      mutate(Isolates = isolate,  # Add isolate column manually
             POS = as.numeric(POS),
             SNP_Frequency_P10 = ifelse(TOTAL_DP > 0, ALT_DP / TOTAL_DP, 0))
    
    # Merge P1 and P10 data
    comparison <- data_p1 %>%
      select(Isolates, POS, REF, ALT, REF_AA, ALT_AA, SNP_Frequency_P1) %>%
      left_join(data_p10 %>%
                  select(Isolates, POS, REF, ALT, REF_AA, ALT_AA, SNP_Frequency_P10),
                by = c("Isolates", "POS", "REF", "ALT", "REF_AA", "ALT_AA"))
    
    # Calculate frequency changes
    comparison <- comparison %>%
      mutate(Frequency_Change = SNP_Frequency_P10 - SNP_Frequency_P1)
    
    # Filter for significant changes
    significant_changes <- comparison %>%
      filter(!is.na(Frequency_Change), abs(Frequency_Change) > 0.2)
    
    # Save the results
    write_csv(comparison, paste0(isolate, "_All_Non_Synonymous_Frequency_Changes.csv"))
    write_csv(significant_changes, paste0(isolate, "_Significant_Non_Synonymous_Frequency_Changes.csv"))
    
    # Plotting all frequency changes
    plot <- ggplot(comparison, aes(x = POS, y = Frequency_Change, color = Frequency_Change > 0)) +
      geom_point() +
      theme_minimal() +
      labs(title = paste("Non-Synonymous SNP Frequency Changes from P1 to P10 for", isolate),
           x = "Genomic Position", y = "Frequency Change") +
      scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "red"))
    
    # Save the plot
    ggsave(paste0(isolate, "_Frequency_Change_Plot.png"), plot = plot)
    
    print(paste("Non-synonymous SNP analysis complete for", isolate))
  } else {
    print(paste("File missing for", isolate, "- skipping analysis"))
  }
}

# Loop through all isolates and process them
for (isolate in isolates) {
  analyze_isolate(isolate)
}

############
#To combine all the results _All_Non_Synonymous_Frequency_Changes.csv into one file

library(readr)
library(dplyr)

# Set the working directory
setwd("/Users/deepachaudhary/Documents/aminoacidanalysis/bam/iVaroutput/parameterchange/0.03_tsv")

# Define genomic regions
genomic_regions <- data.frame(
  Region = c("NP", "P", "M", "F", "HN", "L"),
  Start = c(56, 1804, 3256, 4498, 6321, 8370),
  End = c(1792, 3244, 4487, 6279, 8312, 15073)
)

# Function to determine genomic region based on position
annotate_region <- function(pos) {
  for (i in 1:nrow(genomic_regions)) {
    if (pos >= genomic_regions$Start[i] && pos <= genomic_regions$End[i]) {
      return(genomic_regions$Region[i])
    }
  }
  return("Intergenic") # Return 'Intergenic' if no region matches
}

# Load the combined non-synonymous SNP frequency changes data
data <- read_csv("Combined_All_Non_Synonymous_Frequency_Changes.csv")

# Annotate genomic region
data <- data %>%
  mutate(Genomic_Region = sapply(POS, annotate_region))

# Save the updated data with genomic region annotations
write_csv(data, "Combined_All_Non_Synonymous_Frequency_Changes_With_Region.csv")

print("Data with genomic regions saved successfully.")

#########
#To combine all the results  _Significant_Non_Synonymous_Frequency_Changes.csv
library(readr)
library(dplyr)

# Set the working directory
setwd("/Users/deepachaudhary/Documents/aminoacidanalysis/bam/iVaroutput/parameterchange/0.03_tsv")

# Define genomic regions
genomic_regions <- data.frame(
  Region = c("NP", "P", "M", "F", "HN", "L"),
  Start = c(56, 1804, 3256, 4498, 6321, 8370),
  End = c(1792, 3244, 4487, 6279, 8312, 15073)
)

# Function to determine genomic region based on position
annotate_region <- function(pos) {
  for (i in 1:nrow(genomic_regions)) {
    if (pos >= genomic_regions$Start[i] && pos <= genomic_regions$End[i]) {
      return(genomic_regions$Region[i])
    }
  }
  return("Intergenic") # Return 'Intergenic' if no region matches
}

# List of isolate identifiers
isolates <- c("Iso1", "Iso5", "Iso6", "Iso7", "Iso8", "Iso9", "Iso10")

# Initialize an empty data frame to combine all significant changes
all_significant_changes <- data.frame()

# Loop through each isolate to read, annotate, and combine significant changes
for(isolate in isolates) {
  file_name <- paste(isolate, "_Significant_Non_Synonymous_Frequency_Changes.csv", sep="")
  if (file.exists(file_name)) {
    # Read the significant changes file
    changes <- read_csv(file_name)
    # Annotate genomic region
    changes <- changes %>%
      mutate(Genomic_Region = sapply(POS, annotate_region))
    # Combine with the main data frame
    all_significant_changes <- bind_rows(all_significant_changes, changes)
  } else {
    message(paste("File not found for", isolate, "- skipping."))
  }
}

# Save the combined and annotated data to a new CSV file
write_csv(all_significant_changes, "Combined_Significant_Non_Synonymous_Frequency_Changes_With_Region.csv")

print("All significant changes combined, annotated, and saved.")

##########scatterplot plot
#To visualize the significant non-synonymous SNP frequency changes in scatterplot plot
library(ggplot2)
library(dplyr)
library(readr)

# Set the working directory
setwd("/Users/deepachaudhary/Documents/aminoacidanalysis/bam/iVaroutput/parameterchange/0.03_tsv")

# Load data
data <- read_csv("Combined_Significant_Non_Synonymous_Frequency_Changes_With_Region.csv")

# Convert isolate names
data$Isolates <- recode(data$Isolates,
                        "Iso1" = "Iso1",
                        "Iso5" = "Iso2",
                        "Iso6" = "Iso3",
                        "Iso7" = "Iso4",
                        "Iso8" = "Iso5",
                        "Iso9" = "Iso6",
                        "Iso10" = "Iso7")

# Define colors for each isolate
colors <- c("Iso1" = "blue", "Iso2" = "red", "Iso3" = "green", 
            "Iso4" = "orange", "Iso5" = "purple", "Iso6" = "brown", "Iso7" = "pink")

# Define genomic regions
genomic_regions <- data.frame(
  Region = c("NP", "P", "M", "F", "HN", "L"),
  Start = c(56, 1804, 3256, 4498, 6321, 8370),
  End = c(1792, 3244, 4487, 6279, 8312, 15073),
  Color = c("lightblue", "lightgreen", "lightpink", "lightyellow", "lightcyan", "lightgrey")
)

# Create the scatter plot with genomic regions
snp_plot <- ggplot() +
  # Add genomic region rectangles
  geom_rect(data = genomic_regions, 
            aes(xmin = Start, xmax = End, ymin = -Inf, ymax = Inf, fill = Region), 
            alpha = 0.3) +
  # Plot SNP points
  geom_point(data = data, aes(x = POS, y = Isolates, color = Isolates), size = 3) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = setNames(genomic_regions$Color, genomic_regions$Region)) +
  labs(x = "Genomic Position", y = "Isolate", color = "Isolate", fill = "Genomic Region") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "right") +
  scale_x_continuous(breaks = genomic_regions$Start, labels = genomic_regions$Region)

# Save the plot
ggsave("Significant_nSNPs_Scatter_Plot_with_Regions.jpeg", plot = snp_plot, width = 10, height = 8, dpi = 300)

# Display the plot 
print(snp_plot)
