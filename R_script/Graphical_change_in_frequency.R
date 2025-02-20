# Install and load necessary packages
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
library(tidyverse)
library(gridExtra)

# Verify data loading
file_path <- "Combined_Significant_Non_Synonymous_Frequency_Changes_With_Region.csv"
data <- read.csv(file_path)
print("Data Structure:")
str(data)

# Check for required columns
required_columns <- c("Genomic_Region", "Frequency_Change", "SNP_Frequency_P1", "SNP_Frequency_P10")
missing_columns <- setdiff(required_columns, colnames(data))
if (length(missing_columns) > 0) {
  stop(paste("Missing columns:", paste(missing_columns, collapse = ", ")))
}

# Debugging regional statistics
regional_stats <- calculate_regional_stats(data)
print("Regional Statistics:")
print(regional_stats)

# Debugging plots
tryCatch({
  freq_plot <- plot_frequency_changes(data)
  print(freq_plot)
}, error = function(e) {
  cat("Error in frequency plot:", e$message, "\n")
})

tryCatch({
  property_plot <- plot_property_changes(data)
  print(property_plot)
}, error = function(e) {
  cat("Error in property plot:", e$message, "\n")
})



# Required libraries
library(tidyverse)
library(ggplot2)
library(scales)
library(gridExtra)

# Set the working directory
setwd("~/Documents/aminoacidanalysis/bam/iVaroutput/parameterchange/verification/changeinfrequency")

# Function to read and prepare data
prepare_data <- function(file_path) {
  data <- read.csv(file_path)
  # Convert columns to appropriate types
  data$Isolates <- as.factor(data$Isolates)
  data$Genomic_Region <- as.factor(data$Genomic_Region)
  return(data)
}

# Function to analyze amino acid properties
get_aa_properties <- function(aa) {
  properties <- list(
    'G' = 'Small/Nonpolar',
    'A' = 'Small/Nonpolar',
    'V' = 'Nonpolar',
    'L' = 'Nonpolar',
    'I' = 'Nonpolar',
    'F' = 'Aromatic/Nonpolar',
    'W' = 'Aromatic/Nonpolar',
    'Y' = 'Aromatic/Polar',
    'S' = 'Polar',
    'T' = 'Polar',
    'N' = 'Polar',
    'Q' = 'Polar',
    'C' = 'Polar',
    'M' = 'Nonpolar',
    'P' = 'Nonpolar',
    'H' = 'Basic/Polar',
    'K' = 'Basic/Polar',
    'R' = 'Basic/Polar',
    'D' = 'Acidic/Polar',
    'E' = 'Acidic/Polar'
  )
  return(properties[[aa]])
}

# Function to analyze amino acid changes
analyze_aa_changes <- function(data) {
  data %>%
    mutate(
      REF_Property = sapply(REF_AA, get_aa_properties),
      ALT_Property = sapply(ALT_AA, get_aa_properties),
      Mutation = paste0(REF_AA, POS, ALT_AA)
    ) %>%
    select(Mutation, Genomic_Region, REF_Property, ALT_Property, Frequency_Change)
}

# Function for statistical analysis by region
calculate_regional_stats <- function(data) {
  data %>%
    group_by(Genomic_Region) %>%
    summarise(
      Count = n(),
      Mean_Change = mean(Frequency_Change),
      SD_Change = sd(Frequency_Change),
      Max_Increase = max(Frequency_Change),
      Max_Decrease = min(Frequency_Change),
      .groups = 'drop'
    )
}

# Function to create frequency change visualization
plot_frequency_changes <- function(data) {
  data %>%
    mutate(Mutation = paste0(REF_AA, POS, ALT_AA)) %>%
    ggplot(aes(x = reorder(Mutation, Frequency_Change), y = Frequency_Change, fill = Genomic_Region)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "SNP Frequency Changes Between P1 and P10",
      x = "Mutation",
      y = "Frequency Change"
    )
}

# Function to create property change heatmap
plot_property_changes <- function(data) {
  aa_changes <- analyze_aa_changes(data)
  
  ggplot(aa_changes, aes(x = REF_Property, y = ALT_Property, fill = Frequency_Change)) +
    geom_tile() +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Amino Acid Property Changes vs Frequency Change",
      x = "Original Property",
      y = "New Property"
    )
}

# Function to perform temporal analysis
analyze_temporal_patterns <- function(data) {
  # Compare P1 and P10 frequencies
  data %>%
    group_by(Genomic_Region) %>%
    summarise(
      Mean_P1 = mean(SNP_Frequency_P1),
      Mean_P10 = mean(SNP_Frequency_P10),
      Delta = Mean_P10 - Mean_P1,
      .groups = 'drop'
    )
}

# Main analysis function
main_analysis <- function(file_path) {
  # Read and prepare data
  data <- prepare_data(file_path)
  
  # Perform analyses
  regional_stats <- calculate_regional_stats(data)
  aa_changes <- analyze_aa_changes(data)
  temporal_patterns <- analyze_temporal_patterns(data)
  
  # Create visualizations
  freq_plot <- plot_frequency_changes(data)
  property_plot <- plot_property_changes(data)
  
  # Statistical tests
  t_test_result <- t.test(data$SNP_Frequency_P1, data$SNP_Frequency_P10, paired = TRUE)
  anova_result <- aov(Frequency_Change ~ Genomic_Region, data = data)
  
  # Create results list
  results <- list(
    regional_stats = regional_stats,
    aa_changes = aa_changes,
    temporal_patterns = temporal_patterns,
    t_test = t_test_result,
    anova = summary(anova_result),
    plots = list(
      frequency_plot = freq_plot,
      property_plot = property_plot
    )
  )
  
  return(results)
}

# Run the analysis
results <- main_analysis("Combined_Significant_Non_Synonymous_Frequency_Changes_With_Region.csv")

# Function to generate comprehensive report
generate_report <- function(results) {
  # Print regional statistics
  cat("Regional Statistics:\n")
  print(results$regional_stats)
  cat("\n")
  
  # Print temporal patterns
  cat("Temporal Patterns:\n")
  print(results$temporal_patterns)
  cat("\n")
  
  # Print statistical test results
  cat("Statistical Tests:\n")
  cat("Paired t-test (P1 vs P10):\n")
  print(results$t_test)
  cat("\nANOVA (Regional Differences):\n")
  print(results$anova)
  
  # Display plots
  grid.arrange(
    results$plots$frequency_plot,
    results$plots$property_plot,
    ncol = 1
  )
}

# Generate the report
generate_report(results)
