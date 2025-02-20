#Scatterplot showing each isolate has its unique nSNPs for P1 shown just above its unique nSNPs for P10
#Likewise Shared nSNPs between P1 and P10 plotted just below P10 SNPs in the sameplot

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)

# Set the working directory
setwd("~/Documents/aminoacidanalysis/bam/iVaroutput/parameterchange/verification")

# Load SNP data from Excel files
p1_snp_data <- read_excel("Uni_NSNPs_P1.xlsx")
p10_snp_data <- read_excel("Uni_NSNPs_P10.xlsx")
shared_snp_data <- read_excel("Shared_NSNPs.xlsx")  # Load the shared SNPs

# Convert the data to long format and adjust positions for plotting
p1_snp_long <- p1_snp_data %>%
  pivot_longer(cols = everything(), names_to = "Isolate", values_to = "POS") %>%
  mutate(Isolate = tools::toTitleCase(Isolate),
         SNP_Type = "Unique P1 nSNPs", Y_Position = as.numeric(factor(Isolate)) + 0.5)  # Position for P1 slightly above

p10_snp_long <- p10_snp_data %>%
  pivot_longer(cols = everything(), names_to = "Isolate", values_to = "POS") %>%
  mutate(Isolate = tools::toTitleCase(Isolate),
         SNP_Type = "Unique P10 nSNPs", Y_Position = as.numeric(factor(Isolate)) + 0.2)  # Position for P10 in the middle

shared_snp_long <- shared_snp_data %>%
  pivot_longer(cols = everything(), names_to = "Isolate", values_to = "POS") %>%
  mutate(Isolate = tools::toTitleCase(Isolate),
         SNP_Type = "Shared nSNPs between P1P10", Y_Position = as.numeric(factor(Isolate)))  # Position for Shared SNPs as the baseline

# Combine all three datasets (P1, P10, Shared)
combined_snp_data <- bind_rows(p1_snp_long, p10_snp_long, shared_snp_long)

# Define updated shapes for SNP types
shapes <- c("Unique P1 nSNPs" = 15, "Unique P10 nSNPs" = 17, "Shared nSNPs between P1P10" = 16)  # 15: square, 17: triangle, 16: circle

# Define colors for each isolate
colors <- c("Iso1" = "blue", "Iso2" = "red", "Iso3" = "green", "Iso4" = "orange", 
            "Iso5" = "purple", "Iso6" = "brown", "Iso7" = "pink")

# Define the genomic regions and their colors for the legend
genomic_regions <- data.frame(
  Region = c("NP", "P", "M", "F", "HN", "L"),
  Start = c(56, 1804, 3256, 4498, 6321, 8370),
  End = c(1792, 3244, 4487, 6279, 8312, 15073)
)
genomic_colors <- c("NP" = "lightblue", "P" = "lightgreen", "M" = "lightpink", 
                    "F" = "lightyellow", "HN" = "lightcyan", "L" = "lightgrey")

# Create the scatterplot
snp_plot <- ggplot(combined_snp_data, aes(x = POS, y = Y_Position, color = Isolate, shape = SNP_Type)) +
  geom_point(size = 4) +
  scale_shape_manual(values = shapes) +  # Updated custom shapes for SNP types
  scale_color_manual(values = colors) +  # Custom colors for isolates
  geom_rect(data = genomic_regions, aes(xmin = Start, xmax = End, ymin = -Inf, ymax = Inf, fill = Region), 
            inherit.aes = FALSE, alpha = 0.2) +  # Add genomic regions with legend, no aesthetics inheritance
  scale_fill_manual(values = genomic_colors, name = "Genomic Region") +  # Genomic regions legend
  labs(x = "Genomic Position", y = "Isolate") +  # Removed plot title
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  
  # Adjust Y-axis to align isolate names with Shared SNPs
  scale_y_continuous(breaks = as.numeric(factor(shared_snp_long$Isolate)),
                     labels = shared_snp_long$Isolate)  # Align labels with Shared SNP positions

# Save the plot as a JPEG file
ggsave("NSNPs_P1_P10_Shared_Genomic_Region_Legend_Plot.jpeg", plot = snp_plot, width = 12, height = 8, dpi = 300)

# Print the plot to view
print(snp_plot)


