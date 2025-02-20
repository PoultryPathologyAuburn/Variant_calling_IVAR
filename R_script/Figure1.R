#Scatterplot showing each isolate has its unique SNPs for P1 shown just above its unique SNPs for P10
#Likewise Shared SNPs between P1 and P10 plotted just below P10 SNPs in the sam eplot

#########
# Set working directory
setwd("~/Documents/aminoacidanalysis/bam/iVaroutput/parameterchange/verification")

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)

# Define genomic regions (including HN and other regions)
genomic_regions <- data.frame(
  Region = c("NP", "P", "M", "F", "HN", "L"),
  Start = c(56, 1804, 3256, 4498, 6321, 8370),
  End = c(1792, 3244, 4487, 6279, 8312, 15073)
)

# Load your unique SNP positions files for P1, P10, and the shared SNPs
p1_snp_data <- read_excel("Uni_SNPs_P1.xlsx")
p10_snp_data <- read_excel("Uni_SNPs_P10.xlsx")
shared_snp_data <- read_excel("Shared_SNPs.xlsx")

# Convert the data to long format for P1, P10, and shared SNPs
p1_snp_long <- p1_snp_data %>%
  pivot_longer(cols = everything(), names_to = "Isolate", values_to = "POS") %>%
  mutate(Passage = "P1", Y_Position = 1, Shape = "Unique P1 SNPs")  # Label for P1

p10_snp_long <- p10_snp_data %>%
  pivot_longer(cols = everything(), names_to = "Isolate", values_to = "POS") %>%
  mutate(Passage = "P10", Y_Position = 0, Shape = "Unique P10 SNPs")  # Label for P10

shared_snp_long <- shared_snp_data %>%
  pivot_longer(cols = everything(), names_to = "Isolate", values_to = "POS") %>%
  mutate(Passage = "Shared P1P10", Y_Position = -1, Shape = "Shared SNPs between P1P10")  # Label for shared SNPs

# Convert Isolate names in the data to have the first letter capitalized
p1_snp_long <- p1_snp_long %>% mutate(Isolate = tools::toTitleCase(Isolate))
p10_snp_long <- p10_snp_long %>% mutate(Isolate = tools::toTitleCase(Isolate))
shared_snp_long <- shared_snp_long %>% mutate(Isolate = tools::toTitleCase(Isolate))

# Define colors for each isolate (ensure these match the isolate names in your data)
colors <- c("Iso1" = "blue", "Iso2" = "red", "Iso3" = "green", "Iso4" = "orange", "Iso5" = "purple", "Iso6" = "brown", "Iso7" = "pink")

# Remove any rows with NA values
p1_snp_long <- p1_snp_long %>% filter(!is.na(POS))
p10_snp_long <- p10_snp_long %>% filter(!is.na(POS))
shared_snp_long <- shared_snp_long %>% filter(!is.na(POS))

# Adjust the Y positions for the isolates to group P1, P10, and shared SNPs together for each Isolate
p1_snp_long <- p1_snp_long %>% mutate(Y_Position = as.numeric(factor(Isolate))*3 + Y_Position)
p10_snp_long <- p10_snp_long %>% mutate(Y_Position = as.numeric(factor(Isolate))*3 + Y_Position)
shared_snp_long <- shared_snp_long %>% mutate(Y_Position = as.numeric(factor(Isolate))*3 + Y_Position)

# Combine the data for plotting
combined_snp_data <- bind_rows(p1_snp_long, p10_snp_long, shared_snp_long)

# Plot SNPs for P1, P10, and shared SNPs, grouped by Isolate
snp_plot_combined <- ggplot() +
  geom_rect(data = genomic_regions, aes(xmin = Start, xmax = End, ymin = -Inf, ymax = Inf, fill = Region), alpha = 0.2) +
  geom_point(data = combined_snp_data, aes(x = POS, y = Y_Position, color = Isolate, shape = Shape), size = 2) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = c("NP" = "lightblue", "P" = "lightgreen", "M" = "lightpink", "F" = "lightyellow", "HN" = "lightcyan", "L" = "lightgrey")) +
  labs(x = "Genomic Position", y = "Isolate", color = "Isolate", fill = "Genomic Region", shape = "SNP Type") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = genomic_regions$Start, labels = genomic_regions$Region) +
  scale_y_continuous(breaks = seq(min(combined_snp_data$Y_Position), max(combined_snp_data$Y_Position), by = 3), labels = unique(p1_snp_long$Isolate)) +
  scale_shape_manual(values = c("Unique P1 SNPs" = 15, "Unique P10 SNPs" = 17, "Shared SNPs between P1P10" = 16))

# Save the combined SNP plot as JPEG
ggsave("SNPs_plot_combined_P1_P10_Shared.jpeg", plot = snp_plot_combined, width = 10, height = 8, dpi = 300)

# Print the plot to view
print(snp_plot_combined)

