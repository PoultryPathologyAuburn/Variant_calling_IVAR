
#This code contains
#1)Shared/Unique SNPs + Indels + counts
#2)Mean density stats (Kruskal + Dunn) for SNPs AND nSNPs 
#3)Shared/Unique nSNPs + summary

############################################################
## Project 1
##   Inputs:   Variants/*.tsv
##   Outputs:  Variants/REPRO_Project1_2026/...
##
## PARAMETERS:
##   - PASS == TRUE
##   - SNP vs Indel based on nchar(REF/ALT)
##   - Shared/Unique by POS+REF+ALT
##   - Same region_map
##   - nSNP: REF/ALT length 1, ALT_AA not NA, REF_AA != ALT_AA
##
## NOTE (FIXED):
##   NDV genome orientation: positions 1-121 are 3'UTR (leader),
##   positions 14996-15073 are 5'UTR (trailer).
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(writexl)
  library(tibble)   # tribble()
  library(tidyr)
  library(readxl)   # for Meandensity.xlsx section
  library(FSA)      # dunnTest
})

## ---- PROJECT 1 INPUT PATH ----
setwd("~/Documents/aminoacidanalysis/bam/iVaroutput/parameterchange/M1try/Variants")
data_dir <- getwd()

## ---- SAFE REPRO OUTPUT ROOT ----
repro_root <- file.path(getwd(), "REPRO_Project1_2026")
dir.create(repro_root, showWarnings = FALSE)

## ---- PROJECT 1 isolate map (7 isolates) ----
iso_map <- list(
  iso1  = "Iso1",
  iso5  = "Iso2",
  iso6  = "Iso3",
  iso7  = "Iso4",
  iso8  = "Iso5",
  iso9  = "Iso6",
  iso10 = "Iso7"
)
isolates <- names(iso_map)

## ---- Genomic region map (UTRs FIXED; coordinates UNCHANGED) ----
region_map <- tribble(
  ~Region, ~Start, ~End,
  "3'UTR", 1, 121,          
  "NP", 122, 1591,
  "NP_P", 1592, 1886,
  "P", 1887, 3074,
  "P_M", 3075, 3289,
  "M", 3290, 4384,
  "M_F", 4385, 4543,
  "F", 4544, 6205,
  "F_HN", 6206, 6411,
  "HN", 6412, 8145,
  "HN_L", 8146, 8380,
  "L", 8381, 14995,
  "5'UTR", 14996, 15073     
)

annotate_region <- function(pos) {
  hit <- region_map %>% filter(pos >= Start & pos <= End)
  if (nrow(hit) == 1) return(hit$Region)
  NA_character_
}

cat("\n✅ Inputs:", data_dir, "\n")
cat("✅ Outputs:", repro_root, "\n\n")

######################################################################
## CHUNK 1) SNPs & Indels (Shared, Unique P1, Unique P10) + summary
######################################################################
cat("===== CHUNK 1: Shared/Unique SNPs & Indels =====\n")

out1_dir <- file.path(repro_root, "PASS_Shared_Unique_SNPs_Indels_Project1")
dir.create(out1_dir, showWarnings = FALSE)

snp_results    <- list()
indel_results  <- list()
summary_counts <- data.frame()

for (iso in isolates) {
  file_p1  <- file.path(data_dir, paste0(iso, "p1.tsv"))
  file_p10 <- file.path(data_dir, paste0(iso, "p10.tsv"))
  
  if (file.exists(file_p1) & file.exists(file_p10)) {
    
    df_p1  <- read_tsv(file_p1,  show_col_types = FALSE) %>% filter(PASS == TRUE)
    df_p10 <- read_tsv(file_p10, show_col_types = FALSE) %>% filter(PASS == TRUE)
    
    df_p1$Genomic_Region  <- sapply(df_p1$POS,  annotate_region)
    df_p10$Genomic_Region <- sapply(df_p10$POS, annotate_region)
    
    ## SNPs
    snp_p1  <- df_p1  %>% filter(nchar(REF) == 1 & nchar(ALT) == 1)
    snp_p10 <- df_p10 %>% filter(nchar(REF) == 1 & nchar(ALT) == 1)
    
    shared_snps     <- inner_join(snp_p1, snp_p10, by = c("POS", "REF", "ALT"))
    unique_p1_snps  <- anti_join(snp_p1, snp_p10, by = c("POS", "REF", "ALT"))
    unique_p10_snps <- anti_join(snp_p10, snp_p1,  by = c("POS", "REF", "ALT"))
    
    ## Indels
    indel_p1  <- df_p1  %>% filter(nchar(REF) != 1 | nchar(ALT) != 1)
    indel_p10 <- df_p10 %>% filter(nchar(REF) != 1 | nchar(ALT) != 1)
    
    shared_indels     <- inner_join(indel_p1, indel_p10, by = c("POS", "REF", "ALT"))
    unique_p1_indels  <- anti_join(indel_p1, indel_p10, by = c("POS", "REF", "ALT"))
    unique_p10_indels <- anti_join(indel_p10, indel_p1,  by = c("POS", "REF", "ALT"))
    
    iso_label <- iso_map[[iso]]
    
    snp_results[[paste0(iso_label, "_Shared_SNPs")]]     <- shared_snps
    snp_results[[paste0(iso_label, "_Unique_P1_SNPs")]]  <- unique_p1_snps
    snp_results[[paste0(iso_label, "_Unique_P10_SNPs")]] <- unique_p10_snps
    
    indel_results[[paste0(iso_label, "_Shared_Indels")]]     <- shared_indels
    indel_results[[paste0(iso_label, "_Unique_P1_Indels")]]  <- unique_p1_indels
    indel_results[[paste0(iso_label, "_Unique_P10_Indels")]] <- unique_p10_indels
    
    summary_counts <- rbind(summary_counts, data.frame(
      Isolate          = iso_label,
      Shared_SNPs      = nrow(shared_snps),
      Unique_P1_SNPs   = nrow(unique_p1_snps),
      Unique_P10_SNPs  = nrow(unique_p10_snps),
      Shared_Indels    = nrow(shared_indels),
      Unique_P1_Indels = nrow(unique_p1_indels),
      Unique_P10_Indels= nrow(unique_p10_indels)
    ))
    
  } else {
    warning(paste("Missing files for", iso))
  }
}

write_xlsx(snp_results,   file.path(out1_dir, "Shared_Unique_SNPs_Project1.xlsx"))
write_xlsx(indel_results, file.path(out1_dir, "Shared_Unique_Indels_Project1.xlsx"))
write.csv(summary_counts, file.path(out1_dir, "Shared_Unique_Variant_Counts_Project1.csv"), row.names = FALSE)

cat("✅ Chunk 1 outputs saved in:\n  ", out1_dir, "\n\n", sep="")
print(summary_counts)

######################################################################
## CHUNK 2) Mean density stats (Kruskal-Wallis + Dunn) for SNPs AND nSNPs
######################################################################

cat("\n===== CHUNK 2: Mean density stats (SNP + nSNP) =====\n")

file_path <- file.path(repro_root, "Meandensity_Re.xlsx")
out2_dir  <- out1_dir  

if (!file.exists(file_path)) {
  warning(paste("Meandensity_Re.xlsx NOT found at:", file_path))
  warning("Skipping Chunk 2. (Missing Excel file.)")
} else {
  
  analyze_sheet <- function(sheet_name, label) {
    df <- read_excel(file_path, sheet = sheet_name)
    
    # expects columns: Gene, Iso1..Iso7 (these can be densities or counts)
    long_df <- df %>%
      select(Gene, Iso1:Iso7) %>%
      pivot_longer(cols = starts_with("Iso"),
                   names_to = "Isolate", values_to = "Value") %>%
      mutate(Value = suppressWarnings(as.numeric(Value))) %>%
      filter(!is.na(Value))
    
    # Kruskal–Wallis
    kw <- kruskal.test(Value ~ Gene, data = long_df)
    
    write.table(
      data.frame(Sheet = label, Statistic = as.numeric(kw$statistic), P_Value = as.numeric(kw$p.value)),
      file = file.path(out2_dir, paste0("Kruskal_Statistic_", label, "_Project1.txt")),
      row.names = FALSE, quote = FALSE, sep = "\t"
    )
    
    # Dunn post-hoc (Bonferroni)
    dunn <- FSA::dunnTest(Value ~ Gene, data = long_df, method = "bonferroni")
    write.csv(dunn$res,
              file = file.path(out2_dir, paste0("Dunn_PostHoc_", label, "_Project1.csv")),
              row.names = FALSE)
    
    invisible(TRUE)
  }
  
  ## ---- SNP density sheets ----
  analyze_sheet("SharedSNPs",    "Shared_SNPs")
  analyze_sheet("UniqueP1SNPs",  "UniqueP1_SNPs")
  analyze_sheet("UniqueP10SNPs", "UniqueP10_SNPs")
  
  ## ---- nSNP density sheets ✅ ADDED ----
  analyze_sheet("SharednSNPs",    "Shared_nSNPs")
  analyze_sheet("UniqueP1nSNPs",  "UniqueP1_nSNPs")
  analyze_sheet("UniqueP10nSNPs", "UniqueP10_nSNPs")
  
  cat("✅ Chunk 2 outputs saved in:\n  ", out2_dir, "\n", sep = "")
}





######################################################################
## CHUNK 3) Shared/Unique nSNPs + summary
######################################################################
cat("\n===== CHUNK 3: Shared/Unique nSNPs =====\n")

out3_dir <- file.path(repro_root, "Results_iVarTSV_nSNPs_Project1")
dir.create(out3_dir, showWarnings = FALSE)

xlsx_output <- file.path(out3_dir, "Shared_Unique_nSNPs_Project1.xlsx")
summary_csv <- file.path(out3_dir, "NSNP_P1_P10_Shared_Unique_Summary_Project1.csv")

nsnp_results <- list()
summary_list <- list()

for (iso in isolates) {
  file_p1  <- file.path(data_dir, paste0(iso, "p1.tsv"))
  file_p10 <- file.path(data_dir, paste0(iso, "p10.tsv"))
  
  if (file.exists(file_p1) && file.exists(file_p10)) {
    
    df_p1  <- read_tsv(file_p1,  show_col_types = FALSE) %>% filter(PASS == TRUE)
    df_p10 <- read_tsv(file_p10, show_col_types = FALSE) %>% filter(PASS == TRUE)
    
    nsnp_p1 <- df_p1 %>%
      filter(nchar(REF) == 1, nchar(ALT) == 1, !is.na(ALT_AA), REF_AA != ALT_AA) %>%
      mutate(Genomic_Region = sapply(POS, annotate_region))
    
    nsnp_p10 <- df_p10 %>%
      filter(nchar(REF) == 1, nchar(ALT) == 1, !is.na(ALT_AA), REF_AA != ALT_AA) %>%
      mutate(Genomic_Region = sapply(POS, annotate_region))
    
    shared     <- inner_join(nsnp_p1, nsnp_p10, by = c("POS", "REF", "ALT"))
    unique_p1  <- anti_join(nsnp_p1, nsnp_p10, by = c("POS", "REF", "ALT"))
    unique_p10 <- anti_join(nsnp_p10, nsnp_p1,  by = c("POS", "REF", "ALT"))
    
    iso_label <- iso_map[[iso]]
    
    nsnp_results[[paste0(iso_label, "_Shared_nSNPs")]]     <- shared
    nsnp_results[[paste0(iso_label, "_Unique_P1_nSNPs")]]  <- unique_p1
    nsnp_results[[paste0(iso_label, "_Unique_P10_nSNPs")]] <- unique_p10
    
    summary_list[[iso_label]] <- data.frame(
      Isolate         = iso_label,
      Shared_nSNPs    = nrow(shared),
      Unique_P1_nSNPs = nrow(unique_p1),
      Unique_P10_nSNPs= nrow(unique_p10)
    )
  } else {
    warning(paste("Missing files for", iso))
  }
}

nsnp_summary <- bind_rows(summary_list)
write.csv(nsnp_summary, summary_csv, row.names = FALSE)
write_xlsx(nsnp_results, xlsx_output)

cat("✅ Chunk 3 outputs saved in:\n  ", out3_dir, "\n\n", sep="")
print(nsnp_summary)

cat("\n🎉 ALL DONE. Everything saved under:\n  ", repro_root, "\n", sep = "")







############################################################
## Project 1
## (A) Shared nSNPs per isolate + recurrence + UpSet plots
## (B) Mapping of Shared / Unique P1 / Unique P10 nSNPs
## (C) Recurrence of UNIQUE P1 & UNIQUE P10 nSNPs + UpSet
## (D) Significant Δ-frequency nSNPs (|Δ| > 0.2) + Fisher + BH
##
## SAFE OUTPUT MODE:
##   Inputs:  Variants/*.tsv
##   Output:  Variants/REPRO_Project1_2026/...
##
## PARAMETERS UNCHANGED:
##  - PASS == TRUE
##  - SNP only: nchar(REF)==1 & nchar(ALT)==1
##  - nSNP: REF_AA and ALT_AA not NA, REF_AA != ALT_AA
##
## NOTE (FIXED UTR LABELS):
##  NDV genome orientation:
##   - positions 1–121 are 3'UTR (leader)
##   - positions 14996–15073 are 5'UTR (trailer)
##  (Coordinates unchanged; only labels corrected.)
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)  
  library(UpSetR)
  library(ggplot2)
})

setwd("~/Documents/aminoacidanalysis/bam/iVaroutput/parameterchange/M1try/Variants")
data_dir <- getwd()

repro_root <- file.path(getwd(), "REPRO_Project1_2026")
dir.create(repro_root, showWarnings = FALSE)

iso_map <- list(
  iso1  = "Iso1",
  iso5  = "Iso2",
  iso6  = "Iso3",
  iso7  = "Iso4",
  iso8  = "Iso5",
  iso9  = "Iso6",
  iso10 = "Iso7"
)
isolates <- names(iso_map)

## ---- Region map 
region_map <- tribble(
  ~Region, ~Start, ~End,
  "3'UTR", 1, 121,         
  "NP", 122, 1591,
  "NP_P", 1592, 1886,
  "P", 1887, 3074,
  "P_M", 3075, 3289,
  "M", 3290, 4384,
  "M_F", 4385, 4543,
  "F", 4544, 6205,
  "F_HN", 6206, 6411,
  "HN", 6412, 8145,
  "HN_L", 8146, 8380,
  "L", 8381, 14995,
  "5'UTR", 14996, 15073     
) %>% mutate(Mid = (Start + End)/2)

annotate_region <- function(pos) {
  hit <- region_map %>% filter(pos >= Start & pos <= End)
  if (nrow(hit) == 1) hit$Region else NA_character_
}

############################################################
## Helper: extract nSNPs from a TSV data.frame (PASS only)
############################################################
get_nsnp <- function(df){
  df %>%
    filter(
      PASS == TRUE,
      nchar(REF) == 1, nchar(ALT) == 1,
      !is.na(REF_AA), !is.na(ALT_AA),
      REF_AA != ALT_AA
    ) %>%
    mutate(unique_id = paste(POS, REF, ALT, REF_AA, ALT_AA, sep = "_"))
}

######################################################################
## (A) SHARED nSNPs per isolate + recurrence + UpSet
######################################################################
cat("\n===== (A) Shared nSNPs per isolate + recurrence + UpSet =====\n")

outA <- file.path(repro_root, "Results_Shared_nSNPs_Project1")
dir.create(outA, showWarnings = FALSE)
dir.create(file.path(outA, "per_isolate"), showWarnings = FALSE)
dir.create(file.path(outA, "plots"), showWarnings = FALSE)

extract_shared_nsnp <- function(iso) {
  f1  <- file.path(data_dir, paste0(iso, "p1.tsv"))
  f10 <- file.path(data_dir, paste0(iso, "p10.tsv"))
  if (!file.exists(f1) || !file.exists(f10)) {
    message("Missing files for ", iso); return(NULL)
  }
  
  p1  <- read_tsv(f1,  show_col_types = FALSE) %>% filter(PASS == TRUE)
  p10 <- read_tsv(f10, show_col_types = FALSE) %>% filter(PASS == TRUE)
  
  p1_nsnp  <- get_nsnp(p1)
  p10_nsnp <- get_nsnp(p10)
  
  shared_ids <- intersect(p1_nsnp$unique_id, p10_nsnp$unique_id)
  
  shared <- p1_nsnp %>%
    filter(unique_id %in% shared_ids) %>%
    mutate(Genomic_Region = sapply(POS, annotate_region))
  
  iso_label <- iso_map[[iso]]
  
  out_df <- shared %>%
    select(Genomic_Region, POS, REF, ALT, REF_AA, ALT_AA, ALT_FREQ, unique_id) %>%
    mutate(Isolate = iso_label)
  
  write_tsv(out_df, file.path(outA, "per_isolate",
                              paste0(iso_label, "_Shared_nSNPs_simplified.tsv")))
  message("Processed ", iso_label, " - shared nSNPs: ", nrow(out_df))
  
  out_df
}

shared_all <- bind_rows(lapply(isolates, extract_shared_nsnp))

shared_summary <- shared_all %>%
  group_by(unique_id) %>%
  summarise(
    Isolate_Count = n_distinct(Isolate),
    Isolates = paste(sort(unique(Isolate)), collapse = ","),
    POS = first(POS), REF = first(REF), ALT = first(ALT),
    REF_AA = first(REF_AA), ALT_AA = first(ALT_AA),
    Genomic_Region = first(Genomic_Region),
    ALT_FREQ = mean(ALT_FREQ, na.rm = TRUE),
    .groups = "drop"
  ) %>% arrange(desc(Isolate_Count), POS)

rec_file <- file.path(outA, "Shared_nSNPs_Isolate_Recurrence_Project1.tsv")
write.table(shared_summary, rec_file, sep = "\t", quote = FALSE, row.names = FALSE)

shared_all7 <- shared_summary %>% filter(Isolate_Count == 7)
all7_file <- file.path(outA, "Shared_nSNPs_All7_Isolates_Project1.tsv")
write.table(shared_all7, all7_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("✅ Outputs:\n",
    "• ", rec_file, "\n",
    "• ", all7_file, "\n",
    "• per-isolate files in ", file.path(outA, "per_isolate"), "\n", sep = "")

## ---- UpSet ----
plot_dirA <- file.path(outA, "plots")
df <- read.delim(rec_file, stringsAsFactors = FALSE)

df_bin <- df %>%
  separate_rows(Isolates, sep = ",") %>%
  mutate(Isolates = trimws(Isolates), value = 1) %>%
  distinct(unique_id, Isolates, .keep_all = TRUE) %>%      # de-dup pairs
  pivot_wider(names_from = Isolates, values_from = value,
              values_fill = 0, values_fn = ~1) %>%          # robust if dupes
  as.data.frame()

rownames(df_bin) <- make.unique(df$unique_id)

isolate_order <- paste0("Iso", 1:7)   # bottom -> top order
missing_cols <- setdiff(isolate_order, colnames(df_bin))
if (length(missing_cols)) df_bin[, missing_cols] <- 0
df_bin <- df_bin[, isolate_order, drop = FALSE]

png(file.path(plot_dirA, "Shared_nSNPs_UpSetPlot_Project1_Publication.png"),
    width = 2200, height = 1400, res = 300)
UpSetR::upset(df_bin,
              sets = isolate_order, keep.order = TRUE, order.by = "freq",
              sets.x.label = "Number of Shared nSNPs",
              mainbar.y.label = "nSNP Intersections Across Isolates",
              point.size = 2.5, line.size = 1, text.scale = 1.2)
dev.off()

png(file.path(plot_dirA, "Shared_nSNPs_UpSetPlot_Project1_AuburnTheme.png"),
    width = 2200, height = 1400, res = 300)
UpSetR::upset(df_bin,
              sets = isolate_order, keep.order = TRUE, order.by = "freq",
              sets.x.label = "Number of Shared nSNPs",
              mainbar.y.label = "nSNP Intersections Across Isolates",
              point.size = 3, line.size = 1.2,
              main.bar.color = "#0C2340", matrix.color = "#E87722",
              sets.bar.color = "#0C2340", text.scale = 1.6)
dev.off()

cat("\n📊 UpSet plots saved to: ", plot_dirA, "\n", sep = "")


######################################################################
## (B) MAPPING: Shared / Unique P1 / Unique P10 SNPs and nSNPs across genome
######################################################################
cat("\n===== (B) Mapping Shared/Unique P1/Unique P10 SNPs and nSNPs =====\n")

###########################
# Figure 2 (Project 1): SNP position map (Shared vs Unique P1 vs Unique P10)
# FIXED for Project 1:
#  - 7 isolates (Iso1..Iso7)
#  - UTR labels corrected (3'UTR at 1–121, 5'UTR at 14996–15073)
#  - Saves outputs inside REPRO_Project1_2026

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

# ---- INPUT (Excel) ----
file_path <- file.path(repro_root, "Figures2N3_Project1.xlsx")

# ---- OUTPUT folder (SAFE) ----
out_dir <- file.path(repro_root, "Results_Figure2_Project1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Genomic regions (Project 1; UTR labels FIXED, coordinates unchanged) ----
genomic_regions <- data.frame(
  Region = c("3'UTR", "NP", "NP_P", "P", "P_M", "M", "M_F", "F", "F_HN", "HN", "HN_L", "L", "5'UTR"),
  Start  = c(1, 122, 1592, 1887, 3075, 3290, 4385, 4544, 6206, 6412, 8146, 8381, 14996),
  End    = c(121, 1591, 1886, 3074, 3289, 4384, 4543, 6205, 6411, 8145, 8380, 14995, 15073)
)
genomic_regions$Mid <- (genomic_regions$Start + genomic_regions$End)/2

# ---- Region colors ----
region_colors <- c(
  "NP" = "lightblue", "P" = "lightgreen", "M" = "lightpink",
  "F" = "lightyellow", "HN" = "lightcyan", "L" = "lightgrey",
  "3'UTR" = "#FDE725", "NP_P" = "#5DC863", "P_M" = "#21918C",
  "M_F" = "#3B528B", "F_HN" = "#440154", "HN_L" = "#9E0142", "5'UTR" = "#F46D43"
)

# ---- Isolate colors (Project 1: Iso1..Iso7) ----
isolate_colors_7 <- c(
  "Iso1" = "blue", "Iso2" = "red", "Iso3" = "green", "Iso4" = "orange",
  "Iso5" = "purple", "Iso6" = "brown", "Iso7" = "pink"
)

# ---- Sheet mapping ----
sheet_mapping <- list(
  "Shared SNPs between P1P10" = "SharedSNPs",
  "Unique P1 SNPs"            = "UniqueP1SNPs",
  "Unique P10 SNPs"           = "UniqueP10SNPs"
)

# Helper to safely read a sheet even if it has accidental trailing spaces
read_sheet_fuzzy <- function(path, sheet_name) {
  sheets <- excel_sheets(path)
  hit <- sheets[trimws(sheets) == trimws(sheet_name)]
  if (length(hit) == 0) stop("Sheet not found (after trimming): ", sheet_name)
  read_excel(path, sheet = hit[1])
}

# ---- Read and reshape (expects columns Iso1..Iso7; values are POS) ----
iso_levels <- paste0("Iso", 1:7)

all_data <- lapply(names(sheet_mapping), function(type) {
  sheet <- sheet_mapping[[type]]
  df <- read_sheet_fuzzy(file_path, sheet)
  
  df %>%
    pivot_longer(cols = everything(), names_to = "Isolate", values_to = "POS") %>%
    filter(!is.na(POS)) %>%
    mutate(
      Isolate = trimws(Isolate),
      Type = type
    ) %>%
    filter(Isolate %in% iso_levels)
}) %>% bind_rows()

# ---- Add plotting offsets ----
all_data <- all_data %>%
  mutate(
    Y_Offset  = case_when(
      Type == "Unique P1 SNPs" ~  1,
      Type == "Unique P10 SNPs" ~ 0,
      Type == "Shared SNPs between P1P10" ~ -1
    ),
    Isolate = factor(Isolate, levels = iso_levels),
    Type    = factor(Type, levels = c("Shared SNPs between P1P10", "Unique P1 SNPs", "Unique P10 SNPs")),
    Y_Position = as.numeric(Isolate) * 3 + Y_Offset
  )

# ---- Plot (SNPs) ----
p_snp <- ggplot() +
  geom_rect(
    data = genomic_regions,
    aes(xmin = Start, xmax = End, ymin = -Inf, ymax = Inf, fill = Region),
    alpha = 0.25
  ) +
  geom_point(
    data = all_data,
    aes(x = POS, y = Y_Position, color = Isolate, shape = Type),
    size = 2
  ) +
  scale_shape_manual(values = c(
    "Shared SNPs between P1P10" = 16,
    "Unique P1 SNPs" = 15,
    "Unique P10 SNPs" = 17
  )) +
  scale_color_manual(values = isolate_colors_7) +
  scale_fill_manual(values = region_colors) +
  scale_x_continuous(
    breaks = genomic_regions$Mid,
    labels = genomic_regions$Region,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    breaks = seq(3, 3*7, by = 3),
    labels = iso_levels
  ) +
  labs(
    x = "Genomic Position",
    y = "Isolate",
    color = "Isolates",
    shape = "SNP Types",
    fill = "Genomic Regions",
    title = "Figure 2 (Project 1): SNP Position Map (Shared vs Unique P1 vs Unique P10)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

# ---- Save outputs (SNPs) ----
ggsave(
  filename = file.path(out_dir, "Figure2_SNP_Position_Map_Project1.jpeg"),
  plot = p_snp, width = 12, height = 8, dpi = 300
)
write.csv(all_data, file.path(out_dir, "Figure2_SNP_Position_Data_Long_Project1.csv"), row.names = FALSE)

print(p_snp)
cat("✅ Figure 2 saved to:", out_dir, "\n")


###########################
# Figure 3 (Project 1): nSNP position map (Shared vs Unique P1 vs Unique P10)

outB <- file.path(repro_root, "Results_Figure_nSNP_Map_Project1")
dir.create(outB, showWarnings = FALSE)

all_nsnps <- list()

for (iso in isolates){
  f1  <- file.path(data_dir, paste0(iso, "p1.tsv"))
  f10 <- file.path(data_dir, paste0(iso, "p10.tsv"))
  if (!file.exists(f1) || !file.exists(f10)) next
  
  p1  <- read_tsv(f1,  show_col_types = FALSE)
  p10 <- read_tsv(f10, show_col_types = FALSE)
  
  nsnp_p1  <- get_nsnp(p1)
  nsnp_p10 <- get_nsnp(p10)
  
  shared_ids <- intersect(nsnp_p1$unique_id, nsnp_p10$unique_id)
  
  shared     <- nsnp_p1  %>% filter(unique_id %in% shared_ids)
  unique_p1  <- nsnp_p1  %>% filter(!unique_id %in% shared_ids)
  unique_p10 <- nsnp_p10 %>% filter(!unique_id %in% shared_ids)
  
  iso_label <- iso_map[[iso]]
  
  all_nsnps[[iso_label]] <- bind_rows(
    shared     %>% transmute(Isolate = iso_label, Type = "Shared nSNPs", POS),
    unique_p1  %>% transmute(Isolate = iso_label, Type = "Unique P1 nSNPs", POS),
    unique_p10 %>% transmute(Isolate = iso_label, Type = "Unique P10 nSNPs", POS)
  )
}

plot_df <- bind_rows(all_nsnps) %>%
  mutate(
    Isolate = factor(Isolate, levels = paste0("Iso", 1:7)),
    Type = factor(Type, levels = c("Shared nSNPs", "Unique P1 nSNPs", "Unique P10 nSNPs")),
    Y_Offset = case_when(
      Type == "Unique P1 nSNPs"  ~  1,
      Type == "Unique P10 nSNPs" ~  0,
      TRUE                       ~ -1
    ),
    Y_Position = as.numeric(Isolate) * 3 + Y_Offset
  )

# Use the SAME region_map but make sure it has FIXED UTR labels in your main script:
# region_map must have 3'UTR at 1-121 and 5'UTR at 14996-15073
# region_colors already matches those names.

p_nsnp <- ggplot() +
  geom_rect(
    data = region_map,
    aes(xmin = Start, xmax = End, ymin = -Inf, ymax = Inf, fill = Region),
    alpha = 0.25
  ) +
  geom_point(
    data = plot_df,
    aes(x = POS, y = Y_Position, color = Isolate, shape = Type),  # ✅ FIX: color by isolate
    size = 2
  ) +
  scale_shape_manual(values = c(
    "Shared nSNPs" = 16,
    "Unique P1 nSNPs" = 15,
    "Unique P10 nSNPs" = 17
  )) +
  scale_color_manual(values = isolate_colors_7) +                  # ✅ FIX: apply isolate colors
  scale_fill_manual(values = region_colors) +
  scale_x_continuous(
    breaks = region_map$Mid,
    labels = region_map$Region,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    breaks = seq(3, 3*7, by = 3),
    labels = paste0("Iso", 1:7)
  ) +
  labs(
    x = "Genomic Position",
    y = "Isolate",
    color = "Isolate",                                             # ✅ legend
    shape = "nSNP Type",
    fill  = "Genomic Region",
    title = "Figure 3 (Project 1): nSNP Position Map (Shared vs Unique P1 vs Unique P10)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

png_file <- file.path(outB, "Figure3_nSNP_Position_Map_Project1.png")
pdf_file <- file.path(outB, "Figure3_nSNP_Position_Map_Project1.pdf")
ggsave(png_file, p_nsnp, width = 12, height = 8, dpi = 300)
ggsave(pdf_file, p_nsnp, width = 12, height = 8)
write.csv(plot_df, file.path(outB, "Figure3_nSNP_Position_Data_Long_Project1.csv"), row.names = FALSE)

print(p_nsnp)
cat("\n✅ nSNP mapping plot saved to:\n  ", outB, "\n", sep = "")

cat("\n🎉 DONE (B). Outputs saved under:\n  ", repro_root, "\n", sep = "")



######################################################################
## (C) Recurrence of UNIQUE P1 & UNIQUE P10 nSNPs across isolates + UpSet
######################################################################
cat("\n===== (C) Unique P1 & Unique P10 nSNP recurrence + UpSet =====\n")

outC <- file.path(repro_root, "Results_Unique_nSNPs_Recurrence_Project1")
dir.create(outC, showWarnings = FALSE)
dir.create(file.path(outC, "plots"), showWarnings = FALSE)

get_unique_p1_p10_nsnp <- function(iso){
  f1  <- file.path(data_dir, paste0(iso, "p1.tsv"))
  f10 <- file.path(data_dir, paste0(iso, "p10.tsv"))
  if (!file.exists(f1) || !file.exists(f10)) return(NULL)
  
  p1  <- read_tsv(f1,  show_col_types = FALSE)
  p10 <- read_tsv(f10, show_col_types = FALSE)
  
  nsnp_p1  <- get_nsnp(p1)  %>% mutate(Genomic_Region = sapply(POS, annotate_region))
  nsnp_p10 <- get_nsnp(p10) %>% mutate(Genomic_Region = sapply(POS, annotate_region))
  
  shared_ids <- intersect(nsnp_p1$unique_id, nsnp_p10$unique_id)
  
  unique_p1  <- nsnp_p1  %>% filter(!unique_id %in% shared_ids)
  unique_p10 <- nsnp_p10 %>% filter(!unique_id %in% shared_ids)
  
  iso_label <- iso_map[[iso]]
  
  list(
    p1 = unique_p1  %>% mutate(Isolate = iso_label),
    p10= unique_p10 %>% mutate(Isolate = iso_label)
  )
}

u_list <- lapply(isolates, get_unique_p1_p10_nsnp)
unique_p1_all  <- bind_rows(lapply(u_list, `[[`, "p1"))
unique_p10_all <- bind_rows(lapply(u_list, `[[`, "p10"))

summarize_recurrence <- function(df){
  df %>%
    group_by(unique_id) %>%
    summarise(
      Isolate_Count = n_distinct(Isolate),
      Isolates = paste(sort(unique(Isolate)), collapse = ","),
      POS = first(POS),
      REF = first(REF),
      ALT = first(ALT),
      REF_AA = first(REF_AA),
      ALT_AA = first(ALT_AA),
      Genomic_Region = first(Genomic_Region),
      ALT_FREQ = mean(as.numeric(ALT_FREQ), na.rm = TRUE),
      .groups = "drop"
    ) %>% arrange(desc(Isolate_Count), POS)
}

summary_uP1  <- summarize_recurrence(unique_p1_all)
summary_uP10 <- summarize_recurrence(unique_p10_all)

p1_rec_file  <- file.path(outC, "UniqueP1_nSNPs_Isolate_Recurrence_Project1.tsv")
p10_rec_file <- file.path(outC, "UniqueP10_nSNPs_Isolate_Recurrence_Project1.tsv")
write.table(summary_uP1,  p1_rec_file,  sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary_uP10, p10_rec_file, sep = "\t", quote = FALSE, row.names = FALSE)

write.table(summary_uP1  %>% filter(Isolate_Count == 7),
            file.path(outC, "UniqueP1_nSNPs_All7_Isolates_Project1.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary_uP10 %>% filter(Isolate_Count == 7),
            file.path(outC, "UniqueP10_nSNPs_All7_Isolates_Project1.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("✅ Recurrence tables saved in:\n  ", outC, "\n", sep="")

prepare_upset_matrix <- function(summary_tbl, iso_order){
  mat <- summary_tbl %>%
    separate_rows(Isolates, sep = ",") %>%
    mutate(Isolates = trimws(Isolates), value = 1) %>%
    distinct(unique_id, Isolates, .keep_all = TRUE) %>%   # de-dup pairs
    pivot_wider(names_from = Isolates, values_from = value,
                values_fill = 0, values_fn = ~1) %>%
    as.data.frame()
  
  rownames(mat) <- make.unique(summary_tbl$unique_id)
  
  missing_cols <- setdiff(iso_order, colnames(mat))
  if (length(missing_cols)) mat[, missing_cols] <- 0
  
  mat[, iso_order, drop = FALSE]
}

iso_order <- paste0("Iso", 1:7)

p1_mat  <- prepare_upset_matrix(summary_uP1,  iso_order)
p10_mat <- prepare_upset_matrix(summary_uP10, iso_order)

write.csv(p1_mat,  file.path(outC, "UniqueP1_nSNPs_UpSet_Matrix_Project1.csv"),  row.names = TRUE)
write.csv(p10_mat, file.path(outC, "UniqueP10_nSNPs_UpSet_Matrix_Project1.csv"), row.names = TRUE)

plot_dirC <- file.path(outC, "plots")

png(file.path(plot_dirC, "UniqueP1_nSNPs_UpSet_Project1_Publication.png"),
    width = 2200, height = 1400, res = 300)
UpSetR::upset(p1_mat,
              sets = iso_order, keep.order = TRUE, order.by = "freq",
              mainbar.y.label = "Unique P1 nSNPs",
              sets.x.label    = "nSNPs per Isolate",
              point.size = 2.5, line.size = 1, text.scale = 1.3)
dev.off()

png(file.path(plot_dirC, "UniqueP10_nSNPs_UpSet_Project1_Publication.png"),
    width = 2200, height = 1400, res = 300)
UpSetR::upset(p10_mat,
              sets = iso_order, keep.order = TRUE, order.by = "freq",
              mainbar.y.label = "Unique P10 nSNPs",
              sets.x.label    = "nSNPs per Isolate",
              point.size = 2.5, line.size = 1, text.scale = 1.3)
dev.off()

png(file.path(plot_dirC, "UniqueP1_nSNPs_UpSet_Project1_Auburn.png"),
    width = 2200, height = 1400, res = 300)
UpSetR::upset(p1_mat,
              sets = iso_order, keep.order = TRUE, order.by = "freq",
              mainbar.y.label = "Unique P1 nSNPs",
              sets.x.label    = "nSNPs per Isolate",
              point.size = 2.8, line.size = 1.2,
              main.bar.color = "#0C2340", matrix.color = "#E87722",
              sets.bar.color = "#0C2340", text.scale = 1.6)
dev.off()

png(file.path(plot_dirC, "UniqueP10_nSNPs_UpSet_Project1_Auburn.png"),
    width = 2200, height = 1400, res = 300)
UpSetR::upset(p10_mat,
              sets = iso_order, keep.order = TRUE, order.by = "freq",
              mainbar.y.label = "Unique P10 nSNPs",
              sets.x.label    = "nSNPs per Isolate",
              point.size = 2.8, line.size = 1.2,
              main.bar.color = "#0C2340", matrix.color = "#E87722",
              sets.bar.color = "#0C2340", text.scale = 1.6)
dev.off()

cat("📊 Unique nSNP UpSet plots saved to:\n  ", plot_dirC, "\n", sep = "")

######################################################################
## (D) Significant Δ-frequency nSNPs (|Δ| > 0.2) + Fisher + BH
## Project 1 version (7 isolates) - outputs in REPRO_Project1_2026
######################################################################
cat("\n===== (D) Δ-frequency nSNPs + Fisher tests + BH (Project 1) =====\n")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

output_dir_D <- file.path(repro_root, "Statistical_Results_DeltaFreq_nSNPs")
dir.create(output_dir_D, showWarnings = FALSE, recursive = TRUE)

region_boundaries <- data.frame(
  Gene  = c("NP","P","M","F","HN","L"),
  Start = c(122, 1887, 3290, 4544, 6412, 8381),
  End   = c(1591, 3074, 4384, 6205, 8145, 14995)
)

annotate_with_gene <- function(df, bounds) {
  if (is.null(df) || nrow(df) == 0) return(df)
  df$Gene <- NA_character_
  df$Rel_POS_nt <- NA_integer_
  for (i in seq_len(nrow(bounds))) {
    idx <- df$POS >= bounds$Start[i] & df$POS <= bounds$End[i]
    df$Gene[idx] <- bounds$Gene[i]
    df$Rel_POS_nt[idx] <- df$POS[idx] - bounds$Start[i] + 1
  }
  df
}

perform_fishers_test <- function(alt_p1, total_p1, alt_p10, total_p10) {
  if (is.na(total_p1) || is.na(total_p10)) {
    return(list(p_value = NA_real_, odds_ratio = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_))
  }
  if (is.na(alt_p1))  alt_p1  <- 0
  if (is.na(alt_p10)) alt_p10 <- 0
  
  ref_p1  <- total_p1  - alt_p1
  ref_p10 <- total_p10 - alt_p10
  
  tab <- matrix(c(alt_p1, ref_p1, alt_p10, ref_p10),
                nrow = 2, byrow = TRUE,
                dimnames = list(c("P1","P10"), c("Alt","Ref")))
  
  tryCatch({
    ft <- fisher.test(tab)
    list(
      p_value    = as.numeric(ft$p.value),
      odds_ratio = as.numeric(ft$estimate),
      ci_lower   = as.numeric(ft$conf.int[1]),
      ci_upper   = as.numeric(ft$conf.int[2])
    )
  }, error = function(e) {
    list(p_value = NA_real_, odds_ratio = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_)
  })
}

delta_thresh <- 0.2
min_coverage <- 50
alpha <- 0.05

all_sig_list   <- list()
all_stats_list <- list()

analyze_isolate_delta <- function(orig_iso, new_iso) {
  p1_file  <- file.path(data_dir, paste0(orig_iso, "p1.tsv"))
  p10_file <- file.path(data_dir, paste0(orig_iso, "p10.tsv"))
  if (!file.exists(p1_file) || !file.exists(p10_file)) return(NULL)
  
  message("Processing ", new_iso, " ...")
  
  df_p1 <- read_tsv(p1_file, show_col_types = FALSE) %>%
    filter(PASS == TRUE,
           nchar(REF) == 1, nchar(ALT) == 1,
           !is.na(ALT_AA), !is.na(REF_AA),
           ALT_AA != REF_AA) %>%
    select(REGION, POS, POS_AA, REF, ALT, REF_AA, ALT_AA,
           ALT_FREQ, TOTAL_DP, ALT_DP) %>%
    rename(ALT_FREQ_P1 = ALT_FREQ, TOTAL_DP_P1 = TOTAL_DP, ALT_DP_P1 = ALT_DP)
  
  df_p10 <- read_tsv(p10_file, show_col_types = FALSE) %>%
    filter(PASS == TRUE,
           nchar(REF) == 1, nchar(ALT) == 1,
           !is.na(ALT_AA), !is.na(REF_AA),
           ALT_AA != REF_AA) %>%
    select(REGION, POS, POS_AA, REF, ALT, REF_AA, ALT_AA,
           ALT_FREQ, TOTAL_DP, ALT_DP) %>%
    rename(ALT_FREQ_P10 = ALT_FREQ, TOTAL_DP_P10 = TOTAL_DP, ALT_DP_P10 = ALT_DP)
  
  df_combined <- full_join(
    df_p1, df_p10,
    by = c("REGION","POS","POS_AA","REF","ALT","REF_AA","ALT_AA")
  ) %>%
    mutate(
      ALT_FREQ_P1  = ifelse(is.na(ALT_FREQ_P1),  0, ALT_FREQ_P1),
      ALT_FREQ_P10 = ifelse(is.na(ALT_FREQ_P10), 0, ALT_FREQ_P10),
      ALT_DP_P1    = ifelse(is.na(ALT_DP_P1),    0, ALT_DP_P1),
      ALT_DP_P10   = ifelse(is.na(ALT_DP_P10),   0, ALT_DP_P10),
      Frequency_Change = ALT_FREQ_P10 - ALT_FREQ_P1,
      Isolate = new_iso,
      Direction = case_when(
        Frequency_Change > 0  ~ "Up in P10",
        Frequency_Change < 0  ~ "Down in P10",
        TRUE ~ "No change"
      )
    ) %>%
    filter(!is.na(TOTAL_DP_P1) & !is.na(TOTAL_DP_P10),
           TOTAL_DP_P1 >= min_coverage & TOTAL_DP_P10 >= min_coverage) %>%
    annotate_with_gene(region_boundaries)
  
  message("Performing Fisher's exact tests...")
  statistical_results <- df_combined %>%
    mutate(p_value = NA_real_, odds_ratio = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_)
  
  for (i in seq_len(nrow(statistical_results))) {
    r <- perform_fishers_test(
      statistical_results$ALT_DP_P1[i],  statistical_results$TOTAL_DP_P1[i],
      statistical_results$ALT_DP_P10[i], statistical_results$TOTAL_DP_P10[i]
    )
    statistical_results$p_value[i]    <- r$p_value
    statistical_results$odds_ratio[i] <- r$odds_ratio
    statistical_results$ci_lower[i]   <- r$ci_lower
    statistical_results$ci_upper[i]   <- r$ci_upper
  }
  
  statistical_results$p_adj <- p.adjust(statistical_results$p_value, method = "BH")
  
  statistical_results <- statistical_results %>%
    mutate(
      Statistically_Significant = ifelse(!is.na(p_adj) & p_adj < alpha, "Yes", "No"),
      Effect_Size = case_when(
        abs(Frequency_Change) >= 0.5 ~ "Large",
        abs(Frequency_Change) >= 0.2 ~ "Moderate",
        abs(Frequency_Change) >= 0.1 ~ "Small",
        TRUE ~ "Very Small"
      )
    )
  
  df_sig_magnitude <- statistical_results %>%
    filter(!is.na(Frequency_Change), abs(Frequency_Change) > delta_thresh) %>%
    arrange(Gene, desc(abs(Frequency_Change)), POS)
  
  df_sig_statistical <- statistical_results %>%
    filter(!is.na(p_adj), p_adj < alpha) %>%
    arrange(p_adj, desc(abs(Frequency_Change)))
  
  df_sig_both <- statistical_results %>%
    filter(!is.na(Frequency_Change), abs(Frequency_Change) > delta_thresh,
           !is.na(p_adj), p_adj < alpha) %>%
    arrange(p_adj, desc(abs(Frequency_Change)))
  
  if (nrow(df_sig_magnitude) > 0)
    write_csv(df_sig_magnitude, file.path(output_dir_D, paste0("Magnitude_Significant_Changes_", new_iso, ".csv")))
  if (nrow(df_sig_statistical) > 0)
    write_csv(df_sig_statistical, file.path(output_dir_D, paste0("Statistical_Significant_Changes_", new_iso, ".csv")))
  if (nrow(df_sig_both) > 0)
    write_csv(df_sig_both, file.path(output_dir_D, paste0("Both_Significant_Changes_", new_iso, ".csv")))
  
  write_csv(statistical_results, file.path(output_dir_D, paste0("Complete_Statistical_Results_", new_iso, ".csv")))
  
  message("Summary for ", new_iso, ": Total=", nrow(statistical_results),
          ", |Δ|>", delta_thresh, "=", nrow(df_sig_magnitude),
          ", padj<", alpha, "=", nrow(df_sig_statistical),
          ", both=", nrow(df_sig_both))
  
  list(magnitude_sig = df_sig_magnitude, complete = statistical_results)
}

for (orig in names(iso_map)) {
  res <- analyze_isolate_delta(orig, iso_map[[orig]])
  if (!is.null(res)) {
    all_sig_list[[iso_map[[orig]]]]   <- res$magnitude_sig
    all_stats_list[[iso_map[[orig]]]] <- res$complete
  }
}

if (length(all_sig_list) > 0) {
  all_sig      <- bind_rows(all_sig_list)
  all_complete <- bind_rows(all_stats_list)
  
  write_csv(all_sig,      file.path(output_dir_D, "Combined_Magnitude_Significant_Changes_Project1.csv"))
  write_csv(all_complete, file.path(output_dir_D, "Combined_Complete_Statistical_Results_Project1.csv"))
  
  summary_stats <- all_complete %>%
    group_by(Isolate, Gene) %>%
    summarise(
      n_positions        = n(),
      n_magnitude_sig    = sum(abs(Frequency_Change) > delta_thresh, na.rm = TRUE),
      n_statistical_sig  = sum(p_adj < alpha, na.rm = TRUE),
      n_both_sig         = sum(abs(Frequency_Change) > delta_thresh & p_adj < alpha, na.rm = TRUE),
      mean_abs_change    = mean(abs(Frequency_Change), na.rm = TRUE),
      .groups = "drop"
    )
  
  write_csv(summary_stats, file.path(output_dir_D, "Summary_Statistics_By_Gene_Project1.csv"))
} else {
  message("No significant changes found across isolates (|Δfreq| > ", delta_thresh, ").")
}

message("\n✅ Section (D) complete! Results saved to: ", output_dir_D)






# ============================================================
# Project 1: Change-in-Frequency Map (Iso1–Iso7) with x-jitter
# Input:  REPRO_Project1_2026/Figures2N3_Project1.xlsx (sheet: ChangeInFreq)
# Output: REPRO_Project1_2026/Results_ChangeInFreq_Project1/
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# ---- Project 1 roots ----
project1_root <- "/Users/deepachaudhary/Documents/aminoacidanalysis/bam/iVaroutput/parameterchange/M1try/Variants"
repro_root    <- file.path(project1_root, "REPRO_Project1_2026")

# ---- INPUT Excel (Project 1) ----
file_path <- file.path(repro_root, "Figures2N3_Project1.xlsx")
sheet_name <- "ChangeInFreq"

if (!file.exists(file_path)) {
  stop("Excel file not found: ", file_path)
}

# ---- Read sheet ----
change_freq <- read_excel(file_path, sheet = sheet_name)

# ---- Project 1 isolates: Iso1..Iso7 ----
all_isolates <- paste0("Iso", 1:7)

# ---- Long format (expects columns named Iso1..Iso7 with POS values) ----
change_freq_long <- change_freq %>%
  pivot_longer(cols = everything(), names_to = "Isolate", values_to = "POS") %>%
  filter(!is.na(POS)) %>%
  mutate(
    Isolate = trimws(Isolate),
    POS     = suppressWarnings(as.numeric(POS)),
    Isolate = factor(Isolate, levels = all_isolates)
  ) %>%
  filter(!is.na(POS), Isolate %in% all_isolates)

# ---- Genomic regions (NDV AF077761; gene blocks) ----
region_boundaries <- data.frame(
  Region = c("NP", "P", "M", "F", "HN", "L"),
  Start  = c(122, 1887, 3290, 4544, 6412, 8381),
  End    = c(1591, 3074, 4384, 6205, 8145, 14995)
)

# ---- Colors (Project 1: Iso1..Iso7) ----
region_colors <- c(
  NP = "lightblue",
  P  = "lightgreen",
  M  = "lightpink",
  F  = "lightyellow",
  HN = "lightcyan",
  L  = "lightgrey"
)

isolate_colors <- c(
  Iso1 = "blue",
  Iso2 = "red",
  Iso3 = "green",
  Iso4 = "orange",
  Iso5 = "purple",
  Iso6 = "brown",
  Iso7 = "pink"
)

# ---- Plot (x-jitter only to reduce overlaps) ----
p <- ggplot() +
  geom_rect(
    data = region_boundaries,
    aes(xmin = Start, xmax = End, ymin = -Inf, ymax = Inf, fill = Region),
    alpha = 0.30
  ) +
  geom_point(
    data = change_freq_long,
    aes(x = POS, y = Isolate, color = Isolate),
    size = 2,
    position = position_jitter(width = 20, height = 0)
  ) +
  scale_fill_manual(values = region_colors) +
  scale_color_manual(values = isolate_colors, breaks = all_isolates, limits = all_isolates) +
  scale_y_discrete(limits = all_isolates) +
  labs(
    x = "Genomic Position",
    y = "Isolate",
    fill = "Genomic Region",
    color = "Isolate",
    title = "Project 1: Change-in-Frequency SNP Map (Iso1–Iso7)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

# ---- Save outputs in REPRO folder ----
out_dir <- file.path(repro_root, "Results_ChangeInFreq_Project1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = file.path(out_dir, "Figure_ChangeInFreq_Map_Project1_Iso1toIso7.jpeg"),
  plot = p, width = 10, height = 6, dpi = 300
)

write.csv(change_freq_long,
          file.path(out_dir, "ChangeInFreq_Long_Project1.csv"),
          row.names = FALSE)

# ---- Print plot to RStudio viewer ----
print(p)

cat("\n✅ Saved Change-in-Frequency outputs to:\n  ", out_dir, "\n", sep = "")






