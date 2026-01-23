###########################
## Project 1: CLD with rcompanion for SNP/nSNP DENSITIES (Iso*_SD)
## - Kruskal across genes
## - Dunn post-hoc (Bonferroni)
## - rcompanion::cldList letters
## - Outputs one Excel with KW summary + per-sheet Means/Pairs/CLD

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(FSA)          # dunnTest
  library(rcompanion)   # cldList
  library(writexl)
  library(tibble)
})



####Superscript

in_xlsx  <- file.path(repro_root, "Meandensity_Re.xlsx")   # <-- put file here OR change path
out_dir  <- file.path(repro_root, "Stat_Rerun_rcompanion_Project1")
out_xlsx <- file.path(out_dir, "SNPnSNP_density_rcompanion_CLD_Project1.xlsx")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

alpha <- 0.05



sheet_map <- c(
  SharedSNPs     = "Shared_SNPs",
  UniqueP1SNPs   = "UniqueP1_SNPs",
  UniqueP10SNPs  = "UniqueP10_SNPs",
  SharednSNPs    = "Shared_nSNPs",
  UniqueP1nSNPs  = "UniqueP1_nSNPs",
  UniqueP10nSNPs = "UniqueP10_nSNPs"
)


to_long_sd <- function(df){
  need <- paste0("Iso", 1:7, "_SD")
  missing <- setdiff(need, names(df))
  if (length(missing) > 0) {
    stop("Missing Iso*_SD columns (Project 1 expects Iso1_SD..Iso7_SD): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  
  df %>%
    select(Gene, all_of(need)) %>%
    pivot_longer(cols = all_of(need), names_to = "Isolate", values_to = "Density") %>%
    mutate(
      Gene    = factor(Gene),
      Density = suppressWarnings(as.numeric(Density))
    ) %>%
    filter(!is.na(Density))
}

# ---- helper: build CLD letters from Dunn table ----
cld_from_dunn <- function(dunn_tbl, threshold = 0.05){
  tmp <- dunn_tbl %>%
    transmute(
      comparison = gsub("\\s+", "", Comparison),  # "A-B"
      p.adjusted = P.adj
    )
  
  cld <- rcompanion::cldList(p.adjusted ~ comparison, data = tmp, threshold = threshold)
  cld <- as.data.frame(cld, stringsAsFactors = FALSE)
  
  # normalize names to Gene / Letters
  if ("Group"  %in% names(cld)) names(cld)[names(cld) == "Group"]  <- "Gene"
  if ("Letter" %in% names(cld)) names(cld)[names(cld) == "Letter"] <- "Letters"
  if (!"Gene"    %in% names(cld)) names(cld)[1] <- "Gene"
  if (!"Letters" %in% names(cld)) names(cld)[2] <- "Letters"
  
  cld$Gene    <- as.character(cld$Gene)
  cld$Letters <- as.character(cld$Letters)
  cld
}

# ---- run per sheet ----
sheets_out <- list()
kw_rows    <- list()

for (sh in names(sheet_map)) {
  label <- sheet_map[[sh]]
  
  wide <- read_excel(in_xlsx, sheet = sh)
  dat  <- to_long_sd(wide)
  
  # Kruskal–Wallis
  kw <- kruskal.test(Density ~ Gene, data = dat)
  
  # Dunn (Bonferroni)
  dunn <- FSA::dunnTest(Density ~ Gene, data = dat, method = "bonferroni")$res %>%
    mutate(
      Comparison = as.character(Comparison),
      P.adj = as.numeric(P.adj),
      P.unadj = as.numeric(P.unadj),
      Z = as.numeric(Z)
    )
  
  # CLD
  cld_tbl <- cld_from_dunn(dunn, threshold = alpha) %>% arrange(Gene)
  
  # Means per gene
  means_tbl <- dat %>%
    group_by(Gene) %>%
    summarise(Mean_Density = mean(Density, na.rm = TRUE), .groups = "drop") %>%
    mutate(Gene = as.character(Gene)) %>%
    arrange(Gene)
  
  # Pairwise compact table
  pairs_tbl <- dunn %>%
    transmute(
      Comparison = gsub("\\s+", "", Comparison),
      Z          = Z,
      P.unadj    = P.unadj,
      P.adj      = P.adj,
      Significant= ifelse(P.adj <= alpha, "yes", "no")
    )
  
  # Means + letters
  summary_tbl <- means_tbl %>%
    left_join(cld_tbl, by = "Gene") %>%
    mutate(
      Mean_with_letters =
        ifelse(is.na(Letters) | Letters == "",
               sprintf("%.4f", Mean_Density),
               paste0(sprintf("%.4f", Mean_Density), "^", Letters))
    )
  
  sheets_out[[paste0(label, "_Means")]] <- summary_tbl
  sheets_out[[paste0(label, "_Pairs")]] <- pairs_tbl
  sheets_out[[paste0(label, "_CLD")]]   <- cld_tbl
  
  kw_rows[[length(kw_rows) + 1]] <- tibble(
    Sheet     = label,
    KW_stat   = as.numeric(kw$statistic),
    KW_df     = as.numeric(kw$parameter),
    KW_pvalue = as.numeric(kw$p.value),
    Posthoc   = "Dunn (Bonferroni)",
    CLD_tool  = "rcompanion::cldList",
    Isolates  = "Iso1_SD..Iso7_SD"
  )
}

KW_results <- bind_rows(kw_rows)
sheets_out <- c(list(KW_results = KW_results), sheets_out)

write_xlsx(sheets_out, path = out_xlsx)
cat("✅ Wrote Project 1 rcompanion CLD results to:\n", out_xlsx, "\n", sep = "")
