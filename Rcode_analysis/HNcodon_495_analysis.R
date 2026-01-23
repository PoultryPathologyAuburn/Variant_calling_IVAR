############################################################
# OUTPUT DESCRIPTION (Project 1 – Codon/AA validation)
#
# This script extracts codon- and amino-acid–level frequencies
# at NDV genome positions 7894–7896 directly from BAM files
# using Q30-filtered reads.
#
# Outputs written to:
#   REPRO_Project1_2026/CodonAA_7894_7896_Q30/
#
# For EACH isolate, the following files are produced:
#
# 1) *_aa_GE3pct.csv   <-- PRIMARY RESULT
#    - Amino acid frequencies ≥3%
#    - Matches iVar variant-calling threshold
#    - Used for validation, tables, and interpretation
#
# 2) *_codon_GE3pct.csv
#    - Codon-level frequencies ≥3%
#    - Used to confirm synonymous vs non-synonymous changes
#
# 3) *_aa_1to3pct.csv
#    - Amino acids between 1–3%
#    - Supplementary / emerging minority variants
#
# 4) *_aa_ALL.csv
#    - All observed amino acids (no frequency filtering)
#    - Used only for sanity checks or IGV troubleshooting
#
# 5) *_codon_ALL.csv
#    - All observed codons
#    - Not used for downstream analysis
#
# Interpretation:
# - Focus primarily on *_aa_GE3pct.csv
# - Compare results directly with iVar TSV variant calls
# - Variants <3% are considered supportive but not primary
#
############################################################

suppressPackageStartupMessages({
  library(Rsamtools)
  library(GenomicAlignments)
  library(Biostrings)
  library(GenomicRanges)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})

# ============================================================
# PROJECT 1 PATHS
# ============================================================
project1_root <- "/Users/deepachaudhary/Documents/aminoacidanalysis/bam/iVaroutput/parameterchange/M1try/Variants"


repro_root <- file.path(project1_root, "REPRO_Project1_2026")

# REF used for variant calling 
ref_fa <- file.path(repro_root, "Ref", "LaSota.fasta")

# BAM directory (Project 1)
aligned_dir <- "/Users/deepachaudhary/Documents/aminoacidanalysis/bam/iVaroutput/parameterchange/M1try/WIld_birds_Project1_BAMs"

# Main output folder (inside REPRO)
out_dir <- file.path(repro_root, "CodonAA_7894_7896_Q30")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# SETTINGS
# ============================================================
codon_pos      <- c(7894L, 7895L, 7896L)  # LaSota coords
min_baseq      <- 30L                     # Q30
ivar_threshold <- 3                       # %  (main)
supp_min       <- 1                       # %  (supplementary)
supp_max       <- 3                       # %  (supplementary upper bound, exclusive)

# ============================================================
# OUTPUT ORGANIZATION (folders)
# ============================================================
main_dir  <- file.path(out_dir, "MAIN_RESULTS")           
supp_dir  <- file.path(out_dir, "SUPPLEMENTARY_1to3pct")  # often empty
all_dir   <- file.path(out_dir, "ALL_READS_QC")           # QC only
codon_dir <- file.path(out_dir, "CODON_SUPPORT")          # codon evidence

dir.create(main_dir,  showWarnings = FALSE)
dir.create(supp_dir,  showWarnings = FALSE)
dir.create(all_dir,   showWarnings = FALSE)
dir.create(codon_dir, showWarnings = FALSE)

writeLines(
  c(
    "Codon/AA validation at nt 7894–7896 (Q30).",
    "",
    "MAIN_RESULTS:",
    "  - *_aa_GE3pct.csv  (Amino acids ≥3%; matches iVar threshold; use for tables/figures)",
    "",
    "SUPPLEMENTARY_1to3pct:",
    "  - *_aa_1to3pct.csv (Amino acids 1–3%; minority variants; often empty)",
    "",
    "ALL_READS_QC:",
    "  - *_aa_ALL.csv and *_codon_ALL.csv (all detected variants; QC/reference only)",
    "",
    "CODON_SUPPORT:",
    "  - *_codon_GE3pct.csv and *_codon_1to3pct.csv (codon-level evidence behind AA calls)"
  ),
  file.path(out_dir, "README.txt")
)

# ============================================================
# HELPERS
# ============================================================
`%||%` <- function(x, y) if (is.null(x) || is.na(x)) y else x

codon2aa <- function(codon){
  cc <- toupper(gsub("U","T", codon))
  aa <- GENETIC_CODE[as.character(DNAStringSet(cc))]
  ifelse(is.na(aa), "X", aa)
}

pick_chrom <- function(bam){
  hdr <- scanBamHeader(bam)[[1]]$targets
  chs <- names(hdr)
  cand <- chs[grepl("^AF077761(\\.|$)", chs)]
  if (length(cand)) return(cand[1])
  cand <- chs[grepl("LaSota|AF077761|NC_", chs, ignore.case=TRUE)]
  if (length(cand)) return(cand[1])
  chs[1]
}

ref_to_read_pos <- function(read_start, cigar, ref_pos_vec){
  if (is.na(cigar) || cigar=="*" || length(cigar)==0) return(rep(NA_integer_, length(ref_pos_vec)))
  cigar_chr <- as.character(cigar)[1]
  ops  <- explodeCigarOps(cigar_chr); if (is.list(ops))  ops  <- unlist(ops,  use.names=FALSE)
  lens <- explodeCigarOpLengths(cigar_chr); if (is.list(lens)) lens <- unlist(lens, use.names=FALSE)
  lens <- as.integer(lens)
  if (length(ops)!=length(lens)) return(rep(NA_integer_, length(ref_pos_vec)))
  
  ref_cur <- as.integer(read_start)[1]
  read_cur <- 1L
  out <- rep(NA_integer_, length(ref_pos_vec))
  
  for (k in seq_along(ops)){
    op <- ops[k]; len <- lens[k]
    if (is.na(len) || len<=0L) next
    
    if (op %in% c("M","=","X")){
      ref_start <- ref_cur; ref_end <- ref_cur + len - 1L
      hit <- (ref_pos_vec>=ref_start) & (ref_pos_vec<=ref_end) & is.na(out)
      if (any(hit)){
        off <- ref_pos_vec[hit] - ref_start
        out[hit] <- read_cur + off
      }
      ref_cur <- ref_cur + len
      read_cur <- read_cur + len
      
    } else if (op=="I"){
      read_cur <- read_cur + len
    } else if (op %in% c("D","N")){
      ref_cur <- ref_cur + len
    } else if (op=="S"){
      read_cur <- read_cur + len
    }
  }
  out
}

read_codon <- function(seq_str, qual_raw, read_pos_idx, min_baseq=30L){
  if (any(is.na(read_pos_idx))) return(NA_character_)
  n <- nchar(seq_str)
  if (any(read_pos_idx<1L) || any(read_pos_idx>n)) return(NA_character_)
  
  q <- tryCatch({
    if (is.raw(qual_raw)) as.integer(qual_raw)
    else if (inherits(qual_raw, "XStringQuality")) as.integer(quality(qual_raw))
    else as.integer(qual_raw)
  }, error=function(...) rep(40L, n))
  
  bases <- substring(seq_str, read_pos_idx, read_pos_idx)
  if (any(q[read_pos_idx] < min_baseq) || any(toupper(bases) %in% "N")) return(NA_character_)
  toupper(paste0(bases, collapse=""))
}

count_codons_bam <- function(bam_path, chrom, codon_pos, min_baseq=30L,
                             skip_secondary=TRUE, skip_supplementary=TRUE, skip_qc_fail=TRUE){
  
  rng <- GRanges(seqnames=chrom, ranges=IRanges(min(codon_pos), max(codon_pos)))
  param <- ScanBamParam(which=rng, what=c("pos","cigar","seq","qual","flag"))
  x <- scanBam(bam_path, param=param)[[1]]
  if (length(x$pos)==0) return(list(counts=integer(), used=0L))
  
  flag <- x$flag
  keep <- rep(TRUE, length(flag))
  if (skip_secondary)     keep <- keep & (bitwAnd(flag, 0x100)==0)
  if (skip_supplementary) keep <- keep & (bitwAnd(flag, 0x800)==0)
  if (skip_qc_fail)       keep <- keep & (bitwAnd(flag, 0x200)==0)
  
  pos   <- x$pos[keep]
  cig   <- x$cigar[keep]
  seqs  <- as.character(x$seq[keep])
  quals <- x$qual[keep]
  
  used <- 0L
  tab <- integer(); names(tab) <- character()
  
  for (i in seq_along(pos)){
    if (is.na(cig[i]) || cig[i]=="*") next
    rpos <- ref_to_read_pos(pos[i], cig[i], codon_pos)
    cod  <- read_codon(seqs[i], quals[[i]], rpos, min_baseq)
    if (is.na(cod)) next
    used <- used + 1L
    tab[cod] <- (tab[cod] %||% 0L) + 1L
  }
  list(counts=tab, used=used)
}

# ============================================================
# CHECKS
# ============================================================
if (!file.exists(ref_fa)) stop("REF fasta not found: ", ref_fa)
if (!file.exists(paste0(ref_fa, ".fai"))) stop("REF fasta index (.fai) not found: ", paste0(ref_fa, ".fai"))

bam_files <- list.files(aligned_dir, pattern = "\\.sorted\\.bam$", full.names = FALSE)
if (length(bam_files) == 0) stop("No .sorted.bam files found in: ", aligned_dir)

message("Found ", length(bam_files), " BAMs in: ", aligned_dir)
message("Output root: ", out_dir)
message("REF used: ", ref_fa)

# ============================================================
# LOOP OVER ALL BAMS
# ============================================================
summary_all <- list()

for (bam_file in bam_files) {
  
  bam_path <- file.path(aligned_dir, bam_file)
  bai_path <- paste0(bam_path, ".bai")
  
  # auto sample name from bam filename (iso10p1.sorted.bam -> iso10p1)
  sample_name <- sub("\\.sorted\\.bam$", "", bam_file)
  
  if (!file.exists(bai_path)) {
    warning("Missing BAI index, skipping: ", bam_file)
    next
  }
  
  chrom <- pick_chrom(bam_path)
  message("\n--- Running: ", sample_name, " | contig: ", chrom, " ---")
  
  res <- count_codons_bam(bam_path, chrom, codon_pos, min_baseq)
  if (res$used == 0L) {
    warning("No usable reads (Q", min_baseq, ") at ", paste(codon_pos, collapse="-"),
            " for ", sample_name, " — skipping")
    next
  }
  
  codon_counts <- res$counts
  reads_used   <- res$used
  
  codon_tbl_all <- tibble(
    Sample     = sample_name,
    Codon      = names(codon_counts),
    AA         = codon2aa(names(codon_counts)),
    Reads      = as.integer(codon_counts),
    Reads_used = reads_used
  ) %>%
    mutate(Percent = round(100 * Reads / Reads_used, 4)) %>%
    arrange(desc(Reads))
  
  aa_tbl_all <- codon_tbl_all %>%
    group_by(Sample, AA) %>%
    summarise(Reads = sum(Reads),
              Reads_used = dplyr::first(Reads_used),
              .groups="drop") %>%
    mutate(Percent = round(100 * Reads / Reads_used, 4)) %>%
    arrange(desc(Reads))
  
  codon_tbl_ge3 <- codon_tbl_all %>% filter(Percent >= ivar_threshold)
  aa_tbl_ge3 <- codon_tbl_ge3 %>%
    group_by(Sample, AA) %>%
    summarise(Reads = sum(Reads),
              Reads_used = dplyr::first(Reads_used),
              .groups="drop") %>%
    mutate(Percent = round(100 * Reads / Reads_used, 4)) %>%
    arrange(desc(Reads))
  
  codon_tbl_1to3 <- codon_tbl_all %>% filter(Percent >= supp_min, Percent < supp_max)
  aa_tbl_1to3 <- codon_tbl_1to3 %>%
    group_by(Sample, AA) %>%
    summarise(Reads = sum(Reads),
              Reads_used = dplyr::first(Reads_used),
              .groups="drop") %>%
    mutate(Percent = round(100 * Reads / Reads_used, 4)) %>%
    arrange(desc(Reads))
  
  # ---- write files into organized folders ----
  pos_tag <- paste(codon_pos, collapse = "_")
  base <- paste0(sample_name, "_", pos_tag, "_Q", min_baseq)
  
  # MAIN RESULTS (AA ≥3%)
  write_csv(aa_tbl_ge3, file.path(main_dir, paste0(base, "_aa_GE3pct.csv")))
  
  # SUPPLEMENTARY (AA 1–3%) (may be empty)
  write_csv(aa_tbl_1to3, file.path(supp_dir, paste0(base, "_aa_1to3pct.csv")))
  
  # QC (ALL)
  write_csv(aa_tbl_all,    file.path(all_dir,   paste0(base, "_aa_ALL.csv")))
  write_csv(codon_tbl_all, file.path(all_dir,   paste0(base, "_codon_ALL.csv")))
  
  # CODON evidence
  write_csv(codon_tbl_ge3,  file.path(codon_dir, paste0(base, "_codon_GE3pct.csv")))
  write_csv(codon_tbl_1to3, file.path(codon_dir, paste0(base, "_codon_1to3pct.csv")))
  
  # ---- store summary row for this sample (top AA at ≥3%) ----
  top_aa <- aa_tbl_ge3 %>% slice_max(order_by = Reads, n = 1, with_ties = FALSE)
  summary_all[[sample_name]] <- tibble(
    Sample = sample_name,
    Reads_used = reads_used,
    Top_AA_GE3pct = ifelse(nrow(top_aa) == 0, NA_character_, top_aa$AA[1]),
    Top_AA_percent_GE3pct = ifelse(nrow(top_aa) == 0, NA_real_, top_aa$Percent[1]),
    N_AA_GE3pct = nrow(aa_tbl_ge3),
    N_AA_1to3pct = nrow(aa_tbl_1to3)
  )
}

# ============================================================
# WRITE ONE SUMMARY TABLE (easy to check everything)
# ============================================================
summary_df <- bind_rows(summary_all) %>% arrange(Sample)
summary_csv <- file.path(out_dir, "SUMMARY_all_samples_AA_GE3pct.csv")
write_csv(summary_df, summary_csv)

cat("\n============================================================\n")
cat("✅ DONE!\n")
cat("Main results folder:\n  ", main_dir, "\n", sep="")
cat("Summary table:\n  ", summary_csv, "\n", sep="")
cat("README:\n  ", file.path(out_dir, "README.txt"), "\n", sep="")
cat("============================================================\n")





library(readr)
library(dplyr)

files <- list.files(file.path(out_dir, "MAIN_RESULTS"),
                    pattern = "_aa_GE3pct.csv",
                    full.names = TRUE)

all_aa <- bind_rows(lapply(files, read_csv))
all_aa
list.files(file.path(out_dir, "MAIN_RESULTS"))

files_all <- list.files(file.path(out_dir, "ALL_READS_QC"),
                        pattern = "_aa_ALL.csv",
                        full.names = TRUE)

all_aa_all <- bind_rows(lapply(files_all, read_csv))
all_aa_all


