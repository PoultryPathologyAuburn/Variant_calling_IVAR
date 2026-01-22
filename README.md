# Genomic Variants Analysis of LoNDVs from Wild Birds Using IVAR 

### Project Description
This project aims to analyze genomic variants (SNPs and indels) in low pathogenic Newcastle Disease Virus (LoNDV) isolates from wild birds by passaging in chicken embryo using the IVAR software. The study focuses on identifying genetic changes that occur during viral adaptation in chicken embryos.

### Software and Tools Used and Pipeline Overview 
- FASTQC(v/0.11.9) – Quality control of sequencing data
- Trimmomatic (v0.39)- Perform trimming
- BWA (v0.7.12)– Perform alignment 
- Samtools (v1.13)– BAM file processing
- IVAR (v1.4.3) – Variant calling and filtering
- R – Data processing and visualization

### Pipeline Overview
The following pipeline was used for variant analysis:

1. **Quality Control (QC) of Raw Reads**
   Performed using FastQC
2. **Trimming**
   Adapter trimming using Trimmomatic
   Alignment to Reference Genome (LaSota AF077761)
3. **Alignment**
   Aligned reads using Bowtie2/Hisat2/BWA
   Converted to BAM format using Samtools
   Sorted and indexed BAM files
4. **Variant Calling**
   Variant Calling with IVAR
   Used ivar variants for SNP/Indel calling
   IVAR filtering parameters: minimum quality(Default: 20) and minimum      frequency(Default: 0.03)
5. **Data Visualization and Interpretation**
   Data Visualization and Interpretation were done using R which includes:
   - Total shared SNPs between P1 and P10, unique P1 and P10 SNPs
   - Generated genomic distribution plots of SNPs between P1 and P10
   - Total shared SNPs between P1 and P10, unique P1 and P10 indels
   - Generated genomic distribution plots of indels between P1 and P10
   - Total shared nSNPs between P1 and P10, unique P1 and P10 nSNPs
   - Change in frequency in nSNPs during P1 to P10 transition


## 📁 Repository Structure
Variant_calling_IVAR/
├── Project1_variant_calling_github/
│ ├── Rcode_analysis/
│ │ ├── HNcodon_495_analysis.R
│ │ ├── Project1analysis.R
│ │ └── Superscript_meandensities_SNPs_nSNPs.R
│ ├── ncbi.sh
│ └── updatetrim.sh


### 📂 Folder Description

- `Project1_variant_calling_github/`
  - Contains the main analysis pipeline for Project 1
- `Rcode_analysis/`
  - R scripts for variant analysis, visualization, and statistics
- `HNcodon_495_analysis.R`
  - Focused analysis of HN gene codon 495 mutation
- `Project1analysis.R`
  - Main R script for SNP, indel, and nSNP analysis, mapping of SNPs and nSNPs, change in frequency and mapping of them, reccurence of (shared and unique P1 and P10) variants among all isolates 
- `Superscript_meandensities_SNPs_nSNPs.R`
  - SNP and nSNP density and distribution analysis
- `ncbi.sh`
  - to call variants using iVar tool
- `updatetrim.sh`
  - Trimming script

---

## Outputs

The pipeline generates:

- Tables of:
  - Shared vs unique SNPs (P1 vs P10)
  - Shared vs unique indels
  - Shared vs unique nSNPs
- Genomic distribution plots:
  - SNPs across genome
  - Indels across genome
  - nSNPs across genome
- Frequency change plots for nSNPs (P1 → P10)
- Codon-level mutation summaries

---

##  Reference Genome

- Newcastle Disease Virus LaSota strain  
  **Accession:** AF077761

### References & Documentation
For more details on iVar commands and usage, refer to the official iVar documentation:
(https://andersen-lab.github.io/ivar/html/manualpage.html)
