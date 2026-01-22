
# Genomic Variants Analysis of LoNDVs from Wild Birds Using iVar

## Project Description

This project aims to analyze genomic variants (SNPs and indels) in low pathogenic Newcastle Disease Virus (LoNDV) isolates from wild birds after passaging in chicken embryos using the iVar software. The study focuses on identifying genetic changes that occur during viral adaptation in chicken embryos.

---

## Software and Tools Used

- FastQC (v0.11.9) – Quality control of sequencing data  
- Trimmomatic (v0.39) – Adapter and quality trimming  
- BWA (v0.7.12) – Read alignment  
- Samtools (v1.13) – BAM file processing  
- iVar (v1.4.3) – Variant calling and filtering  
- R – Data processing and visualization  

---

## Pipeline Overview

1. **Quality Control (QC) of Raw Reads**
   - Performed using FastQC

2. **Trimming**
   - Adapter trimming using Trimmomatic
   - Alignment to reference genome (LaSota AF077761)

3. **Alignment**
   - Reads aligned using Bowtie2 / Hisat2 / BWA
   - Converted to BAM format using Samtools
   - Sorted and indexed BAM files

4. **Variant Calling**
   - Variant calling performed using iVar
   - Used `ivar variants` for SNP/Indel calling
   - Filtering parameters:
     - Minimum base quality: 20
     - Minimum allele frequency: 0.03

5. **Data Visualization and Interpretation (R)**
   - Total shared SNPs between P1 and P10, unique P1 and P10 SNPs
   - Genomic distribution plots of SNPs
   - Total shared and unique indels between P1 and P10
   - Genomic distribution plots of indels
   - Total shared and unique nSNPs between P1 and P10
   - Change in frequency of nSNPs from P1 → P10

---

## **Repository Structure**

```text
Variant_calling_IVAR/
...
├── Project1_variant_calling_github/
│   ├── Rcode_analysis/
│   │   ├── HNcodon_495_analysis.R
│   │   ├── Project1analysis.R
│   │   └── Superscript_meandensities_SNPs_nSNPs.R
│   ├── ncbi.sh
│   └── updatetrim.sh
├── R_script/
├── Script/
└── README.md
```

---

### **Folder Description**

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

## **Outputs**

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

##  **Reference Genome**

- Newcastle Disease Virus LaSota strain  
  **Accession:** AF077761

### **References & Documentation**
For more details on iVar commands and usage, refer to the official iVar documentation:
(https://andersen-lab.github.io/ivar/html/manualpage.html)
