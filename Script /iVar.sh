#!/bin/bash

# Load modules bwa, samtools
source /apps/profiles/modules_asax.sh.dyn
module load anaconda/3-2024.02
module load bwa/0.7.12
module load samtools/1.13
#module load ivar/1.4.3

# Indexing reference
bwa index -p ref /home/shared/hauck_research/Deepa_NDV_updated/troubleshooting/ref/AF077761/AF077761.fasta

# Perform alignment for each pair of reads
bwa mem -t 4 ref "/home/aubdxc001/hauck_research/Deepa_NDV_updated/NDV23trim/f_paired_iso8p10_S48.fq.gz" "/home/aubdxc001/hauck_research/Deepa_NDV_updated/NDV23trim/r_paired_iso8p10_S48.fq.gz" > iso8p10.sam

# Convert SAM to BAM
samtools view -Sb iso8p10.sam > iso8p10.bam

# Sort BAM file
samtools sort iso8p10.bam -o iso8p10.sorted.bam

# Index the sorted BAM file
samtools index iso8p10.sorted.bam

# Generate mpileup and call variants using ivar
samtools mpileup -aa -A -d 600000 -B -Q 0 iso8p10.sorted.bam | ivar variants -p iso8p10 -q 20 -t 0.03 -r /home/shared/hauck_research/Deepa_NDV_updated/troubleshooting/ref/AF077761/AF077761.fasta -g /home/shared/hauck_research/Deepa_NDV_updated/troubleshooting/ref/AF077761/sequence.gff3
