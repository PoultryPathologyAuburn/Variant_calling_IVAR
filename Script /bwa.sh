#!/bin/bash

# load modules bwa, samtools
#source /opt/asn/etc/asn-bash-profiles-special/modules.sh
source /apps/profiles/modules_asax.sh.dyn
module load bwa/0.7.12
module load samtools/1.13
module load gatk/4.4.0.0

#indexing ref
#bwa index -p ref AF077761.fasta


# Perform alignment for each pair of reads
bwa mem -t 4 ref "/home/aubdxc001/hauck_research/Deepa_NDV_updated/NDV23trim/f_paired_iso1p10_S42.fq.gz" "/home/aubdxc001/hauck_research/Deepa_NDV_updated/NDV23trim/r_paired_iso1p10_S42.fq.gz" > iso1p10.sam

#Convert SAM to BAM:
samtools view -bS -o iso1p10.bam iso1p10.sam 
 
#Sort the BAM File
samtools sort -o iso1p10_sorted.bam iso1p10.bam

#Index the Sorted BAM File
samtools index iso1p10_sorted.bam

# Add read groups to the sorted BAM file
gatk --java-options "-Xmx4G" AddOrReplaceReadGroups \
 -I iso1p10_sorted.bam \
 -O iso1p10_rg.sorted.bam \
 -RGID 1 -RGLB 1 -RGPL ILLUMINA -RGPU 1 -RGSM iso1p10

# Index the Sorted BAM File
samtools index iso1p10_rg.sorted.bam