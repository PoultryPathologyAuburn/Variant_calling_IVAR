#!/bin/bash

# Load modules
source /apps/profiles/modules_asax.sh.dyn
module load anaconda/3-2024.02
module load bwa/0.7.12
module load samtools/1.13


# Define directories
BASE_DIR="/home/aubdxc001/hauck_research/Deepa_NDV_updated/M1retry"
REF_DIR="${BASE_DIR}/Ref"
TRIM_DIR="${BASE_DIR}/NDVtrim"
BAM_DIR="${BASE_DIR}/BAMs"
TSV_DIR="${BASE_DIR}/Variants"
REF_FASTA="${REF_DIR}/LaSota.fasta"
REF_GFF="${REF_DIR}/LaSota.gff3"
REF_INDEX="${REF_DIR}/ref_index"

# Create output folders if they don't exist
mkdir -p "$BAM_DIR" "$TSV_DIR"

# Index reference (run only once)
bwa index -p "$REF_INDEX" "$REF_FASTA"
samtools faidx "$REF_FASTA"

# Define sample names and corresponding S numbers
declare -A SAMPLE_NUMS
SAMPLE_NUMS["iso1p1"]="33"; SAMPLE_NUMS["iso1p10"]="42"
SAMPLE_NUMS["iso2p1"]="34"; SAMPLE_NUMS["iso2p10"]="43"
SAMPLE_NUMS["iso3p1"]="35"; SAMPLE_NUMS["iso3p10"]="44"
SAMPLE_NUMS["iso5p1"]="36"; SAMPLE_NUMS["iso5p10"]="45"
SAMPLE_NUMS["iso6p1"]="37"; SAMPLE_NUMS["iso6p10"]="46"
SAMPLE_NUMS["iso7p1"]="38"; SAMPLE_NUMS["iso7p10"]="47"
SAMPLE_NUMS["iso8p1"]="39"; SAMPLE_NUMS["iso8p10"]="48"
SAMPLE_NUMS["iso9p1"]="40"; SAMPLE_NUMS["iso9p10"]="49"
SAMPLE_NUMS["iso10p1"]="41"; SAMPLE_NUMS["iso10p10"]="50"

# Process each isolate and passage
for SAMPLE in "${!SAMPLE_NUMS[@]}"; do
    SNUM="${SAMPLE_NUMS[$SAMPLE]}"
    FWD="${TRIM_DIR}/f_paired_${SAMPLE}_S${SNUM}.fq.gz"
    REV="${TRIM_DIR}/r_paired_${SAMPLE}_S${SNUM}.fq.gz"

    echo "🔹 Processing $SAMPLE..."

    # Alignment
    bwa mem -t 4 "$REF_INDEX" "$FWD" "$REV" > "${SAMPLE}.sam"

    # Convert and sort BAM
    samtools view -Sb "${SAMPLE}.sam" | samtools sort -o "${BAM_DIR}/${SAMPLE}.sorted.bam"

    # Index BAM
    samtools index "${BAM_DIR}/${SAMPLE}.sorted.bam"

    # Variant calling
    samtools mpileup -aa -A -d 600000 -B -Q 0 "${BAM_DIR}/${SAMPLE}.sorted.bam" \
    | ivar variants -p "${TSV_DIR}/${SAMPLE}" -q 20 -t 0.03 -r "$REF_FASTA" -g "$REF_GFF"

    echo "✅ $SAMPLE completed."
done
