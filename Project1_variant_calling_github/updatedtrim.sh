#!/bin/bash
source /apps/profiles/modules_asax.sh.dyn
module load trimmomatic/0.39

#######################
# SLURM for 36 samples: 2 CPU, 10min, 200mb
#######################

# Define the directory path
directory="/home/shared/hauck_research/Deepa_NDV_updated/RawfileNDV"

# Enter the directory
cd "$directory"

# Empty variable for the list
sample_list=()

# Loop through files matching the pattern
for file in ./*_L001_R1_001.fastq.gz; do
    # Check if any files match the pattern
    if [ -e "$file" ]; then
        # Get the base name without _L001_R1_001.fastq.gz
        base_name=$(basename "$file" _L001_R1_001.fastq.gz)
        # Add the base name to the sample_list array
        sample_list+=("$base_name")
    fi
done

# Loop through the list with the samples and perform Trimmomatic
for base_name in "${sample_list[@]}"; do
        # Check
        echo "CURRENT FILE: $base_name"
        java -jar /apps/x86-64/apps/spack_0.19.1/spack/opt/spack/linux-rocky8-zen3/gcc-11.3.0/trimmomatic-0.39-iu723m2xenra563gozbob6ansjnxmnfp/bin/trimmomatic-0.39.jar \
        PE -threads 6 -phred33 \
        ${base_name}_L001_R1_001.fastq.gz ${base_name}_L001_R2_001.fastq.gz \
        f_paired_${base_name}.fq.gz f_unpaired_${base_name}.fq.gz \
        r_paired_${base_name}.fq.gz r_unpaired_${base_name}.fq.gz \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
done

cd ..

# New directory for the outputs
mkdir -p NDVtrim

# Move outputs to the new directory
mv "$directory"/*.fq.gz NDVtrim
