#!/bin/bash

# Script is written as PBS array jobs

# Packages required for Salmon quant
# Install packages in a conda environment
# fastqc/0.11.9
# multiqc/1.12
# trim-galore/0.6.7
# cutadapt/3.7
# salmon/1.10.1

# First create subdirectories for each sample and store paired-end fastq files in the respective subdirectories
# Ensure that the name of the subdirectory matches the name of the fastq file

base=$(ls -d R0* | head -n $PBS_ARRAY_INDEX | tail -n 1)
cd ./${base} 

# Run fastqc and multiqc
mkdir ${base}_pretrim_fastqc
fastqc -d ./ -o ${base}_pretrim_fastqc/ ${base}_1.fq.gz ${base}_2.fq.gz
mkdir ${base}_pretrim_multiqc
multiqc -o ${base}_pretrim_multiqc/ ${base}_pretrim_fastqc/

# Run trim_galore to remove adaptors from data
trim_galore --gzip --paired ${base}_1.fq.gz ${base}_2.fq.gz
mkdir ${base}_trimmed_fastqc
fastqc -d ./ -o ${base}_trimmed_fastqc/ ${base}_1_val_1.fq.gz ${base}_2_val_2.fq.gz
mkdir ${base}_trimmed_multiqc
multiqc -o ${base}_trimmed_multiqc/ ${base}_trimmed_fastqc/

# Refer to READMe on steps to create a salmon index before running salmon quant
salmon quant -l A -i ./salmon_index/gencode_hg38_v43_salmonIndex \
        -1 ${base}_1_val_1.fq.gz -2 ${base}_2_val_2.fq.gz \
        -o new_salmon_quant -p 4 --gcBias --writeUnmappedNames
