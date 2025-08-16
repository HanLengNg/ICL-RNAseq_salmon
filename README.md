# ICL-RNAseq_salmon
A bash script to perform Salmon on raw RNAseq fastq files

# Main information
The script is written to run on a PBS cluster that will submit multiple jobs in parallel for each sample

Create a conda environment that contain the following packages:
- fastqc/0.11.9
- multiqc/1.12
- trim-galore/0.6.7
- cutadapt/3.7
- salmon 1.10.1

Before running salmon quant, please generate a salmon index (as outlined below).

### Salmon index steps
This script is used to generate a salmon index with version 1.10
For more information on a decoy, refer to the salmon [documentation](https://salmon.readthedocs.io/en/latest/salmon.html).

To generate a decoy for salmon index, the transcript and genome assembly can be downloaded from: www.gencodegenes.org/

I followed the script outlined here to generate a decoy: https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/

```
grep "^>" <(gunzip -c GRCh38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoy_hg38.txt
sed -i.bak -e 's/>//g' decoy_h38.txt
cat gencode.v43.transcripts.fa.gz GRCh38.primary_assembly.genome.fa.gz > gentrome_hg38.fa.gz

salmon index -t gentrome_hg38.fa.gz \
  -d decoy_hg38.txt \
  -p 8 -i gencode_hg38_v43_salmonIndex --gencode
```
## Salmon quant
Please create subdirectories to store each sample paired-end Illumina fastq files.

Ensure that the name of the subdirectories are the same name as the Illumina fastq files, i.e. R001/R001_1.fq.gz & R001/R001_2.fq.gz

All fastqc, multiqc and salmon quant analyses are created in the subdirectories where the Illumina fastq files are kept.

Following salmon quant, DEseq2 can be subsequently performed (See ICL-DESeq2_script on my github page). Refer to the documentation [here](https://github.com/HanLengNg/ICL-DESeq2_script) for DESeq2 analysis.
