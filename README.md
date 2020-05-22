# MicroGMT: A mutation tracker for SARS-CoV-2 and other microbial genome sequences

## Description
MicroGMT is a python based package, which takes either raw sequence reads or assembled genome sequence as input and compares against database sequences to identify and characterize small indels and point mutations in the microbial genoems. Although our default setting is optimized for SARS-CoV-2 virus, the package can be also applied to any other microbial genomes. <br>

## Installation
No installation required. Simply download the repository and unpack by "tar".

## Requirements
* [Python 3](https://www.python.org/). Required python packages: argparse, os, subprocess, sys, time.
* [snpEff](http://snpeff.sourceforge.net/)
* [SAMtools](http://samtools.sourceforge.net/)
* [JAVA](https://www.java.com/en/)

* If the inputs are fasta formatted assembled genomes/contigs, you will also need:<br>
&#160;1. [minimap2](https://github.com/lh3/minimap2)<br>
&#160;2. [BCFtools](https://samtools.github.io/bcftools/)

* If the inputs are fasta formatted assembled genomes/contigs, you will also need:<br>
&#160;1. [GATK/3.8-1-0-Java-1.8.0](https://gatk.broadinstitute.org/hc/en-us)
&#160;2. [picard/2.18.27-Java-1.8](https://broadinstitute.github.io/picard/)
&#160;3. [BWA/0.7.17-intel-2018b](http://bio-bwa.sourceforge.net/)

## Usage
* The main function of MicroGMT contains two steps: sequence_to_vcf.py, which aligns the input file to the reference genome and identify variants; and annotate_vcf.py, which annotates the variants and output summary tables.
* MicroGMT also provide additional utility scripts:<br>
&#160;1. combine_summary_tables.py: combine summary tables from different MicroGMT runs.<br>
&#160;2. remove_from_summary_tables.py: remove unwanted strains/IDs from the summary table.<br>
&#160;3. analysis_utilities.py: format the summary table for easy access with [R](https://www.r-project.org/), and find unique mutations.<br>
&#160;4. Find_new_seqs.sh: find new strains/IDs from a fasta formatted file of assembled sequences that are not already in the summary tables.<br>
&#160;5. Find_regiosn_for_new_seqs.sh: find new strains/IDs from the region file accompanying the fasta formatted file of assembled sequences that are not already in the summary tables.

## Inputs

## Outputs

## Tutorial
### 1. Workflow for SARS-CoV-2 sequences
#### Inputs
The input file is a fasta formatted database file containing 29896 SARS-CoV-2 sequences downloaded from [GISAID](https://www.gisaid.org/) on May 20, 2020. It is named "sequences0520.fasta".<br>
This file used strain ID as the fasta header. However, since the strain IDs contain "/", we need to substitute them with something else ("_" in our example):
```bash
sed -i 's/\//_/g' sequences0520.fasta
```
We also need to substitube "/" in the region file, and substitute blanks (" ") in the region file:
```bash
awk -F "\t" '{print $1"\t"$6}' metadata0520.tsv > metadata0520.short.tsv
sed -i 's/\//_/g' metadata0520.short.tsv
sed -i 's/ /_/g' metadata0520.short.tsv
```
  
#### 





## Arguments
### sequence_to_vcf.py
``` bash
usage: sequence_to_vcf.py [-h] -r REF_GENOME -i {assembly,contig,fastq} -o
                          OUT_DIR [-fs FASTA_SEQS] [-fq1 FASTQ1] [-fq2 FASTQ2]
                          [-fq FASTQ] [-l LOG] [-n NAME] [-gatk PATH_TO_GATK]
                          [-picard PATH_TO_PICARD] [-kb] [-ki] [-p PRIOR]
                          [-m MBQ] [-a {asm5,asm10,asm20}] [-t THREAD]

optional arguments:
  -h, --help            show this help message and exit

Mandatory inputs:
  -r REF_GENOME         Fasta formatted reference genome file
  -i {assembly,contig,fastq}
                        Type of input file
  -o OUT_DIR            Output directory

Additional arguments for inputs:
  -fs FASTA_SEQS        Fasta format assembly or contig file.
  -fq1 FASTQ1           Fastq file 1. For paired-end fastq data.
  -fq2 FASTQ2           Fastq file 2. For paired-end fastq data.
  -fq FASTQ             Fastq file. For single-end fastq data.

Optional arguments:
  -l LOG                Name of the log file [Sequence_to_vcf.log]
  -n NAME               Name of the input sample. Does not work with 'assembly' option. [test]
  -gatk PATH_TO_GATK    Absolute path to GenomeAnalysisTK.jar. Only required for 'fastq' option.
  -picard PATH_TO_PICARD
                        Absolute path to picard.jar. Only required for 'fastq' option.
  -kb                   Keep BAM files.
  -ki                   Keep index files. Only works with 'fastq' option.
  -p PRIOR              Prior for bcftools variant caller (expected substitution rate). 0 means the prior is disabled. 
                        Only works for 'assembly' or 'contig' option. [0]'.
  -m MBQ                Minimum base quality for variant caller. Only works with 'fastq' option. [10]
  -a {asm5,asm10,asm20}
                        Sequence divergence: asm5/asm10/asm20 for ~0.1/1/5 percentages. 
                        Only works with 'assembly' option. [asm5]
  -t THREAD             Number of threads. [10]
```

### annotate_vcf.py
``` bash
usage: annotate_vcf.py [-h] -i IN_DIR -o OUT_DIR [-r REFERENCE]
                       [-p TABLE_PREFIX] [-c] [-l LOG] [-f {a,b,both}]
                       [-rg REGION_FILE] [-na] [-eff PATH_TO_SNPEFF]

optional arguments:
  -h, --help           show this help message and exit

Mandatory inputs:
  -i IN_DIR            Input directory
  -o OUT_DIR           Output directory

Optional arguments:
  -r REFERENCE         Name of reference. [NC_045512]
  -p TABLE_PREFIX      Prefix of summary tables for annotated vcf files. 
                       Do not include path, except for folder name(s) inside output directory!
  -c                   Create CSV format snpEff summary files.
  -l LOG               Name of the log file [Annotate_vcf.log]
  -f {a,b,both}        Format of summary tables [a]
  -rg REGION_FILE      Name of the region file (Optional)
  -na                  Skip vcf annotation step, just make summary tables from annotated vcfs.
  -eff PATH_TO_SNPEFF  Absolute path to snpEff.jar. Required if annotae vcf files.
```

### combine_summary_tables.py
``` bash
usage: combine_summary_tables.py [-h] -i1 IN_TABLE1 -i2 IN_TABLE2 -d OUT_DIR
                                 -f {a,b} [-p OUT_PREF] [-l LOG]

optional arguments:
  -h, --help     show this help message and exit

Mandatory inputs:
  -i1 IN_TABLE1  Input summary table 1
  -i2 IN_TABLE2  Input summary table 2
  -d OUT_DIR     Output directory
  -f {a,b}       Format of summary tables

Optional arguments:
  -p OUT_PREF    Prefix of the output summary tables.
                 Do not include path, except for folder name(s) inside output directory!
  -l LOG         Name of the log file [Combine_summary_tables.log]
```

### remove_from_summary_tables.py
``` bash
usage: remove_from_summary_tables.py [-h] -i IN_TABLE -r REMOVE_LIST -d
                                     OUT_DIR -f {a,b} [-p OUT_PREF] [-l LOG]

optional arguments:
  -h, --help      show this help message and exit

Mandatory inputs:
  -i IN_TABLE     Input summary table
  -r REMOVE_LIST  List of IDs to be removed from table
  -d OUT_DIR      Output directory
  -f {a,b}        Format of summary tables

Optional arguments:
  -p OUT_PREF     Prefix of the output summary tables.
                  Do not include path, except for folder name(s) inside output directory!
  -l LOG          Name of the log file [Remove_from_summary_tables.log]
```

## Questions and bug report
Please direct all questions and bug reports to Yue Xing at: yue.july.xing@gmail.com

