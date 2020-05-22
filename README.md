# MicroGMT
Mutation tracker for microbial genomes

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
* The main function of MicroGMT contains two steps: sequence_to_vcf.py, which aligns the input file to the reference genome and identify variants; and annotate_vcf.py, which annotates the variants and output summary tables.<br>
* MicroGMT also provide additional utility scripts:
&#160;1. combine_summary_tables.py: combine summary tables from different MicroGMT runs.
&#160;2. remove_from_summary_tables.py: remove unwanted strains/IDs from the summary table.
&#160;3. analysis_utilities.py: format the summary table for easy access with [R](https://www.r-project.org/), and find unique mutations.
&#160;4. Find_new_seqs.sh: find new strains/IDs from a fasta formatted file of assembled sequences that are not already in the summary tables.
&#160;5. Find_regiosn_for_new_seqs.sh: find new strains/IDs from the region file accompanying the fasta formatted file of assembled sequences that are not already in the summary tables.

## Arguments
### Step 1: sequence_to_vcf.py, which aligns the input file to the reference genome and identify variants
``` bash
usage: sequence_to_vcf.py [-h] -r REF_GENOME -i {assembly,contig,fastq} -o
                          OUT_DIR [-fs FASTA_SEQS] [-fq1 FASTQ1] [-fq2 FASTQ2]
                          [-fq FASTQ] [-l LOG] [-n NAME] [-gatk PATH_TO_GATK]
                          [-picard PATH_TO_PICARD] [-kb] [-ki] [-p PRIOR]
                          [-m MBQ] [-a {asm5,asm10,asm20}] [-t THREAD]
```
#### Optional arguments:

|   Parameter                    |     Default value     |    Explanation                             | Restrictions |
| :----------------------------: | :-------------------: | :----------------------------------------- | :----------- |
| -h, --help | - | show this help message and exit | - |

#### Mandatory arguments:

|   Parameter                    |     Default value     |    Explanation                             | Restrictions |
| :----------------------------: | :-------------------: | :----------------------------------------- | :----------- |
| -r REF_GENOME | - | Fasta formatted reference genome file  | - |
| -i {assembly,contig,fastq} | - | Type of input file | - |
| -o OUT_DIR | - | Output directory | - |

#### Additional arguments for inputs:
|   Parameter                    |     Default value     |    Explanation                             | Restrictions |
| :----------------------------: | :-------------------: | :----------------------------------------- | :----------- | 
| -fs FASTA_SEQS | - | Fasta format assembly or contig file.  | Can't be used together with -fq1, fq2 or -fq. |
| -fq1 FASTQ1 | - | Fastq file 1. For paired-end fastq data.  | Can't be used together with -fs or -fq. Must be used with -fq2. |
| -fq2 FASTQ2 | - | Fastq file 2. For paired-end fastq data.  | Can't be used together with -fs or -fq. Must be used with -fq1. |
| -fq FASTQ | - | Fastq file. For single-end fastq data.  | Can't be used together with -fq1, fq2 or -fs. |

#### Arguments for simulating rearranged genomes:

|   Parameter                    |     Default value     |    Explanation                             | Restrictions |
| :----------------------------: | :-------------------: | :----------------------------------------- | :----------- |
| -l LOG | Sequence_to_vcf.log | Name of the log file | - |
| -n NAME | test | Name of the input sample. Does not work with 'assembly' option. | - |
| -gatk PATH_TO_GATK | - | Absolute path to GenomeAnalysisTK.jar. | Only required for 'fastq' option. |
| -picard PATH_TO_PICARD | - | Absolute path to picard.jar. | Only required for 'fastq' option. |
| -kb | - | Keep BAM files. | - |
| -ki | - | Keep index files. | Only works with 'fastq' option. |
| -p PRIOR | 0 | Prior for bcftools variant caller (expected substitution rate). 0 means the prior is disabled. | Only works for 'assembly' or 'contig' option. |
| -m MBQ | 10 | Minimum base quality for variant caller. | Only works with 'fastq' option. |
| -a {asm5,asm10,asm20} | asm5 | Sequence divergence: asm5/asm10/asm20 for ~0.1/1/5 percentages. | Only works with 'assembly' option. |
| -t THREAD | 10 | Number of threads. | - |

### Step 2: annotate_vcf.py, which annotates the variants and output summary tables
``` bash
usage: annotate_vcf.py [-h] -i IN_DIR -o OUT_DIR [-r REFERENCE]
                       [-p TABLE_PREFIX] [-c] [-l LOG] [-f {a,b,both}]
                       [-rg REGION_FILE] [-na] [-eff PATH_TO_SNPEFF]

Vcf file annotation

optional arguments:
  -h, --help           show this help message and exit

Mandatory inputs:
  -i IN_DIR            Input directory
  -o OUT_DIR           Output directory

Optional arguments:
  -r REFERENCE         Name of reference. [NC_045512]
  -p TABLE_PREFIX      Prefix of summary tables for annotated vcf files. Do not include path, except for folder name(s) inside output directory!
  -c                   Create CSV format snpEff summary files.
  -l LOG               Name of the log file [Annotate_vcf.log]
  -f {a,b,both}        Format of summary tables [a]
  -rg REGION_FILE      Name of the region file (Optional)
  -na                  Skip vcf annotation step, just make summary tables from annotated vcfs.
  -eff PATH_TO_SNPEFF  Absolute path to snpEff.jar. Required if annotae vcf files.



## Questions and bug report
Please direct all questions and bug reports to Yue Xing at: yue.july.xing@gmail.com

