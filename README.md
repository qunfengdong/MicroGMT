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
``` bash

```

## Questions and bug report
Please direct all questions and bug reports to Yue Xing at: yue.july.xing@gmail.com

