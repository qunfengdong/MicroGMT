# MicroGMT: Microbial Genomics Mutation Tracker
A mutation tracker for SARS-CoV-2 and other microbial genome sequences

## Updates
### Important update for SARS-CoV-2:
To handle the multiple CDS and the -1 ribosomal frameshift in ORF1ab of SARS-CoV-2, attributes were added to the gene ID and gene name of ORF1ab to denote which CDS the mutation is on in the output vcfs and summary tables.

There are two CDS for ORF1ab gene: 1) CDS joining (266..13468,13468..21555), on which the -1 ribosomal frameshift occurs during translation, produces pp1ab; 2) CDS of (266..13483) produces pp1a.

* For mutations occur on mature peptides produced by both pp1a and pp1ab, or by pp1a only, the gene ID and name in output vcfs and summary tables are: GU280_gp01_pp1a and ORF1ab_pp1a.
* For mutations occur on mature peptides produced by pp1ab only, the gene ID and name in output vcfs and summary tables are: GU280_gp01_pp1ab and ORF1ab_pp1ab.

### Version 1.3 (June 6 2020) update
* Added a utility script to check the running environment.
* Added options to fine tuning the BWA alignment for raw reads inputs.
* Added mutation frequency summary outputs for step 2.

### Version 1.2 (May 31 2020) update
* Added 4 utility scripts: get_ids.sh, add_custom_annotation.py, mask_sequences.py and sequence_ID_extractor.py.
* Revised analysis_utilities.py to accommodate the added features of the new utility scripts.

### Version 1.1 (May 31 2020) update
To handle the multiple CDS and the -1 ribosomal frameshift in ORF1ab of SARS-CoV-2, attributes were added to the gene ID and gene name of ORF1ab to denote which CDS the mutation is on in the output vcfs and summary tables (see "Important update for SARS-CoV-2" above）.

## Description
MicroGMT is a python based package, which takes either raw sequence reads or assembled genome sequence as input and compares against database sequences to identify and characterize small indels and point mutations in the microbial genoems. Although our default setting is optimized for SARS-CoV-2 virus, the package can be also applied to any other microbial genomes.

## Citation
Xing Y, Li X, Gao X, Dong Q (2020) MicroGMT: A mutation tracker for SARS-CoV-2 and other microbial genome sequences. Accepted. Front. Microbiol. doi: 10.3389/fmicb.2020.01502. 

## Installation
No installation required. Simply download the repository and unpack by "tar".

## Requirements
* [Python 3](https://www.python.org/). Required python modules: argparse, os, subprocess, sys, time. Python module pandas is also required if users want to further output mutation frequency tables from summary tables.
* [snpEff 4.3t](http://snpeff.sourceforge.net/)
* [SAMtools 1.6 or above](http://samtools.sourceforge.net/)
* [JAVA 1.8 or above](https://www.java.com/en/)

* If the inputs are fasta genome assembly sequences/contigs, you will also need:<br>
&#160;1. [minimap2](https://github.com/lh3/minimap2)<br>
&#160;2. [BCFtools 1.6 or above](https://samtools.github.io/bcftools/)

* If the inputs are fastq formatted raw reads, you will also need:<br>
&#160;1. [GATK 3.8](https://gatk.broadinstitute.org/hc/en-us): NOT GATK 4 because GATK 3.8 has the Indel Realigner to use.<br>
&#160;2. [picard 2 or above](https://broadinstitute.github.io/picard/)<br>
&#160;3. [BWA 0.7 or above](http://bio-bwa.sourceforge.net/)

## Usage
* The main function of MicroGMT contains two steps: sequence_to_vcf.py (step 1), which aligns the input file(s) to the reference genome and identify variants; and annotate_vcf.py (step 2), which annotates the variants and output summary tables. Note: annotate_vcf.py (step 2) can take in the step 1 outputs from multiple runs, as long as they are in one input folder. So you can first process all samples by step 1, and then run them all together by step 2.
* MicroGMT also provide additional utility scripts:<br>
&#160;1. Check_environment.sh: check running environment<br>
&#160;2. combine_summary_tables.py: combine summary tables from different MicroGMT runs.<br>
&#160;3. remove_from_summary_tables.py: remove unwanted strains/sequences from the summary table.<br>
&#160;4. analysis_utilities.py: reformat the summary table for easy access with [R](https://www.r-project.org/) or other tools, or find unique mutations (unqiue mutations are defined by only one strain/sequence has that mutation at a specific locus).<br>
&#160;5. Find_new_seqs.sh: find new strains/IDs from a fasta formatted file of assembled sequences that are not already in existing summary tables.<br>
&#160;6. Find_regiosn_for_new_seqs.sh: extract region information from a region file for a list of strains/sequences.
&#160;7. get_ids.sh: extract the strain/sequence IDs from a summary table.
&#160;8. add_custom_annotation.py: add custom annotations to summary tables according to geomic coordinates of the annotation features.
&#160;9. mask_sequences.py: mask sequences according to user-supplied geomic coordinates.
&#160;10. sequence_ID_extractor.py: extract and summarize mutation information for specific strains/sequences.

## Inputs
### The main functions:
* The fasta genome reference file. For SARS-CoV-2, its fasta genome reference file is located at <path_to_MicroGMT>/NC_045512_source_files/NC_045512.fa. It is downloaded from https://www.ncbi.nlm.nih.gov/nuccore/nc_045512. The accession number ".2" is deleted from the fasta header.
* The annotation database. For SARS-CoV-2, its annotation database is pre-built and is the default database. For user-supplied genomes, please see "Pre-built annotation database for SARS-CoV-2 and build own annotation databases for user-supplied genomes" for building own databases.
* A fasta formatted database sequence file, which can contain one or multiple fasta genome assembly sequences from multiple strains (i.e. fasta genome assembly sequences downloaded from [NCBI](https://www.ncbi.nlm.nih.gov/) or [GISAID](https://www.gisaid.org/. The fasta headers should be strain/sequence IDS.). Called "fasta assembly file" below). Or a fastq formatted single/paired end raw sequence file  (called "fastq raw reads file" below). Or a fasta formatted contig sequence file from one sample (i.e. a collection of short fasta sequences from one sample, called "fasta contig file" below. **Caution: the contig sequence file option is not tested. Use at your own risk**).
* Optional: A tab delimited regional file contain region information of the samples. Format: "strain/sequence_ID	region(without blanks)"

**Note: Please don't use duplicate fasta headers and sample ID prefixes for different strains/samples. MicroGMT is not checking duplicate IDs.**

### The utility scripts:
All are optional depending on which script to use. Please see "Quick start" and "Tutorial" for more details.
* Summary tables.
* Fasta assembly file(s).
* Regional file.
* Custom annotation file.
* Sequence mask file.
* An id list containg the strain/sequence IDs in the summary tables. One strain/sequence ID per line. Needed for some utility scripts (please see "Quick start" and "Tutorial" sections). It is produced by sequence_to_vcf.py automatically for fasta formatted database sequence file inputs. Users can modify it by manually adding/deleting IDs from it, or use our utility scripts (see "Quick start" and "Tutorial" sections). For the pre-built summary tables for SARS-CoV-2, it is provided with the summary tables. It can also be produced by get_ids.sh.

## Outputs
### sequence_to_vcf.py:
* An id list containg the strain/sequence IDs processed. One strain/sequence ID per line.
* Vcf format variant calling files. One for each strain/sequence. For fasta assembly file, the strain/sequence IDs are the fasta headers (before space, if any). For fastq raw reads file, the strain/sequence ID is the prefix provided when running MicroGMT.
* Log file. For raw read sequence input, the log file contains the alignment quality information.

### annotate_vcf.py:
* Vcf formatted variant calling files with variants annotated. File names end by "anno.vcf". Can be skipped. You may produce annotated vcfs from multiple runs, put them in a same folder, then run this step with "-sa" to skip annotate vcfs, and make summary tables for all the input annotated vcfs from multiple runs.
* Tab delimited summary file produced by snpEff. File names end by "snpEff_summary.genes.txt". Need to provide file prefix by "-p".
* Csv format snpEff summary file (optional, different than "snpEff_summary.genes.txt" file). File names end by "snpEff_summary.csv".
* Tab delimited mutation summary tables of all vcf files in the input folder. Columns represent IDs of strains/sequences. Rows represent mutation loci. The summary tables have two formats: Format 1, one locus per line, each cell has the gene ID, gene name with mutation information for that locus; format 2, one locus per line with the mutated gene ID and name, each cell has the mutation information. Different summary files are provided for each format: all information ("all"), the gene ID and name the mutation locates ("gene", only for format 1), effect of the mutation ("effect"), the mutation on DNA sequence level ("gene_mut"), the gene ID and name the mutation locates along with the DNA level mutation ("gene_name_mut"only for format 1), mutation type ("mut_type"), CDS change ("cds_change"), and amino acid change ("prot_change"). In the cells, if the strain/sequence has no mutation at a specific loci, that cell is labelled by "R". If the region files is provided as input, the column headers will have both strain/sequence ID and region information, separated by "|". Please see the provided sample output summary tables in "test_dataset" folder as examples. **Note: To distinguish from the two CDS produced by ORF1ab, for SARS-CoV-2 output summary tables, if a mutation takes place in the mature peptide region produced by pp1ab only, the gene ID and name will be "GU280_gp01_pp1ab" and "ORF1ab_pp1ab"; if a mutation takes place in the mature peptide region produced by both pp1a an pp1ab, or only by pp1a, the gene ID and name will be "GU280_gp01_pp1a" and "ORF1ab_pp1a".**
* Frequency summary tables for loci and strains/sequences, which summarizes the output summary tables. Can be skipped.
* Log file.

### Utility scripts:
Optional outputs include the following. For more details about the outputs of each utility script, please see "Quick start" and "Tutorial". For sample outputs please look at "test_dataset" folder.
* Summary tables
* Reformatted summary tables and other useful tables produced based on summary tables
* Regional file and fasta assembly file with selected strains/IDs
* ID lists to extract strains/IDs from fasta assembly file and region file.
* Log file.

## The pre-built annotation database for SARS-CoV-2
The annotation database is built by snpEff. For SARS-CoV-2, the annotation database is pre-built in <path_to_MicroGMT>/database and is the default database in variant annotaion. It is built by revised NC_045512's GenBank file downloaded from https://www.ncbi.nlm.nih.gov/nuccore/nc_045512 to handle the multiple CDS and the -1 ribosomal frameshift of ORF1ab. The version number of the genome is 2. Please see below about how this database was built.

Other useful files for SARS-CoV-2 are stored in <path_to_MicroGMT>/NC_045512_source_files. These files are based on NCBI's fasta sequence and annotation files of NC_045512.

## Build own annotation databases for user-supplied genomes
For user-supplied genomes, you can find out if the genome is supported by snpEff:
```bash
java -jar <path_to_snpEff>/snpEff.jar databases
```
If supported, you can download it by:
```bash
java -jar <path_to_snpEff>/snpEff.jar download -v <genome_name> \
	-c <path_to_MicroGMT>/snpEff.config -dataDir <path_to_MicroGMT>/database
```
Please make sure to use "-c" and "-dataDir" to direct the download to MicroGMT directory!

Then you may use "-r <database_name>" in annotate_vcf.py to use the downloaded database.

**Caution: make sure the chromosome name in the downloaded database is the same with that in your fasta reference genome file. If they don't match, no annotation will be produced for vcf outputs and summary tables. Check if accession number is in the fasta header if they don't match.**

<br>

If the genome is not supported, you need to build your own database. The following steps are modified from [snpEff manual](http://snpeff.sourceforge.net/SnpEff_manual.html#databases) to create the database. Here we use SARS-CoV-2 as an example to show the process:
#### 1. Configure the new genome in the configration file provided by MicroGMT: <path_to_MicroGMT>/snpEff.config:
Open the file:
```bash
vi <path_to_MicroGMT>/snpEff.config
```
&#160;Add your genome information into the file.
```bash
# SARS-CoV-2, version NC_045512.2
NC_045512.genome : SARS-CoV-2
```

#### 2. If the genome uses a non-standard codon table: Add codon table parameter. No need for SARS-CoV-2.

#### 3. Get genome annotations. 
Four different formats are accepted: GTF, GFF, RefSeq table from UCSC, and GenBank file. The SARS-CoV-2's annotation file we used is GenBank file downloaded from https://www.ncbi.nlm.nih.gov/nuccore/nc_045512. Rename it by "genes.gbk". Create a folder named "NC_045512" under <path_to_MicroGMT>/database/. Finally, put "genes.gbk" under <path_to_MicroGMT>/database/NC_045512. For other annotation file formats, you will also need the fasta reference genome file. Please see snpEff's manual about how to use GTF, GFF or RefSeq table from UCSC to create database. **Caution: The GenBank file of SARS-CoV-2 was modified to handle the multiple CDS and the -1 ribosomal frameshift of ORF1ab. If your genome has similar issues (i.e. one gene with multiple overlapping CDS), you need to revise your annotation file accordingly.**

* Here we will illustrate how we revised the annotation file "<path_to_MicroGMT>/database/NC_045512/genes.gbk" downloaded from NCBI for SARS-CoV-2:

There are two CDS for ORF1ab gene: 1) CDS joining (266..13468,13468..21555), on which the -1 ribosomal frameshift occurs during translation, produces pp1ab; 2) CDS of (266..13483) produces pp1a. They are paritally overlapped. If the annotation is not revised, the mutation position on CDS and peptides will be shifted for positions after 13468 because snpEff (one software used in MicroGMT) fused the two CDS together. The DNA and amino acid changes, and the mutation loci on the genome will still be correct.

This issue is cause by that the two CDS share the same gene ID and name, so snpEff fused them together. To prevent this, for the first CDS, we added "_pp1ab" to its gene name and locus tag (line 86-87 in "genes.gbk"); for the second CDS, we added "_pp1a" to its gene name and locus tag (line 322-323 in "genes.gbk"). After building the database with revised annotation, the resulting summary tables for running SARS-CoV-2 sequences will have the following attributes: For mutations occur on mature peptides produced by both pp1a and pp1ab, or by pp1a only, the gene ID and name in output vcfs and summary tables are: GU280_gp01_pp1a and ORF1ab_pp1a. For mutations occur on mature peptides produced by pp1ab only, the gene ID and name in output vcfs and summary tables are: GU280_gp01_pp1ab and ORF1ab_pp1ab.

If you have any questions on how to build your own databases on genomes having similar issues, you are more than welcome to email me at yue.july.xing@gmail.com and I'm more than glad to help.

#### 4. Create the database:
```bash
java -jar <path_to_snpEff>/snpEff.jar \
	build -genbank -c <path_to_MicroGMT>/snpEff.config \
	-dataDir <path_to_MicroGMT>/database -v NC_045512
``` 
Please make sure to use "-c" and "-dataDir" to direct the download to MicroGMT directory!

You may also add more annotation information to create the database. Please see [snpEff's manual](http://snpeff.sourceforge.net/SnpEff_manual.html#databases) for more information on building the annotation database.

You will also need the fasta format reference genome sequence file for running MicroGMT. For example, the SARS-CoV-2's reference genome sequence is downloaded from https://www.ncbi.nlm.nih.gov/nuccore/nc_045512. Remove the version number in fasta header: change ">NC_045512.2" to ">NC_045512".

**Caution: make sure the chromosome name in the downloaded database is the same with that in your fasta reference genome file. If they don't match, no annotation will be produced for vcf outputs and summary tables. Check if accession number is in the fasta header if they don't match.**

Another example is in Tutorial section.

## Quick start
### Check environment
```bash
<path_to_MicroGMT>/Check_environment.sh \
	<path_to_snpEff> <path_to_GATK3.8> <path_to_PICARD>
```
Output will be print to the screen to let you know if any tool is not properly installed, or if the tool version is not correct.

### Running MicroGMT for fasta formatted database sequences
#### For SART-CoV-2:
Input: fasta_assembly_file, region_file

Output: mutation annotated vcf files, mutation summary tables. See "output" section above for details.

```bash
# Step 1
python <path_to_MicroGMT>/sequence_to_vcf.py \
  -r <path_to_MicroGMT>/NC_045512_source_files/NC_045512.fa \
  -i assembly -fs <fasta_assembly_file> \
  -o <out_dir_1>

# Step 2
python <path_to_MicroGMT>/annotate_vcf.py \
  -i <out_dir_1> -c -o <out_dir_2> \
  -rg <region_file> -f both -cf \
  -eff <path_to_snpEff>
```
You may add "-sa" in step 2 to skip annotate vcf files, but only output summary tables and/or mutation frequency summaries. Then the input direcotry must contain mutation annotated vcf files by MicroGMT/snpEff. Same for user-supplied genomes and raw reads inputs.

You may remove "-cf" to not further output mutation frequency summaries from summary tables. That will increase speed. Same for user-supplied genomes and raw reads inputs.

#### For user-supplied genomes:
Input: fasta_assembly_file, fasta_reference_sequence_file, reference_genome_database, region_file

Output: mutation annotated vcf files, mutation summary tables. See "output" section above for details.

```bash
# Step 1
python <path_to_MicroGMT>/sequence_to_vcf.py \
  -r <fasta_reference_sequence_file> \
  -i assembly -fs <fasta_assembly_file> \
  -o <out_dir_1>

# Step 2
python <path_to_MicroGMT>/annotate_vcf.py \
  -i <out_dir_1> -c -o <out_dir_2> \
  -rg <region_file> -f both \
  -r <name_of_reference_genome_database> \
  -eff <path_to_snpEff>
```

### Running MicroGMT for fastq formatted raw read sequences
#### For SART-CoV-2:
Input: fastq_raw_reads file(s) (single or paired-end), region_file

Output: mutation annotated vcf files, mutation summary tables. See "output" section above for details.

* For step 1, to run one sample, do the following.
```bash
python <path_to_MicroGMT>/sequence_to_vcf.py \
  -r <path_to_MicroGMT>/NC_045512_source_files/NC_045512.fa \
  -i fastq -fq1 <fastq_raw_reads_R1_file> -fq2 <fastq_raw_reads_R2_file> \
  -o <out_dir_1> \
  -gatk <path_to_gatk> \
  -picard <path_to_picard> \
  -l <log_name> -n <output_prefix> -ki
```
* For step 1, to run multiple samples at one time, do the following, in which the <fastq_prefix_list> file is a text file containg the prefix of each fastq raw reads file to be processed. One sample per line. In this example, if the fastq files are named test_1.fq and test_2.fq, the prefix is "test". The "-fq1" and "-fq2" paramters take in the full names of the fastq files including absolute or relative path.
```bash
# Step 1
cat <fastq_prefix_list> | while read line
do
  python <path_to_MicroGMT>/sequence_to_vcf.py \
  -r <path_to_MicroGMT>/NC_045512_source_files/NC_045512.fa \
  -i fastq -fq1 <dir_to_fqfiles>/${line}_1.fq -fq2 <dir_to_fqfiles>/${line}_2.fq \
  -o <out_dir_1> \
  -gatk <path_to_gatk> \
  -picard <path_to_picard> \
  -l ${line}.log -n ${line} -ki
done
```
&#160;You can also change "-fq1" and "-fq2" to "-fq" to run single end fastq_raw_reads samples.

&#160;A mapping quality summary will be produced in the log file.

&#160;You can filter low quality bases by adding "-m" option. You can fine tune the BWA alignment by revising the "Prepare_fastq_config.txt" file in MicroGMT. It contains parameter settings for BWA. Currently these settings are default. Just change it if you need, and save it.

* Step 2 takes in all the vcf files in a folder produced by step 1 at one time. So you can first process all the samples by step 1, and then process them all together by step 2:
```bash
# Step 2
python <path_to_MicroGMT>/annotate_vcf.py \
  -i <out_dir_1> -c -o <out_dir_2> \
  -rg <region_file_for_raw_read_files> -f both \
  -eff <path_to_snpEff>
```

#### For user-supplied genomes:
Input: fastq_raw_reads file(s) (single or paired-end), fasta_reference_sequence_file, reference_genome_database, region_file

Output: mutation annotated vcf files, mutation summary tables. See "output" section above for details.

* Step 1
```bash
cat <fastq_prefix_list> | while read line
do
  python <path_to_MicroGMT>/sequence_to_vcf.py \
  -r <fasta_reference_sequence_file> \
  -i fastq -fq1 <dir_to_fqfiles>/${line}_1.fq -fq2 <dir_to_fqfiles>/${line}_2.fq \
  -o <out_dir_1> \
  -gatk <path_to_gatk> \
  -picard <path_to_picard> \
  -l ${line}.log -n ${line} -ki
done
```
&#160;Or:
```bash
python <path_to_MicroGMT>/sequence_to_vcf.py \
  -r <fasta_reference_sequence_file> \
  -fq1 <fastq_raw_reads_R1_file> -fq2 <fastq_raw_reads_R2_file> \
  -o <out_dir_1> \
  -gatk <path_to_gatk> \
  -picard <path_to_picard> \
  -l <log_name> -n <output_prefix> -ki
```
&#160;Note: you can also change "-fq1" and "-fq2" to "-fq" to run single end fastq_raw_reads samples.

* Step 2
```bash
python <path_to_MicroGMT>/annotate_vcf.py \
  -i <out_dir_1> -c -o <out_dir_2> \
  -rg <region_file_for_raw_read_files> -f both \
  -r <name_of_reference_genome_database> \
  -eff <path_to_snpEff>
```

### Running MicroGMT for fasta formatted contig sequences
**Warning: This option is not tested. Use at your own risk!**
#### For SART-CoV-2:
Input: contig_sequences_file, region_file

Output: mutation annotated vcf files, mutation summary tables. See "output" section above for details.

```bash
# Step 1
python <path_to_MicroGMT>/sequence_to_vcf.py \
  -r <path_to_MicroGMT>/NC_045512_source_files/NC_045512.fa \
  -i contig -fs <contig_sequences_file> \
  -o <out_dir_1>

# Step 2
python <path_to_MicroGMT>/annotate_vcf.py \
  -i <out_dir_1> -c -o <out_dir_2> \
  -rg <region_file> -f both \
  -eff <path_to_snpEff>
```

#### For user-supplied genomes:
Input: contig_sequences_file, fasta_reference_sequence_file, reference_genome_database, region_file

Output: mutation annotated vcf files, mutation summary tables. See "output" section above for details.

```bash
# Step 1
python <path_to_MicroGMT>/sequence_to_vcf.py \
  -r <fasta_reference_sequence_file> \
  -i contig -fs <fasta_contig_file> \
  -o <out_dir_1>

# Step 2
python <path_to_MicroGMT>/annotate_vcf.py \
  -i <out_dir_1> -c -o <out_dir_2> \
  -rg <region_file> -f both \
  -r <name_of_reference_genome_database> \
  -eff <path_to_snpEff>
```

### Extract strain/IDs from a new fasta database sequence file to make a list and a fasta file for sequences that are not in an existing summary tables
Input: new_database_sequences_file, id_list_for_existing_summary_tables

Output: new_id_list_for_sequences_in_the_new_database_sequence_file, file_containg_ids_for_sequences_not_in_existing_summary_tables, fasta_file_for_sequences_not_in_existing_summary_tables

For the last two outputs, please make sure there's no files with the same names exist in the output directory before running the command!

```bash
<path_to_MicroGMT>/Find_new_seqs.sh \
	<new_database_sequences_file> <id_list_for_existing_summary_tables> \
	<name_of_new_id_list_for_sequences_in_the_new_database_sequence_file> \
	<name_of_file_containg_ids_for_sequences_not_in_existing_summary_tables> \
	<name_of_fasta_file_for_sequences_not_in_existing_summary_tables>
```

### Extract regions for above strain/IDs from a big region file
Input: file_containg_ids_for_sequences_not_in_existing_summary_tables, input_region_information_file

Output: region_information_file_containg_ids_for_sequences_not_in_existing_summary_tables

For the output, please make sure there's no file with the same name exist in the output directory before running the command!

```bash
<path_to_MicroGMT>/Find_regiosn_for_new_seqs.sh \
	<file_containg_ids_for_sequences_not_in_existing_summary_tables> <input_region_information_file> \
	<name_of_region_information_file_containg_ids_for_sequences_not_in_existing_summary_tables>
```

### Remove strains/sequences from summary tables
Input: summary table (**Only \<prefix>.all.form1.txt or \<prefix>.all.form2.txt is required**), list of strain/sequence ids to be removed

Output: **All** kinds of summary tables of form1 or form 2 (depending on your input table form), with strains/sequences removed

Remove strains from format 1 summary tables:
```bash
python <path_to_MicroGMT>/remove_from_summary_tables.py \
  -i <prefix>.all.form1.txt -r <remove_id_list> \
  -f a -d <remove_out_dir>
```
Remove strains from format 2 summary tables:
```bash
python <path_to_MicroGMT>/remove_from_summary_tables.py \
  -i <prefix>.all.form2.txt -r <remove_id_list> \
  -f b -d <remove_out_dir>
```

### Combine summary tables
Input: summary tables (**Only \<prefix>.all.form1.txt or \<prefix>.all.form2.txt is required. The input summary tables need to be in a same form**)

Output: **All** kinds of combined summary tables of form1 or form 2 (depending on your input table form)

Combine format 1 summary tables:
```bash
python <path_to_MicroGMT>/combine_summary_tables.py \
  -d <combine_out_dir> -f a \
  -i1 <prefix_for_input_table_1>.all.form1.txt \
  -i2 <prefix_for_input_table_2>.all.form1.txt
```
Combine format 2 summary tables:
```bash
python <path_to_MicroGMT>/combine_summary_tables.py \
  -d <combine_out_dir> -f b \
  -i1 <prefix_for_input_table_1>.all.form2.txt \
  -i2 <prefix_for_input_table_2>.all.form2.txt
```
### Extract strain/sequence IDs from summary tables
Input: a form2 summary table

Output: a strain/sequence id list for that summary table (it's also the strain/sequence id list for the form1 summary table of the same dataset)

```bash
<path_to_MicroGMT>/get_ids.sh \
	<form2_summary_table> <name_of_id_list>
```
For the output, please make sure there's no file with the same name exist in the output directory before running the command!

### Mask regions on the genome for summary tables
Input: form2 summary table (**Only \<prefix>.all.form2.txt is required**), mask file (tab-delimited, "chr  mask_start  mask_end" one per line)

Output: **All** kinds of masked form2 summary tables.

```bash
python <path_to_MicroGMT>/mask_sequences.py \
  -i <prefix>.all.form2.txt -d <out_dir> -m <input_mask_file>
```

### Add custom annotations
#### For SARS-CoV-2:
Input: form2 summary table (**Only \<prefix>.all.form2.txt is required**)

Output: **All** kinds of form2 summary tables with custom annotations added as a column.

```bash
python <path_to_MicroGMT>/add_custom_annotation.py \
  -i <prefix>.all.form2.txt -d <out_dir> \
  -a <path_to_MicroGMT>/NC_045512_source_files/NC_045512_cus_anno.txt
```
You may also use your own custom annotation file (tab-delimited, "chr  feature_start  feature_end  feature_name" one per line). The annotation file supports multiple annotation lines for one region. Different annotations will be combined by "|" in the output. 

#### For user-supplied genomes:
Input: form2 summary table (**Only \<prefix>.all.form2.txt is required**), custom annotation file (tab-delimited, "chr  feature_start  feature_end  feature_name" one per line) The annotation file supports multiple annotation lines for one region. Different annotations will be combined by "|" in the output. 

Output: **All** kinds of form2 summary tables with custom annotations added as a column.

```bash
python <path_to_MicroGMT>/add_custom_annotation.py \
  -i <prefix>.all.form2.txt -d <out_dir> -a <custom_annotation_file>
```

### Downstream utilities: analysis_utilities.py
#### Reformat summary tables:
* Input: form2 summary table or custom annotated form2 summary table (input one summary table at a time)
* Output: tab delimited file with one mutation per line ("chr pos gene_id gene_name mutation strain_ID region" for form2 summary table, and "chr pos gene_id gene_name custom_annotation mutation strain_ID region" for custom annotated form2 summary table).

```bash
# form2 summary table as input
python <path_to_MicroGMT>/analysis_utilities.py \
  -i <form2_summary_table> -o <name_of_output_table> \
  -t a -a n

# custom annotated form2 summary table as input
python <path_to_MicroGMT>/analysis_utilities.py \
  -i <custom_annotated_form2_summary_table> -o <name_of_output_table> \
  -t a -a y
```

#### Find unqiue mutations (unqiue mutations are defined by only one strain/ID has that mutation at a specific locus):
* Input: form2 summary table or custom annotated form2 summary table (input one summary table at a time)
* Output: tab delimited file with one mutation per line ("strain_ID region chr pos gene_id gene_name mutation" for form2 summary table, and "strain_ID region chr pos gene_id gene_name custom_annotation mutation" for custom annotated form2 summary table).

```bash
# form2 summary table as input
python <path_to_MicroGMT>/analysis_utilities.py \
  -i <form2_summary_table> -o <name_of_output_table> \
  -t b -a n

# custom annotated form2 summary table as input
python <path_to_MicroGMT>/analysis_utilities.py \
  -i <custom_annotated_form2_summary_table> -o <name_of_output_table> \
  -t b -a y
```

### Extract mutation informtation for a specific strain/sequence ID
* Input: form2 summary table or custom annotated form2 summary table (input one summary table at a time)
* Output: tab delimited file with mutation informtation for a user-supplied strain/sequence ID. Long form: "ID  region  chr  pos  gene_id  gene_name  custom_annotation (if input custom annotated form2 summary table)  mutation  mutation_type" + other strains/sequences' mutation information on the mutation loci of this strain/sequence. Short form: "ID  region  chr  pos  gene_id  gene_name  custom_annotation (if input custom annotated form2 summary table)  mutation  mutation_type". There are three kinds of mutation_types: Unique_mutation (On that locus, the mutation for the extracted strain/sequence is not in any other strains in the summary table); New_mutation (On that locus, there are mutations exist for other strains. But the mutation for the extracted strain/sequence is different than all other existing mutations on that locus in the summary table.); Known_mutation (On that locus, there are mutations exist for other strains. And the mutation for the extracted strain/sequence is the same as one of the existing mutations on that locus in the summary table.).

```bash
# form2 summary table as input
python <path_to_MicroGMT>/sequence_ID_extractor.py \
  -i <form2_summary_table> -o <name_of_output_file> \
  -id <user-supplied_strain_or_sequence_ID> -f <form_of_output> -a n

# custom annotated form2 summary table as input
python <path_to_MicroGMT>/sequence_ID_extractor.py \
  -i <custom_annotated_form2_summary_table> -o <name_of_output_file> \
  -id <user-supplied_strain_or_sequence_ID> -f <form_of_output> -a y
  ```

## Tutorial
### 1. Example workflow for SARS-CoV-2 sequences
#### Fasta assembly file as input
**Preprocessing**

Download the desired fasta assembly sequences from [NCBI](https://www.ncbi.nlm.nih.gov/) and put them in one fasta file (named "sequences.fasta" in this example). The fasta headers should be strain/sequence IDs.

Create the regional information file. It is named "metadata.tsv" in this example. The formate is "strain/sequence_ID  region", separated by tab.

If you are using files downloaded from GISAID, follow the steps below for preprocessing: Download the fasta assembly sequences and metadata containing region information. The fasta assembly file used strain ID as the fasta header. Since the strain IDs contain "/", we need to substitute them with something else ("_" in our example). We also need to extract region information from metadata to make the region file, and substitute blanks (" ") in the region file:
```bash
sed -i 's/\//_/g' sequences.fasta

awk -F "\t" '{print $1"\t"$6}' metadata.tsv > metadata.short.tsv
sed -i 's/\//_/g' metadata.short.tsv
sed -i 's/ /_/g' metadata.short.tsv
```

**Exclude strains/IDs that are already exist in the pre-built summary tables for the new inputs (optional)**

If you have a new fasta assembly file and would like to compare with the existing summary tables to remove existing strains/IDs from it first, we provided utility scripts to do this job conveniently. **This is especially useful for excluding strains/IDs already exist in the pre-built summary tables built from the big input fasta assembly file downloaded from GISAID.**

All you need are the new fasta assembly file (sequences.fasta from last step) and the id.list file containing all the strains/IDs from the existing summary tables, which is provided along with the pre-built summary tables. For user supplied genomes, this id.list file is also produced by sequence_to_vcf.py.

MicroGMT will output a list containing all strains/IDs in the new fasta assembly file (sequences.list in this example), a list containing strains/IDs in the new fasta assembly file that are not exist in the pre-built summary tables (ids_to_add.list in this example) and a new fasta assembly file without strains/IDs in the pre-built summary tables (ids_to_add.fasta in this example).
```bash
<path_to_MicroGMT>/Find_new_seqs.sh \
	sequences.fasta id.list \
	sequences.list ids_to_add.list \
	ids_to_add.fasta
```
You can also extract region information for these IDs from the region file (metadata.short.tsv from last step). It is not necessary, but will increase speed a little bit. You can still use the big region file (metadata.short.tsv from last step). The output region file name is ids_to_add.tsv for this example.
```bash
<path_to_MicroGMT>/Find_regiosn_for_new_seqs.sh \
	ids_to_add.list metadata.short.tsv \
	ids_to_add.tsv
```

**Make summary tables**

Use files from last step to make summary tables:
```bash
python <path_to_MicroGMT>/sequence_to_vcf.py \
	-r <path_to_MicroGMT>/NC_045512_source_files/NC_045512.fa \
  -i assembly -fs ids_to_add.fasta \
	-o <make_out_dir_1>
	
python <path_to_MicroGMT>/annotate_vcf.py \
	-i <make_out_dir_1> -c -o <make_out_dir_2> \
	-rg ids_to_add.tsv -f both \
	-eff <path_to_snpEff>
```
The outputs are all the summary tables of format 1 and format 2 for ids_to_add.fasta, log files, as well as the id.list file which contains all the strain IDs in the ids_to_add.fasta file.

**Remove strains/IDs from summary tables  (optional)**

We noticed that strains may be removed from the GISAID SARS-CoV-2 database. So we designed this utility script to remove unwanted strains from summary tables. You will need a list of strains/IDs that need to be removed. Here we will demostrate how to use it to remove unwanted strains from the pre-built summary tables:

If you have a list of IDs in file A (sequences.list from last step), the existing summary tables (the pre-built summary tables in this example), and you want to identify unwanted strains (that is, strains in the pre-built summary tables but not in file A), you may use the following commands. **This is especially useful for excluding strains/IDs already exist in the pre-built summary tables built from the big input fasta assembly file downloaded from GISAID.**
```bash
cat <path_to_summary_tables>/id.list | while read line
do
  if grep -q "^${line}$" sequences.list
  then
    echo "Sequence ${line} won't be removed, skip."
  else
    echo ${line} >> remove.list
  fi
done
```
Remove strains from format 1 summary tables:
```bash
python <path_to_MicroGMT>/remove_from_summary_tables.py \
	-i <path_to_summary_tables>/all0520.all.form1.txt \
	-r remove.list -p removed -l removed_a \
	-f a -d <remove_out_dir>
```
Remove strains from format 2 summary tables:
```bash
python <path_to_MicroGMT>/remove_from_summary_tables.py \
	-i <path_to_summary_tables>/all0520.all.form2.txt \
	-r remove.list -p removed -l removed_b \
	-f b -d <remove_out_dir>
```

**Combine summary tables  (optional)**

We will demonstrate how to combine the summary tables from "Make summary tables" and "Remove strains/IDs from summary tables" sessions above. 

Combine format 1 summary tables:
```bash
python <path_to_MicroGMT>/combine_summary_tables.py \
	-d <combine_out_dir> -f a -p combined -l combined_a \
	-i1 <remove_out_dir>/removed.all.form1.txt \
	-i2 <make_out_dir_2>/out_summary.all.form1.txt
```
Combine format 2 summary tables:
```bash
python <path_to_MicroGMT>/combine_summary_tables.py \
	-d <combine_out_dir> -f b -p combined -l combined_b \
	-i1 <remove_out_dir>/removed.all.form2.txt \
	-i2 <make_out_dir_2>/out_summary.all.form2.txt
```

**Make a new strain/ID list for use next time (optional)**

We will demostrate an optional step of making a new strain/ID list for use next time (final.list in this example). This list contains all strain/IDs in the final output summary tables. Users can use it as the input list file for removing or adding strains/IDs to the new summary tables in the future.

Please make sure there is no file named "final.list" in your directory before we start.
```bash
cat id.list ids_to_add.list > tmp.list
cat tmp.list | while read line
do
  if grep -q "^${line}$" remove.list
  then
    echo "Sequence ${line} is removed."
  else
    echo ${line} >> final.list
  fi
done
rm -f tmp.list
```

Alternatively, you may just use get_ids.sh:
```bash
<path_to_MicroGMT>/get_ids.sh \
	out_summary.all.form2.txt final.list
```

#### Fastq raw reads files as input
Here we produce summary tables for the simulated fastq raw reads files from 5 strains. The prefix for the fastq files are in the file "ids_for_5_strains.list". The IDs in the summary tables are the prefix for the fastq files in this example.
```bash
cat ids_for_5_strains.list | while read line
do
  python <path_to_MicroGMT>/sequence_to_vcf.py \
  -r <path_to_MicroGMT>/NC_045512_source_files/NC_045512.fa \
  -i fastq -fq1 <dir_to_fqfiles>/${line}_1.fq -fq2 <dir_to_fqfiles>/${line}_2.fq \
  -o <out_dir_1> \
  -gatk <path_to_gatk> \
  -picard <path_to_picard> \
  -l ${line}.log -n ${line} -ki
done
	
python <path_to_MicroGMT>/annotate_vcf.py \
  -i <out_dir_1> -c -o <out_dir_2> \
  -rg 5_strains_region_file.tsv -f both \
  -eff <path_to_snpEff>
```

You may combine the summary tables produced from both the fasta assembly file and the fastq raw reads files.

#### Add custom annotations (optional)
For SARS-CoV-2, a custom annotation file is provided at <path_to_MicroGMT>/NC_045512_source_files/NC_045512_cus_anno.txt. It contains mature peptide and stem loop information of SARS-CoV-2.

```bash
python <path_to_MicroGMT>/add_custom_annotation.py \
  -i <prefix>.all.form2.txt -d <out_dir> \
  -a <path_to_MicroGMT>/NC_045512_source_files/NC_045512_cus_anno.txt
```

#### Downstream utilities (optional)
Reformat custom annotated summary tables:
```bash
python <path_to_MicroGMT>/analysis_utilities.py \
  -i <custom_annotated_form2_summary_table> -o <name_of_output_table> \
  -t a -a y
```
Find unqiue mutations in custom annotated summary tables (unqiue mutations are defined by only one strain/ID has that mutation at a specific locus):
```bash
python <path_to_MicroGMT>/analysis_utilities.py \
  -i <custom_annotated_form2_summary_table> -o <name_of_output_table> \
  -t b -a y
```

#### Extract mutation informtation for a specific strain/sequence ID from custom annotated summary tables (optional)
```bash
python <path_to_MicroGMT>/sequence_ID_extractor.py \
  -i <custom_annotated_form2_summary_table> -o <name_of_output_file> \
  -id <user-supplied_strain_or_sequence_ID> -f l -a y
```

#### Test dataset
If you want to try on a small test dataset first, please use the provided test datasets. Simply unzip it to use. It contains database sequences and simulated raw read sequnces of 5 randomly selected strains from NCBI database. For fasta assembly sequences, start from "Make summary tables". The formatting steps in "Preprocessing" section is already done. 

### 2. Workflow for sequences of E.coli K12 strains
#### Build the annotation database
Here we illustrate how to build a new database for E.coli K12 strains by the Genbank file of the reference sequence.
* Download the full Genbank file of E.coli K12 reference NC_000913 at https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3/.
* Rename it by "genes.gbk".
* Create a new directory with the name of reference under <path_to_MicroGMT>/database/: here we create a directory named "NC_000913" under <path_to_MicroGMT>/database/.
* Put "genes.gbk" in <path_to_MicroGMT>/database/NC_000913.
* Open "snpEff.config" under <path_to_MicroGMT> and add the following. Note: do NOT open "snpEff.config" under the snpEff directory!
```bash
# E.coli_k12 genome, version NC_000913.3
NC_000913.genome : E.coli_k12
````
* Run the following command to build the database:
```bash
java -jar <path_to_snpEff>/snpEff.jar build -genbank -c <path_to_MicroGMT>/snpEff.config \
-dataDir <path_to_MicroGMT>/database -v NC_000913
```
You will also need the fasta formatted reference sequence file. It is downloaded from https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3/ as well. **Important: Your fasta reference sequence file should not contain version number in the header line.** Delete the version number in header line: change ">NC_000913.3" to ">NC_000913"!

#### Fasta assembly file as input
The input fasta assembly test file contains 3 E.coli fasta sequences downloaded from NCBI. We just illustrate the core steps here.
```bash
# Step 1
python <path_to_MicroGMT>/sequence_to_vcf.py \
  -r <fasta_reference_sequence_file> \
  -i assembly -fs <fasta_assembly_file> \
  -o <out_dir1>

# Step 2
python <path_to_MicroGMT>/annotate_vcf.py \
  -i <out_dir1> -c -p <output_prefix> -r NC_000913 \
  -o <out_dir2> -f both -eff <path_to_snpEff>

# Find unique mutations
python <path_to_MicroGMT>/analysis_utilities.py \
  -i <out_dir2>/<output_prefix>.all.form2.txt \
  -o <output_table_name> -t b
```

#### Fastq raw reads files as input
You may use the fasta sequence of an E.coli K12 strain to simulate fastq raw read sequences and use them as test data. Example using [ART-illumia](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) for simulation is shown here. It is using HiSeq2500, paired-end, read length of 150 bps, fragment size of 300 bps with standard deviation of 20 bps, fold coverage 0f 50x, and the base quality between 18 to 38.
```bash
./art_illumina -i <fasta_sequence_of_E.coli_K12_strain> \
	-ss HS25 -p -l 150 -f 50 --mflen 300 -s 20 -na \
	-qL 18 -qU 38 \
	-o <prefix_of_fastq_raw_reads_file>
```

The core MicroGMT steps are listed here:
```bash
# Step 1
python <path_to_MicroGMT>/sequence_to_vcf.py \
  -r <fasta_reference_sequence_file> -i fastq -fq1 <fastq_raw_reads_R1_file> \
  -fq2 <fastq_raw_reads_R2_file> -o <out_dir_1> \
  -gatk <path_to_gatk> -picard <path_to_picard> -ki

# Step 2
python <path_to_MicroGMT>/annotate_vcf.py \
  -i <out_dir_1> -c -p <output_prefix> \
  -o <out_dir_2> \
  -f both -r NC_000913 \
  -eff <path_to_snpEff>

# Find unique mutations
python <path_to_MicroGMT>/analysis_utilities.py \
  -i <out_dir2>/<output_prefix>.all.form2.txt \
  -o <output_table_name> -t b
```

#### Mask regions on the genome for summary tables
For example, mask repeat regions.
```bash
python <path_to_MicroGMT>/mask_sequences.py \
  -i <prefix>.all.form2.txt -d <out_dir> -m <input_mask_file>
```

## Other things you need to know:
* "-"s in input fasta sequences are interpreted as "N"s by MicroGMT. If they represent gaps, they should be removed from fasta sequences.
* The coding of indels in DNA sequence and CDS sequence are slightly different for the same mutations identified by fasta formatted inputs and fastq formatted inputs. But the mutation type, effect, position on DNA and the mutation of amino acids are annotated the same. So, if you want to find unqiue indels across both fasta and fastq formatted inputs, it is suggested to transform them to a same format first.
* For the summary tables of amino acid changes, if there's a mutation but no amino acid change, it is written as blank ("") in the summary table. When reformat by analysis_utilities.py, it is not included in the output table.
* If your fasta assembly file is named by the the sequence ID or one of the sequence IDs in the fasta header, pleas do not direct your output into the same folder as the fasta assembly file. Otherwise your fasta assembly file will be deleted!
* MicroGMT is designed for tracking indels and SNPs among closely related strains instead of detecting large-scale complex genomic rearrangements and duplications. In addition, if you supply fastq raw reads file as input, the accuracy of mutation detection can be slightly affected by unmasked repetitive regions in the reference genome due to the difficult nature of aligning short sequence reads to the reference genome. You may mask repeat regions by mask_sequences.py.

## Run time
The total time used for the core steps for one SARS-CoV-2 sample is:
* 4s for one sequence in the fasta assembly file
* 14s for one fastq raw reads input with 50x depth of coverage

## Arguments
### Check_environment.sh
```bash
<path_to_MicroGMT>/Check_environment.sh \
	<path_to_snpEff> <path_to_GATK3.8> <path_to_PICARD>
```
### sequence_to_vcf.py
``` bash
usage: sequence_to_vcf.py [-h] -r REF_GENOME -i {assembly,contig,fastq} -o
                          OUT_DIR [-fs FASTA_SEQS] [-fq1 FASTQ1] [-fq2 FASTQ2]
                          [-fq FASTQ] [-l LOG] [-n NAME] [-gatk PATH_TO_GATK]
                          [-picard PATH_TO_PICARD] [-kb] [-ki] [-p PRIOR]
                          [-m MBQ] [-a {asm5,asm10,asm20}] [-t THREAD]

Sequence file(s) to vcf file(s)

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
  -n NAME               Name of the input sample. Does not work with 'assembly'                                                                              option. [test]
  -gatk PATH_TO_GATK    Absolute path to GenomeAnalysisTK.jar. Only required for                                                                              'fastq' option.
  -picard PATH_TO_PICARD
                        Absolute path to picard.jar. Only required for 'fastq' o                                                                             ption.
  -kb                   Keep BAM files.
  -ki                   Keep index files. Only works with 'fastq' option.
  -p PRIOR              Prior for bcftools variant caller (expected substitution                                                                              rate). 0 means the prior is disabled. Only works for 'assembly' or 'contig' opt                                                                             ion. [0]'.
  -m MBQ                Minimum base quality for variant caller. Only works with                                                                              'fastq' option. [10]
  -a {asm5,asm10,asm20}
                        Sequence divergence: asm5/asm10/asm20 for ~0.1/1/5 perce                                                                             ntages. Only works with 'assembly' option. [asm5]
  -t THREAD             Number of threads. [10]
```

### annotate_vcf.py
``` bash
usage: annotate_vcf.py [-h] -i IN_DIR -o OUT_DIR [-r REFERENCE]
                       [-p TABLE_PREFIX] [-c] [-l LOG] [-f {a,b,both}]
                       [-rg REGION_FILE] [-sa] [-eff PATH_TO_SNPEFF] [-cf]

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
  -sa                  Skip vcf annotation step, just make summary tables from annotated vcf files in the input directory. The input vcf files must be MicroGMT/snpEff annotated.
  -eff PATH_TO_SNPEFF  Absolute path to snpEff.jar. Required if annotae vcf files.
  -cf                  Calculate frequency summaries for the summary tables.
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

### Find_new_seqs.sh
```bash
<path_to_MicroGMT>/Find_new_seqs.sh \
  <new_database_sequences_file> <id_list_for_existing_summary_tables> \
  <name_of_new_id_list_for_sequences_in_the_new_database_sequence_file> \
  <name_of_file_containg_ids_for_sequences_not_in_existing_summary_tables> \
  <name_of_fasta_file_for_sequences_not_in_existing_summary_tables>
```
 
### Find_regiosn_for_new_seqs.sh 
```bash
<path_to_MicroGMT>/Find_regiosn_for_new_seqs.sh \
  <file_containg_ids_for_sequences_not_in_existing_summary_tables> <input_region_information_file> \
  <name_of_region_information_file_containg_ids_for_sequences_not_in_existing_summary_tables>
```
### get_ids.sh
```bash
<path_to_MicroGMT>/get_ids.sh \
	<form2_summary_table> <name_of_id_list>
```

### mask_sequences.py
```bash
usage: mask_sequences.py [-h] -i IN_TABLE -m MASK -d OUT_DIR [-p OUT_PREF]
                         [-l LOG]
			 
optional arguments:
  -h, --help   show this help message and exit

Mandatory inputs:
  -i IN_TABLE  Input summary table (Only ".all.form2.txt" tables)
  -m MASK      Input mask file
  -d OUT_DIR   Output directory

Optional arguments:
  -p OUT_PREF  Prefix of the output annotated summary tables. Do not include path, except for folder name(s) inside output directory!
  -l LOG       Name of the log file [Mask_sequences.log]
```

### add_custom_annotation.py
```bash
usage: add_custom_annotation.py [-h] -i IN_TABLE -a CUSTOM_ANNOT -d OUT_DIR
                                [-p OUT_PREF] [-l LOG]

optional arguments:
  -h, --help       show this help message and exit

Mandatory inputs:
  -i IN_TABLE      Input summary table (Only ".all.form2.txt" tables)
  -a CUSTOM_ANNOT  Input custom annotation table
  -d OUT_DIR       Output directory

Optional arguments:
  -p OUT_PREF      Prefix of the output annotated summary tables. Do not include path, except for folder name(s) inside output directory!
  -l LOG           Name of the log file [Add_custom_annotation.log]
```

### analysis_utilities.py
```bash
usage: analysis_utilities.py [-h] -i IN_TABLE -o OUT_TABLE -t {a,b} -a {y,n}
                             [-l LOG]

optional arguments:
  -h, --help    show this help message and exit

Mandatory inputs:
  -i IN_TABLE   Input summary table (format 2 table)
  -o OUT_TABLE  Processed output table
  -t {a,b}      Type of analysis (a: format change, b: find unique mutations)
  -a {y,n}      The summary table include custom annotation or not? (y: yes, n: no)

Optional arguments:
  -l LOG        Directory and name of the log file [Analysis_utilities.log]

```
### sequence_ID_extractor.py
```bash
usage: sequence_ID_extractor.py [-h] -i IN_TABLE -o OUT_TABLE -id SEQ_ID -f
                                {l,s} -a {y,n} [-l LOG]

optional arguments:
  -h, --help    show this help message and exit

Mandatory inputs:
  -i IN_TABLE   Input summary table (format 2 table)
  -o OUT_TABLE  Processed output table
  -id SEQ_ID    Strain/sequence ID to extract
  -f {l,s}      The output table format (l: long format, s: short format)
  -a {y,n}      The summary table include custom annotation or not? (y: yes, n: no)

Optional arguments:
  -l LOG        Directory and name of the log file [ID_extraction.log]
```

## Questions and bug report
Please direct all questions and bug reports to Yue Xing at: yue.july.xing@gmail.com

