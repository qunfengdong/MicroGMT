#!/bin/bash

# fasta assembly input
mkdir faout
mkdir faout/out2

python <path_to_MicroGMT>/sequence_to_vcf.py \
	-r <path_to_MicroGMT>/NC_045512_source_files/NC_045512.fa \
	-i assembly -fs sequences.fasta \
	-o faout

python <path_to_MicroGMT>/annotate_vcf.py \
	-i faout -c -p fasummary \
	-o faout/out2 \
	-rg region.tsv -f both \
	-eff <path_to_SNPEff>

python <path_to_MicroGMT>/add_custom_annotation.py \
	-i faout/out2/fasummary.all.form2.txt \
	-d faout/out2 \
	-a <path_to_MicroGMT>/NC_045512_source_files/NC_045512_cus_anno.txt

python <path_to_MicroGMT>/analysis_utilities.py \
	-i faout/out2/annotated.all.form2.txt \
	-o faout/out2/annotated.all.form2.ref.txt \
	-t a -a y

python <path_to_MicroGMT>/sequence_ID_extractor.py \
	-i faout/out2/annotated.all.form2.txt \
	-o faout/out2/annotated.all.form2.MT734046.1.long.txt  \
	-id MT734046.1 -f l -a y

python <path_to_MicroGMT>/sequence_ID_extractor.py \
	-i faout/out2/annotated.all.form2.txt \
	-o faout/out2/annotated.all.form2.MT734046.1.short.txt  \
	-id MT734046.1 -f s -a y

# fastq raw read inputs
mkdir fqout
mkdir fqout/out2

cat ids | while read line
do
	python <path_to_MicroGMT>/sequence_to_vcf.py \
	-r <path_to_MicroGMT>/NC_045512_source_files/NC_045512.fa \
	-i fastq -fq1 ${line}_1.fq \
	-fq2 ${line}_2.fq \
	-o fqout \
	-gatk <path_to_GATK> \
	-picard <path_to_PICARD> \
	-l ${line}.log \
	-n ${line} -ki -p 0
done

python <path_to_MicroGMT>/annotate_vcf.py \
	-i fqout -c -p fqsummary \
	-o fqout/out2 \
	-rg region.tsv -f both \
	-eff <path_to_SNPEff>
