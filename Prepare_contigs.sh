#!/bin/bash
# MicroGMT Version 1.3.2  (Sep 2020)

ref=$1
seqs_file=$2
fid=$3
out_dir=$4
# mutation rate (use bigger for greater sensitivity), use with -m [1.1e-3]
prior=$5
log=$6
th=$7

cd $out_dir

minimap2 -t ${th} -a ${ref} -R "@RG\tID:${fid}\tSM:${fid}" ${seqs_file} | \
samtools view -@ ${th} -Su - | \
samtools sort -@ ${th} - -o ${fid}.bam

bcftools mpileup --threads ${th} -B -Q 0 -f ${ref} ${fid}.bam | \
bcftools call --threads ${th} -O v -o ${fid}.vcf --ploidy 1 -mv -P ${prior} -

if [ -f ${fid}.vcf ]
then
	echo "Processing fasta assembly inputs successful!" >> $log
else
	echo "Error: Processing fasta assembly inputs not successful."
	echo "Error: Processing fasta assembly inputs not successful." >> $log
	exit 1
fi

