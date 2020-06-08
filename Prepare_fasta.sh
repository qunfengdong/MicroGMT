#!/bin/bash
# MicroGMT Version 1.3  (June 2020)

ref=$1
seqs_file=$2
seqs=$3
out_dir=$4
#~0.1/1/5% sequence divergence: asm5/asm10/asm20
asm=$5
# mutation rate (use bigger for greater sensitivity), use with -m [1.1e-3]
prior=$6
log=$7
th=$8
keep_bam=$9

cd $out_dir

cp ${seqs_file} ${out_dir}/tmp_${seqs}
samtools faidx tmp_${seqs}
awk '{print $1}' tmp_${seqs}.fai > id.list

cat id.list | while read fid
do
	#echo ${fid}
	samtools faidx tmp_${seqs} ${fid} > ${fid}.fa

	minimap2 -t ${th} -R "@RG\tID:${fid}\tSM:${fid}" -ax ${asm} ${ref} ${fid}.fa | \
	samtools view -@ ${th} -Su - | \
	samtools sort -@ ${th} - -o ${fid}.bam

	bcftools mpileup --threads ${th} -B -Q 0 -f ${ref} ${fid}.bam | \
	bcftools call --threads ${th} -O v -o ${fid}.vcf --ploidy 1 -mv -P ${prior} -

	if [ -f ${fid}.vcf ]
	then
		rm -f ${fid}.fa
		if [[ $keep_bam == "F" ]]
		then
			rm -f ${fid}.bam
		fi

	else
		echo  "Error: Processing fasta assembly inputs not successful."
		echo  "Error: Processing fasta assembly inputs not successful." >> $log
		exit 1
	fi
done

rm -f tmp_${seqs}.fai
rm -f tmp_${seqs}