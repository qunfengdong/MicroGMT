#!/bin/bash
# MicroGMT Version 1.3  (June 2020)

ref=$1
fq1=$2
fq2=$3
fid=$4
prior=$5
mbq=$6
BAQ=$7
PATH_TO_GATK=$8
out_dir=$9
log=${10}
fastq=${11}
pairing=${12}
th=${13}
mdir=${14}
path_to_picard=${15}
md=${16}


cd $out_dir

a=$(cat ${mdir}/Prepare_fastq_config.txt)

if [[ $pairing == "paired" ]]
then
bwa mem -t ${th} -R "@RG\tID:${fid}\tSM:${fid}" \
${a} \
${ref} ${fq1} ${fq2} | \
samtools view -@ ${th} -Su -q 1 - | \
samtools sort -@ ${th} - | samtools rmdup - ${fid}.bam
else
bwa mem -t ${th} -R "@RG\tID:${fid}\tSM:${fid}" \
${a} \
${ref} ${fastq} | \
samtools view -@ ${th} -Su -q 1 - | \
samtools sort -@ ${th} - | samtools rmdup - ${fid}.bam
fi

#samtools index ${fid}.bam

gatk --java-options "-Xmx4g" SortSam -I ${fid}.bam -O ${fid}.s.bam -SO coordinate

gatk --java-options "-Xmx4g" MarkDuplicatesWithMateCigar -I ${fid}.s.bam -O ${fid}.f.bam \
-M=${fid}.marked_dup_metrics.txt --MINIMUM_DISTANCE ${md} \
--CREATE_INDEX true

#echo "Samtools flagstat for BAM file of ${fid}:" >> $log
samtools flagstat ${fid}.f.bam > ${fid}.flagstat

gatk --java-options "-Xmx4g" HaplotypeCaller  \
-R ${ref} -I ${fid}.f.bam \
-mbq ${mbq} \
-ploidy 1 \
-O ${fid}.vcf

if [ -f ${fid}.vcf ]
then
	rm -f ${ref}
	rm -f ${fid}.bam
	rm -f ${fid}.bam.bai
	rm -f ${fid}.s.bam
	rm -f ${fid}.s.bam.bai
	#rm -f ${fid}.marked_dup_metrics.txt
	#rm -f ${fid}_IndelRealigner.intervals
	mv ${fid}.f.bam ${fid}.bam
	mv ${fid}.f.bai ${fid}.bam.bai
	echo "Processing fastq inputs successful." >> $log
else
	echo "Error: Processing fastq inputs not successful."
	echo "Error: Processing fastq inputs not successful." >> $log
	exit 1
fi
