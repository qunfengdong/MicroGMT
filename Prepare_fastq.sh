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

samtools index ${fid}.bam

java -jar ${PATH_TO_GATK}/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R ${ref} -I ${fid}.bam \
-o ${fid}_IndelRealigner.intervals

java -jar ${PATH_TO_GATK}/GenomeAnalysisTK.jar \
--filter_bases_not_stored -T IndelRealigner \
-R ${ref} -I ${fid}.bam \
-targetIntervals ${fid}_IndelRealigner.intervals \
-o ${fid}.f.bam

echo "Samtools flagstat for BAM file of ${fid}:" >> $log
samtools flagstat ${fid}.f.bam >> $log

java -jar ${PATH_TO_GATK}/GenomeAnalysisTK.jar -T HaplotypeCaller \
-R ${ref} -I ${fid}.f.bam \
-mbq ${mbq} \
-ploidy 1 \
-o ${fid}.vcf

if [ -f ${fid}.vcf ]
then
	rm -f ${ref}
	rm -f ${fid}.bam
	rm -f ${fid}.bam.bai
	rm -f ${fid}_IndelRealigner.intervals
	mv ${fid}.f.bam ${fid}.bam
	mv ${fid}.f.bai ${fid}.bam.bai
	echo "Processing fastq inputs successful." >> $log
else
	echo "Error: Processing fastq inputs not successful."
	echo "Error: Processing fastq inputs not successful." >> $log
	exit 1
fi