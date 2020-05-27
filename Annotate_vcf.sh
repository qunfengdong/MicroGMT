#!/bin/bash

in_dir=$1
out_dir=$2
snpeff_dir=$3
ref=$4
tdir=$5
log=$6
csvS=$7

cd $out_dir

if [[ $csvS == "T" ]]
then
	for fl in $(ls ${in_dir}/*.vcf)
	do
		nm=$(echo ${fl} | rev | cut -d"." -f2- | rev)
		nm=$(echo ${nm} | rev | cut -d"/" -f1 | rev)
		java -Xmx4g -jar ${snpeff_dir}/snpEff.jar \
		-c ${tdir}/snpEff.config \
		-dataDir ${tdir}/database \
		-csvStats ${nm}.snpEff_summary.csv \
		${ref} ${fl} \
		> ${nm}.anno.vcf
	
		if [[ $? != 0 ]]
		then
			echo  "Error: Vcf annotation not successful!"
			echo  "Error: Vcf annotation not successful!" >> $log
			exit 1
		fi
	done
	rm -f snpEff_summary.html
	rm -f snpEff_summary.genes.txt
else
	for fl in $(ls ${in_dir}/*.vcf)
	do
		nm=$(echo ${fl} | rev | cut -d"." -f2- | rev)
		nm=$(echo ${nm} | rev | cut -d"/" -f1 | rev)
		java -Xmx4g -jar ${snpeff_dir}/snpEff.jar \
		-c ${tdir}/snpEff.config \
		-dataDir ${tdir}/database ${ref} \
		${fl} \
		> ${nm}.anno.vcf
	
		if [[ $? != 0 ]]
		then
			echo  "Error: Vcf annotation not successful!"
			echo  "Error: Vcf annotation not successful!" >> $log
			exit 1
		fi
	done
	rm -f snpEff_summary.html
	rm -f snpEff_summary.genes.txt
fi
