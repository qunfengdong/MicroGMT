#!/bin/bash
# MicroGMT Version 1.3  (June 2020)

eff=$1

# check python
a=$(python --version)
if [[ $? != 0 ]]
then
	echo "Error: Python is not installed properly."
else
	if [[ $a == *3.* ]]
	then
	    echo "Python version passed."
	else
	    if [[ $a == "" ]]
	    then
	    	echo "Python version is printed on screen. If your python version is 3 you are good to go."
	    else
	    	echo "Error: Python version is not passed. Need version 3."
	    fi
	fi
fi

# check snpeff
a=$(java -jar ${eff}/snpEff.jar -version)
if [[ $? != 0 ]]
then
	echo "Error: SnpEff is not installed properly."
else
	if [[ $a = *4.3t* ]]
	then
	    echo "SnpEff version passed."
	else
	    echo "Error: SnpEff version is not passed. Need version 4.3t. MicroGMT is not tested with SnpEff 5.0."
	fi
fi

# check samtools
if [[ -x "$(command -v samtools)" ]]
then
	a=$(samtools --version)
	if [[ $a = *1.6* ]] || [[ $a = *1.7* ]] || [[ $a = *1.8* ]] || [[ $a = *1.9* ]]  || [[ $a = *1.10* ]]
	then
	    echo "SAMtools version passed."
	else
	    echo "Error: SAMtools version is not passed. Need version 1.6 or above."
	fi
else
	echo "Error: SAMtools is not installed properly."
fi

# check minimap2
if ! [ -x "$(command -v minimap2)" ]
then
	echo "Warning: minimap2 is not installed properly. It is required for fasta genome assembly sequences/contigs inputs."
else
	echo "minimap2 is installed. Passed."
fi

# check bcftools
a=$(bcftools --version)
if [[ $? != 0 ]]
then
	echo "Warning: Bcftools is not installed properly. It is required for fasta genome assembly sequences/contigs inputs."
else
	if [[ $a = *1.6* ]] || [[ $a = *1.7* ]] || [[ $a = *1.8* ]] || [[ $a = *1.9* ]]  || [[ $a = *1.10* ]]
	then
	    echo "Bcftools version passed."
	else
	    echo "Warning: Bcftools version is not passed. Need version 1.6 or above. It is required for fasta genome assembly sequences/contigs inputs."
	fi
fi

# check gatk
a=$(gatk --version)
if [[ $? != 0 ]]
then
	echo "Warning: GATK 4 is not installed properly. It is required for fastq raw reads inputs."
else
	if [[ $a = *v4* ]]
	then
	    echo "GATK version passed."
	else
	    echo "Warning: GATK version is not passed. Need version 4. It is required for fastq raw reads inputs."
	fi
fi

# check bwa
bwa mem > MicroGMT_checkversion.tmp 2>&1 
if ! [ -x "$(command -v bwa)" ]
then
	echo "Warning: BWA is not installed properly. It is required for fastq raw reads inputs."
elif grep -q "unrecognized command 'mem'" MicroGMT_checkversion.tmp
then
	echo "Warning: BWA version is not passed. Need version 0.7 or above. It is required for fastq raw reads inputs."
else
	echo "BWA version passed."
fi
rm -f MicroGMT_checkversion.tmp

# check java
java -version
java -version > MicroGMT_checkversion.tmp 2>&1
if grep -q "1.8" MicroGMT_checkversion.tmp
then
	echo "JAVA version passed."
else
	echo "JAVA version printed on screen. If your JAVA version is 1.8 or above you are good to go."
fi
rm -f MicroGMT_checkversion.tmp
