#!/bin/bash

fl=$1
id_list=$2
add_id_list=$3
add_id_new_list=$4
add_new_out=$5

samtools faidx ${fl}
awk '{print $1}' ${fl}.fai > ${add_id_list}
cat ${add_id_list} | while read line
do
	if grep -q "^${line}$" ${id_list}
	then
		echo "Sequence ${line} was already processed, skip."
	else
    		echo $line >> ${add_id_new_list}
	fi
done

cat ${add_id_new_list} | while read fid
do
	samtools faidx ${fl} ${fid} >> ${add_new_out}
done