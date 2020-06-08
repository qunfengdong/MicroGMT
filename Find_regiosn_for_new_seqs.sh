#!/bin/bash
# MicroGMT Version 1.3  (June 2020)

add_id_new_list=$1
reg_file=$2
reg_add_file=$3

cat ${add_id_new_list} | while read line
do
	grep $line ${reg_file} >> ${reg_add_file}
done
