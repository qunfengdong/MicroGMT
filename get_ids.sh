#!/bin/bash

file1=$1
file2=$2

ids=$(awk 'NR==1{print}' ${file1})
ids2=$(echo $ids | sed 's/[ ][ ]*/,/g')

OIFS=$IFS
IFS=','

rm -f ${file2}

for i in $ids2
do
    echo $i >> ${file2}
done

IFS=$OIFS

sed -i 's/|.*//g' ${file2}
sed -i '1,4d' ${file2}
