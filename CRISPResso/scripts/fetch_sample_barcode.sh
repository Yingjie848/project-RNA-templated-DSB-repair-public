#!/bin/bash

fastq_path=$1

fastq_filename=$(echo $fastq_path | sed 's/.*\///'); 

totalReads=$(gunzip -c $fastq_path | grep -c ^@)
top1=$(gunzip -c $fastq_path | grep ^@ | cut -d' ' -f2 | cut -d':' -f4 | sort | uniq -c | sort -nr | head -1)
barcode=$(echo $top1 | cut -d' ' -f2)
top1frequency=$(echo $top1 | cut -d' ' -f1)

echo -e "$fastq_filename\t$barcode\t$totalReads\t$top1frequency" >> sample_barcode.txt
