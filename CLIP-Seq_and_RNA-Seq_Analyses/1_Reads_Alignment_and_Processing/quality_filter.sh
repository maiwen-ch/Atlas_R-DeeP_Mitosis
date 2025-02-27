#!/bin/bash
DIR=/home/v096l/varshni/unprocessed_fastq/run230322_NB552269_0349_AH3CT5BGXT

mkdir tmp
zcat $DIR/AS-954459-LR-68266_R1.fastq.gz | fastx_trimmer -l 15 | fastq_quality_filter -q 10 -p 100 | awk 'FNR%4==1 { print $1 }' | sed 's/@//' > tmp/data.qualFilteredIDs.list