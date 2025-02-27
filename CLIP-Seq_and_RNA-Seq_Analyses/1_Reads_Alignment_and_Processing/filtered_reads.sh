#!/bin/bash
DIR_FASTQ=/home/v096l/varshni/unprocessed_fastq/run230322_NB552269_0349_AH3CT5BGXT
DIR_LIST=/home/v096l/varshni/quality_filter/tmp

seqtk subseq $DIR_FASTQ/AS-954459-LR-68266_R1.fastq.gz $DIR_LIST/data.qualFilteredIDs.list | sed 's/ /#/g' | gzip > data.filtered.fastq.gz