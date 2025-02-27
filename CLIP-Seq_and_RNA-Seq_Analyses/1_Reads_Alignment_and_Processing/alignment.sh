#!/bin/bash
module load STAR/2.5.3a
DIR_GENOME=/omics/odcf/analysis/OE0284_projects/clipseq/users/gencode_v43
DIR_FASTQ=/home/v096l/varshni/remove_rRNA/fastq_riboRNA_removed
for i in $(ls $DIR_FASTQ/*.fastq.gz)
do
l=$(basename ${i} | cut -c 1-7)
STAR --runThreadN 64 \
--runMode alignReads \
--genomeDir $DIR_GENOME \
--outFilterMismatchNoverReadLmax 0.04 \
--outFilterMismatchNmax 999 \
--outFilterMultimapNmax 30 \
--alignEndsType Extend5pOfRead1 \
--sjdbGTFfile $DIR_GENOME/gencode.v43.annotation.gtf \
--sjdbOverhang 152 \
--outReadsUnmapped Fastx \
--outSJfilterReads Unique \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--readFilesIn $DIR_FASTQ/${l}_wo_rRNA.fastq.gz \
--outFileNamePrefix ${l}_
done
module unload STAR/2.5.3a
