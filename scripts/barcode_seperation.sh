#!/bin/bash
cur_dir=`dirname $0`
idir=${cur_dir}/../data/raw/
odir=${cur_dir}/../data/fasta/
mkdir -p ${odir}
for f in ${idir}*.fastq.gz; do
    fcore=`basename $f`
    fcore=${fcore%%.*}
    cutadapt -a CTGTAGGCACCATCAAT -m 0 --discard-untrimmed -o ${odir}${fcore}.fastq.gz $f > ${odir}${fcore}.report
    zcat $f | \
    cutadapt -g CTGTAGGCACCATCAAT -m 0 --discard-untrimmed - 2> ${odir}${fcore}_barcode_5p.report | \
    cutadapt -a AGATCGGAAGAG -m 0 --discard-untrimmed -o ${odir}${fcore}_barcode.fastq - > ${odir}${fcore}_barcode_3p.report
done
