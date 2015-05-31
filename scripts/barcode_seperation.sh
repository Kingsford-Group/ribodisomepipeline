#!/bin/bash
cur_dir=`dirname $0`
idir=${cur_dir}/../data/raw/
odir=${cur_dir}/../data/fasta/prelim/
mkdir -p ${odir}
for f in ${idir}*.fastq.gz; do
    fcore=`basename $f`
    fcore=${fcore%%.*}
    cutadapt -a CTGTAGGCACCATCAAT -o ${odir}${fcore}_prefilter.fastq.gz $f > ${odir}${fcore}.report &
    cutadapt -g CTGTAGGCACCATCAAT $f 2> ${odir}${fcore}_barcode_5p.report | \
    cutadapt -a AGATCGGAAGAG -o ${odir}${fcore}_barcode.fastq.gz - > ${odir}${fcore}_barcode_3p.report &
done
