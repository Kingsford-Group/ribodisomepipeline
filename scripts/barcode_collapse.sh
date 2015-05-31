#!/bin/bash
cur_dir=`dirname $0`
data_dir=${cur_dir}/../data/fasta/prelim/
src_dir=${cur_dir}/../src/
for f in ${data_dir}*_prefilter.fastq.gz; do
    fcore=`basename $f`
    fcore=${fcore%%_prefilter*}
    ${src_dir}collapse_barcode $f ${data_dir}${fcore}_barcode.fastq.gz ${data_dir}${fcore}_nodup.fastq.gz
done
