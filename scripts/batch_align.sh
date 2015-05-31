#!/bin/bash
# after barcode collapse
cur_dir=`dirname $0`
data_dir=${cur_dir}/../data/fasta/prelim/
align_dir=${cur_dir}/../alignment/prelim/
mkdir -p ${align_dir}
for f in ${data_dir}*_nodup.fastq.gz; do
    ./read_align.sh $f ${align_dir}
done
