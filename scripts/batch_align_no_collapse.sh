#!/bin/bash
cur_dir=`dirname $0`
data_dir=${cur_dir}/../data/fasta/prelim/
src_dir=${cur_dir}/../src/
align_dir=${cur_dir}/../alignment/prelim/
mkdir -p ${align_dir}
filter out empty reads
for f in ${data_dir}*_prefilter.fastq.gz; do
    fcore=`basename $f`
    fcore=${fcore%%_prefilter*}
    ${src_dir}filter_zero_seq $f ${data_dir}${fcore}_nonempty.fastq.gz
done
# align reads with STAR
for f in ${data_dir}*_nonempty.fastq.gz; do
    ./read_align.sh $f ${align_dir}
done
