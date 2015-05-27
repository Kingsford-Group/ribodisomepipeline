#!/bin/bash
cur_dir=`dirname $0`
idir=${cur_dir}/../data/raw/
odir=${cur_dir}/../data/fasta/
mkdir -p ${odir}
for f in ${idir}*.fastq.gz; do
    fname=`basename $f`
    zcat $f | \
    fastx_clipper -Q33 -a AGATCGGAAGAG -l 40 -c -n –v -z -o ${odir}${fname}
done
