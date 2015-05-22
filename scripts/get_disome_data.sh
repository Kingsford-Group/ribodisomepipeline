#!/bin/bash
src_dir=`dirname $0`
odir=${src_dir}/../data/raw/
mkdir -p ${odir}
scp hwang@fozzie.compbio.cs.cmu.edu:/usr1/home/shared/FASTQ/CWRU/150505/* ${odir}
