#!/bin/bash
src_dir=`dirname $0`
work_dir=${src_dir}/../
data_dir=${work_dir}data/
sra_dir=${data_dir}sra/
fasta_dir=${data_dir}fasta/

mkdir -p ${fasta_dir}

scp hwang@fozzie.compbio.cs.cmu.edu:/usr1/home/shared/FASTQ/CWRU/150505/* ${fasta_dir}
