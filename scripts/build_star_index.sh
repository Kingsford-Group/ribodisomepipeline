#!/bin/bash
nproc=15
cur_dir=`dirname $0`
bin_dir=${cur_dir}/../bin
export PATH=${bin_dir}:$PATH
ref_dir=${cur_dir}/../ref
contaminant_fa=${ref_dir}/rna_coding.fasta
transcript_fa=${ref_dir}/protein_coding_100_filtered.fasta
rrna_idx=${cur_dir}/../StarIndex/contaminant/
transcript_idx=${cur_dir}/../StarIndex/transcript/
mkdir -p ${rrna_idx}
mkdir -p ${transcript_idx}
STAR --runThreadN $nproc --runMode genomeGenerate --genomeDir ${rrna_idx} --genomeFastaFiles ${contaminant_fa} --genomeSAindexNbases 8 --genomeChrBinNbits 11
STAR --runThreadN $nproc --runMode genomeGenerate --genomeDir ${transcript_idx} --genomeFastaFiles ${transcript_fa} --genomeSAindexNbases 11 --genomeChrBinNbits 12

