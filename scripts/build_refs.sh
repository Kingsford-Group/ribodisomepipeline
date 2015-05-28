#!/bin/bash
src_dir=`dirname $0`
work_dir=${src_dir}/../
ref_dir=${work_dir}ref/
mkdir -p ${ref_dir}
get_data=true
if [ "${get_data}" = true ]; then
    echo "downloading references..."
    wget -P ${ref_dir} -N http://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding.fasta.gz
    wget -P ${ref_dir} -N http://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_genomic_1000.fasta.gz
    wget -P ${ref_dir} -N http://downloads.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_trans.fasta.gz
    wget -P ${ref_dir} -N http://downloads.yeastgenome.org/sequence/S288C_reference/rna/rna.README
    wget -P ${ref_dir} -N http://downloads.yeastgenome.org/sequence/S288C_reference/rna/rna_coding.fasta.gz
    echo "unzipping data..."
    gunzip -f ${ref_dir}*.gz
fi
echo "filtering yeast transcriptome..."
python ${src_dir}/filter_yeast_transcript.py ${ref_dir} 100 #according to Joel
