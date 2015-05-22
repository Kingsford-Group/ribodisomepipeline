#!/bin/bash
src_dir=`dirname $0`
work_dir=${src_dir}/../
ref_dir=${work_dir}ref/
mkdir -p ${ref_dir}
get_data=true
nc_url=ftp://ftp.ensembl.org/pub/release-78/fasta/saccharomyces_cerevisiae/ncrna/Saccharomyces_cerevisiae.R64-1-1.ncrna.fa.gz
trna_url=http://gtrnadb.ucsc.edu/download/tRNAs/eukaryotic-tRNAs.fa.gz
if [ "${get_data}" = true ]; then
    echo "downloading references..."
    wget -P ${ref_dir} -N http://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding.fasta.gz
    wget -P ${ref_dir} -N http://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_genomic_1000.fasta.gz
    wget -P ${ref_dir} -N http://downloads.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_trans.fasta.gz
    wget -P ${ref_dir} -N ${nc_url}
    wget -P ${ref_dir} -N ${trna_url}
    echo "unzipping data..."
    gunzip -f ${ref_dir}*.gz
fi
echo "filtering yeast transcriptome..."
python ${src_dir}/filter_yeast_transcript.py ${ref_dir} 100 #according to Joel
echo "merging rrna and trna..."
rrna_fa=${ref_dir}${nc_url##*/}
rrna_fa=${rrna_fa%.gz}
trna_fa=${ref_dir}${trna_url##*/}
trna_fa=${trna_fa%.gz}
python ${src_dir}/build_contaminant.py ${rrna_fa} ${trna_fa} Saccharomyces_cerevisiae ${ref_dir}yeast_contaminant.fa
