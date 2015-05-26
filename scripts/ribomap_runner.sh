#!/bin/bash
ribomap_dir=/home/hw1/scratch/software_testing/ribomap/
work_dir=`dirname $0`/../
transcript_fa=${work_dir}ref/protein_coding_100_filtered.fasta
contaminant_fa=${work_dir}ref/yeast_contaminant.fa
cds_range=${work_dir}ref/cds_range.txt
offset=offset.txt
nproc=30
#==============================
# before barcode collapse
#==============================
riboseq_fq=${work_dir}data/fasta/Lib-5-5-15_2_AGTTCC_R1_nonempty.fastq.gz
rnaseq_fq=${work_dir}data/fasta/Lib-5-5-15_4_GTAGAG_R1_nonempty.fastq.gz
ribo_cmd="${ribomap_dir}scripts/run_ribomap.sh --nproc ${nproc} --rnaseq_fq ${rnaseq_fq} --riboseq_fq ${riboseq_fq} --transcript_fa ${transcript_fa} --contaminant_fa ${contaminant_fa} --cds_range ${cds_range} --offset ${offset} --adapter N --work_dir ${work_dir} --sailfish_dir ${work_dir}sm_quant/Lib-5-5-15_4_GTAGAG_R1_nonempty/"
${ribo_cmd} --output_dir ${work_dir}ribomap
#==============================
# after barcode collapse
#==============================
