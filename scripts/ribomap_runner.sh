#!/bin/bash
ribomap_dir=/home/hw1/scratch/software_testing/ribomap/
work_dir=`dirname $0`/../
transcript_fa=${work_dir}ref/protein_coding_100_filtered.fasta
contaminant_fa=${work_dir}ref/rna_coding.fasta
cds_range=${work_dir}ref/cds_range.txt
offset=offset.txt
nproc=30
align_dir=${work_dir}/alignment/prelim/
data_dir=${work_dir}data/fasta/prelim/
# wrapper to call ribomap
# run_ribomap $riboseq_fq $rnaseq_fq
run_ribomap ()
{
    riboseq_fq=$1
    rnaseq_fq=$2
    rna_core=`basename ${rnaseq_fq}`
    rna_core=${rna_core%%.*}
    ribo_cmd="${ribomap_dir}scripts/run_ribomap.sh --nproc ${nproc} --rnaseq_fq ${rnaseq_fq} --riboseq_fq ${riboseq_fq} --transcript_fa ${transcript_fa} --contaminant_fa ${contaminant_fa} --cds_range ${cds_range} --offset ${offset} --adapter N --work_dir ${work_dir} --sailfish_dir ${work_dir}sm_quant/${rna_core}/ --alignment_dir ${align_dir}"
    ${ribo_cmd} --output_dir ${work_dir}ribomap
}

#==============================
# no chx
# before barcode collapse
#==============================
riboseq_fq=${data_dir}Lib-5-5-15_2_AGTTCC_R1_nonempty.fastq.gz
rnaseq_fq=${data_dir}Lib-5-5-15_4_GTAGAG_R1_nonempty.fastq.gz
run_ribomap ${riboseq_fq} ${rnaseq_fq}
#==============================
# no chx
# after barcode collapse
#==============================
riboseq_fq=${data_dir}Lib-5-5-15_2_AGTTCC_R1_nodup.fastq.gz
rnaseq_fq=${data_dir}Lib-5-5-15_4_GTAGAG_R1_nodup.fastq.gz
run_ribomap ${riboseq_fq} ${rnaseq_fq}
#==============================
# 10 chx
# before barcode collapse
#==============================
riboseq_fq=${data_dir}Lib-5-5-15_6_TCCCGA_R1_nonempty.fastq.gz
rnaseq_fq=${data_dir}Lib-5-5-15_8_CGATGT_R1_nonempty.fastq.gz
run_ribomap ${riboseq_fq} ${rnaseq_fq}
#==============================
# 10 chx
# after barcode collapse
#==============================
riboseq_fq=${data_dir}Lib-5-5-15_6_TCCCGA_R1_nodup.fastq.gz
rnaseq_fq=${data_dir}Lib-5-5-15_8_CGATGT_R1_nodup.fastq.gz
run_ribomap ${riboseq_fq} ${rnaseq_fq}
