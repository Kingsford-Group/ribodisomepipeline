#!/bin/bash
fastq_list=$1
if [ -z "${fastq_list}" ] || [ ! -f "${fastq_list}" ]; then
    echo "usage: sh gerashchenko_analysis_pipeline.sh fastq_list"
    exit
fi
ribomap_dir=/home/hw1/scratch/software_testing/ribomap/
work_dir=`dirname $0`/../
transcript_fa=${work_dir}ref/protein_coding_100_filtered.fasta
contaminant_fa=${work_dir}ref/rna_coding.fasta
cds_range=${work_dir}ref/cds_range.txt
offset=offset_linear.txt
nproc=30
align_dir=${work_dir}/alignment/gerashchenko/
fig_dir=${work_dir}/figures/gerashchenko/codon_usage/
# wrapper to call ribomap
# run_ribomap $riboseq_fq $rnaseq_fq
run_ribomap ()
{
    riboseq_fq=$1
    rnaseq_fq=$2
    rna_core=`basename ${rnaseq_fq}`
    rna_core=${rna_core%%.*}
    ribo_cmd="${ribomap_dir}scripts/run_ribomap.sh --nproc ${nproc} --rnaseq_fq ${rnaseq_fq} --riboseq_fq ${riboseq_fq} --transcript_fa ${transcript_fa} --contaminant_fa ${contaminant_fa} --cds_range ${cds_range} --offset ${offset} --adapter N --work_dir ${work_dir} --sailfish_dir ${work_dir}sm_quant/${rna_core}/ --alignment_dir ${align_dir}"
    ${ribo_cmd} --output_dir ${work_dir}ribomap/gerashchenko/
}

# for riboseq_fq in $(cat ${fastq_list}) ; do
#     # Gerashchenko does not have rna_seq so just gonna pass in a dumb on
#     run_ribomap ${riboseq_fq} ${riboseq_fq}
# done

for f in ${work_dir}ribomap/gerashchenko/*.base; do 
    python single_peak_hist.py ${cds_range} ${transcript_fa} $f ${fig_dir}
done
