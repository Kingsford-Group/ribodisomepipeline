#!/bin/bash
# pre-compiled STAR indices should be provided otherwise STAR won't run
riboseq_fq=$1
if [ -z "${riboseq_fq}" ]; then
    echo "Usage: ./read_align.sh ribo-seq.fq.gz"
    exit
fi
#=============================
# default parameters
#=============================
cur_dir=`dirname $0`
bin_dir=${cur_dir}/../bin
export PATH=${bin_dir}:$PATH
align_dir=${cur_dir}/../data/alignments/
mkdir -p ${align_dir}
contaminant_idx=/home/hw1/scratch/scratch2/ribomap-playground/yeast/StarIndex/contaminant/
transcript_idx=/home/hw1/scratch/scratch2/ribomap-playground/yeast/StarIndex/transcript/
adapter=CTGTAGGCACCATCAAT
nproc=30
#============================================
# step 1: filter rrna
#============================================
echo "filtering contaminated sequences in riboseq"
ribo_core=`basename ${riboseq_fq}`
ribo_core=${ribo_core%%.*}
common_params="--runThreadN ${nproc} --seedSearchLmax 10 --outFilterMultimapScoreRange 0 --outFilterMultimapNmax 255 --outFilterMismatchNmax 1 --clip3pAdapterSeq ${adapter} --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterIntronMotifs RemoveNoncanonical"
oriboprefix=${align_dir}${ribo_core}_rrna_
ribo_nrrna_fa=${oriboprefix}Unmapped.out.mate1
STAR --genomeDir ${contaminant_idx} --readFilesIn ${riboseq_fq} --outFileNamePrefix ${oriboprefix} --outStd SAM --outReadsUnmapped Fastx --outSAMmode NoQS "--readFilesCommand zcat <" ${common_params} > /dev/null
#============================================
# step 2: map to transcriptome
#============================================
echo "aligning riboseq to the transcriptome"
oriboprefix=${align_dir}${ribo_core}_transcript_
STAR --genomeDir ${transcript_idx} --readFilesIn ${ribo_nrrna_fa} --outFileNamePrefix ${oriboprefix} --outSAMtype BAM Unsorted --outSAMmode NoQS --outSAMattributes NH NM ${common_params}
