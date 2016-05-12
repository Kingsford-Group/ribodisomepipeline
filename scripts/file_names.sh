#!/bin/bash
work_dir=/home/hw1/scratch/ribodisomepipeline/
joel_align=${work_dir}alignment/prelim/
nchx_dbam=${joel_align}Lib-5-5-15_3_CACGAT_R1_nodup_transcript_Aligned.out.bam
nchx_sbam=${joel_align}Lib-5-5-15_2_AGTTCC_R1_nodup_transcript_Aligned.out.bam
chx_dbam=${joel_align}Lib-5-5-15_7_ATCACG_R1_nodup_transcript_Aligned.out.bam
chx_sbam=${joel_align}Lib-5-5-15_6_TCCCGA_R1_nodup_transcript_Aligned.out.bam
guydosh_align=${work_dir}alignment/Guydosh14/
wt_dbam=${guydosh_align}Guydosh14_wt_disome_transcript_Aligned.out.bam
wt_sbam=${guydosh_align}Guydosh14_wt_CHX_transcript_Aligned.out.bam
dom34_dbam=${guydosh_align}Guydosh14_dom34KO_disome_transcript_Aligned.out.bam
dom34_sbam=${guydosh_align}Guydosh14_dom34KO_CHX_transcript_Aligned.out.bam
rlen_hist_dir=${work_dir}/uniquely_mapped/

