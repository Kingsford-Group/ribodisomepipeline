#!/bin/bash
cur_dir=`dirname $0`
bam_dir=${cur_dir}/../alignment/prelim/
hist_dir=${cur_dir}/../rlen_hist/prelim/
fig_dir=${cur_dir}/../figures/prelim/
src_dir=${cur_dir}/../src/
cds_range=${cur_dir}/../ref/cds_range.txt
mkdir -p ${hist_dir}
mkdir -p ${fig_dir}
# for f in ${bam_dir}*.bam; do
#     fcore=`basename $f`
#     fcore=${fcore%%_transcript_Aligned*}
#     ${src_dir}read_len_hist $f ${hist_dir}${fcore}.hist
# done
for f in ${hist_dir}*.hist; do
    python read_len_hist.py ${cds_range} $f ${fig_dir}rlen_hist/
    #python meta_profile.py ${cds_range} $f ${fig_dir}pos_hist/
done
