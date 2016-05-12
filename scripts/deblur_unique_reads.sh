#!/bin/bash
deblur_py_dir=/home/hw1/scratch/ribodeblur/scripts/
work_dir=`dirname $0`/../
rlen_hist=${work_dir}src/read_len_hist_unique
cds_range=${work_dir}ref/cds_range.txt

#-----------------------------
# generate rlen_hist file
# get_rlen_hist $bam $odir
get_rlen_hist()
{
    bam=$1
    rlen_odir=$2
    mkdir -p ${rlen_odir}
    file_core=`basename $bam`
    file_core=${file_core%%_transcript_Aligned.out.bam}
    echo "${file_core}"
    hist_fn=${rlen_odir}${file_core}.hist
    ${rlen_hist} $bam ${hist_fn}
}

#-----------------------------
# deblur noisy read pileups
# deblur_reads $bam $output_dir
deblur_reads ()
{
    bam=$1
    deblur_odir=$2
    file_core=`basename $bam`
    file_core=${file_core%%_transcript_Aligned.out.bam}
    hist_fn=${deblur_odir}${file_core}.hist
    vblur_fn=${deblur_odir}${file_core}.vblur
    echo "${file_core}"
    mkdir -p ${deblur_odir}
    # step 1: read length grouping
    echo "making rlen.hist..."
    ${rlen_hist} $bam ${hist_fn}
    # # step 2: meta analysis
    # echo "plotting meta profiles..."
    # python ${deblur_py_dir}meta_profile.py ${hist_fn} ${cds_range} ${deblur_odir}
    # step 3: deblur
    python ${deblur_py_dir}train_vblur_from_meta.py ${hist_fn} ${cds_range} ${deblur_odir}
    time python ${deblur_py_dir}recover_asite_profile.py ${hist_fn} ${vblur_fn} ${cds_range} ${deblur_odir}
    # python ${deblur_py_dir}profile_evaluation.py ${hist_fn} ${cds_range} ${vblur_fn} ${deblur_odir}${fq_core}.eps ${deblur_odir}
}
