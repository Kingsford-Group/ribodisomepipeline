#!/bin/bash
cur_dir=`dirname $0`
work_dir=${cur_dir}/../
data_dir=${work_dir}data/
sra_dir=${data_dir}sra/Guydosh14/
fasta_dir=${data_dir}raw/Guydosh14/
trim_dir=${data_dir}fasta/Guydosh14/
align_dir=${work_dir}/alignment/Guydosh14/
hist_dir=${work_dir}rlen_hist/Guydosh14/
fig_dir=${work_dir}figures/Guydosh14/
src_dir=${work_dir}src/
cds_range=${work_dir}ref/cds_range.txt

ribomap_dir=/home/hw1/scratch/software_testing/ribomap/
transcript_fa=${work_dir}ref/protein_coding_100_filtered.fasta
contaminant_fa=${work_dir}ref/rna_coding.fasta
cds_range=${work_dir}ref/cds_range.txt
offset=offset_linear.txt
nproc=30

# mkdir -p ${sra_dir}
# mkdir -p ${fasta_dir}
# mkdir -p ${trim_dir}
# mkdir -p ${hist_dir}
# mkdir -p ${fig_dir}
# #=============================
# # step 1: download sra
# #=============================
# echo "downloading data..."
# # wild-type mRNA-Seq
# wget -P ${sra_dir} -N ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX386%2FSRX386982/SRR1042851/SRR1042851.sra
# # dom34KO mRNA-Seq
# wget -P ${sra_dir} -N ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX386%2FSRX386983/SRR1042852/SRR1042852.sra
# # wild-type CHX
# wget -P ${sra_dir} -N ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX386%2FSRX386984/SRR1042853/SRR1042853.sra
# # dom34KO CHX
# wget -P ${sra_dir} -N ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX386%2FSRX386985/SRR1042854/SRR1042854.sra 
# # wild-type short footprints
# wget -P ${sra_dir} -N ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX387%2FSRX387003/SRR1042876/SRR1042876.sra
# # dom34KO short footprints
# wget -P ${sra_dir} -N ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX387%2FSRX387004/SRR1042877/SRR1042877.sra
# # wild-type disome footprints
# wget -P ${sra_dir} -N ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX387%2FSRX387005/SRR1042878/SRR1042878.sra
# wget -P ${sra_dir} -N ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX387%2FSRX387005/SRR1042879/SRR1042879.sra
# # dom34KO disome footprints
# wget -P ${sra_dir} -N ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX387%2FSRX387006/SRR1042880/SRR1042880.sra
# wget -P ${sra_dir} -N ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX387%2FSRX387006/SRR1042881/SRR1042881.sra
# #==============================
# # step 2: convert sra to fastq
# #==============================
# echo "converting sra to fastq..."
# for f in ${sra_dir}*.sra; do
#     fastq-dump $f -O ${fasta_dir} --fasta --gzip
# done
# #==============================
# # step 3: rename files
# #==============================
# echo "renaming fastq files to be more informative..."
# mv ${fasta_dir}SRR1042851.fasta.gz ${fasta_dir}Guydosh14_wt_mRNA.fasta.gz
# mv ${fasta_dir}SRR1042852.fasta.gz ${fasta_dir}Guydosh14_dom34KO_mRNA.fasta.gz
# mv ${fasta_dir}SRR1042853.fasta.gz ${fasta_dir}Guydosh14_wt_CHX.fasta.gz
# mv ${fasta_dir}SRR1042854.fasta.gz ${fasta_dir}Guydosh14_dom34KO_CHX.fasta.gz
# mv ${fasta_dir}SRR1042876.fasta.gz ${fasta_dir}Guydosh14_wt_short.fasta.gz
# mv ${fasta_dir}SRR1042877.fasta.gz ${fasta_dir}Guydosh14_dom34KO_short.fasta.gz
# mv ${fasta_dir}SRR1042878.fasta.gz ${fasta_dir}Guydosh14_wt_disome1.fasta.gz
# mv ${fasta_dir}SRR1042879.fasta.gz ${fasta_dir}Guydosh14_wt_disome2.fasta.gz
# mv ${fasta_dir}SRR1042880.fasta.gz ${fasta_dir}Guydosh14_dom34KO_disome1.fasta.gz
# mv ${fasta_dir}SRR1042881.fasta.gz ${fasta_dir}Guydosh14_dom34KO_disome2.fasta.gz
# #==============================
# # step 4: trim adapter
# #==============================
# # following original paper setting:
# # RNA-seq reads: trim to 35nt
# # -u -16: trim 16 bases from 3' ends of reads
# # http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1279566
# cutadapt -u -16 -o ${trim_dir}Guydosh14_wt_mRNA.fasta.gz ${fasta_dir}Guydosh14_wt_mRNA.fasta.gz  
# cutadapt -u -16 -o ${trim_dir}Guydosh14_dom34KO_mRNA.fasta.gz ${fasta_dir}Guydosh14_dom34KO_mRNA.fasta.gz 
# # ribo-seq reads: 
# # -a CTGTAGGCACCATCAAT: trim universal adapter from 3' ends of reads
# # -m 1: filter out empty length reads after trimming (min_len = 1)
# cutadapt -a CTGTAGGCACCATCAAT -m 1 -o ${trim_dir}Guydosh14_wt_CHX.fasta.gz ${fasta_dir}Guydosh14_wt_CHX.fasta.gz
# cutadapt -a CTGTAGGCACCATCAAT -m 1 -o ${trim_dir}Guydosh14_dom34KO_CHX.fasta.gz ${fasta_dir}Guydosh14_dom34KO_CHX.fasta.gz
# cutadapt -a CTGTAGGCACCATCAAT -m 1 -o ${trim_dir}Guydosh14_wt_short.fasta.gz ${fasta_dir}Guydosh14_wt_short.fasta.gz
# cutadapt -a CTGTAGGCACCATCAAT -m 1 -o ${trim_dir}Guydosh14_dom34KO_short.fasta.gz ${fasta_dir}Guydosh14_dom34KO_short.fasta.gz 
# cutadapt -a CTGTAGGCACCATCAAT -m 1 ${fasta_dir}Guydosh14_wt_disome1.fasta.gz > ${trim_dir}Guydosh14_wt_disome.fasta
# cutadapt -a CTGTAGGCACCATCAAT -m 1 ${fasta_dir}Guydosh14_wt_disome2.fasta.gz >> ${trim_dir}Guydosh14_wt_disome.fasta
# gzip -f ${trim_dir}Guydosh14_wt_disome.fasta
# cutadapt -a CTGTAGGCACCATCAAT -m 1 ${fasta_dir}Guydosh14_dom34KO_disome1.fasta.gz > ${trim_dir}Guydosh14_dom34KO_disome.fasta
# cutadapt -a CTGTAGGCACCATCAAT -m 1 ${fasta_dir}Guydosh14_dom34KO_disome2.fasta.gz >> ${trim_dir}Guydosh14_dom34KO_disome.fasta
# gzip -f ${trim_dir}Guydosh14_dom34KO_disome.fasta
# #==============================
# # step 5: star align
# #==============================
# for f in ${trim_dir}*.gz ; do
#     ./read_align.sh $f ${align_dir}
# done
# #==============================
# # step 6: read len hist
# #==============================
# for f in ${align_dir}*.bam; do
#     fcore=`basename $f`
#     fcore=${fcore%%_transcript_Aligned*}
#     ${src_dir}read_len_hist $f ${hist_dir}${fcore}.hist
# done
# for f in ${hist_dir}*.hist; do
#     python read_len_hist.py ${cds_range} $f ${fig_dir}
#     python meta_profile.py ${cds_range} $f ${fig_dir}
# done

#==============================
# step 7: ribomap result
#==============================
# wrapper to call ribomap
# run_ribomap $riboseq_fq $rnaseq_fq
run_ribomap ()
{
    riboseq_fq=$1
    rnaseq_fq=$2
    rna_core=`basename ${rnaseq_fq}`
    rna_core=${rna_core%%.*}
    ribo_cmd="${ribomap_dir}scripts/run_ribomap.sh --nproc ${nproc} --rnaseq_fq ${rnaseq_fq} --riboseq_fq ${riboseq_fq} --transcript_fa ${transcript_fa} --contaminant_fa ${contaminant_fa} --cds_range ${cds_range} --offset ${offset} --adapter N --work_dir ${work_dir} --sailfish_dir ${work_dir}sm_quant/${rna_core}/ --alignment_dir ${align_dir}"
    ${ribo_cmd} --output_dir ${work_dir}ribomap/Guydosh14/
}

#==============================
# WT
#==============================
riboseq_fq=${trim_dir}Guydosh14_wt_CHX.fasta.gz
rnaseq_fq=${trim_dir}Guydosh14_wt_mRNA.fasta.gz
run_ribomap ${riboseq_fq} ${riboseq_fq}
#==============================
# Dom34KO
#==============================
riboseq_fq=${trim_dir}Guydosh14_dom34KO_CHX.fasta.gz
rnaseq_fq=${trim_dir}Guydosh14_dom34KO_mRNA.fasta.gz
run_ribomap ${riboseq_fq} ${riboseq_fq}

# for f in ${work_dir}ribomap/gerashchenko/*.base; do
#     python single_peak_hist.py ${cds_range} ${transcript_fa} $f ${fig_dir}
# done


# IMPORTANT NOTE:
# bam names are simply updated to a new name version
# DID NOT rerun pipeline until Salmon quant (v0.6.0)
# because current version of STAR won't run on old version of param settings
# ERROR MSG:
# Mar 23 15:36:07 ..... Started STAR run

# EXITING: FATAL INPUT ERROR: unrecoginzed parameter name "sjdbInsertSave" in input "genomeParameters.txt"
# SOLUTION: use correct parameter name (check the manual)
#
# this means STAR index needs to be regenerated
# quick fix is skip STAR align for now
