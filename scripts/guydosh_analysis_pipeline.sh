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
mkdir -p ${sra_dir}
mkdir -p ${fasta_dir}
mkdir -p ${trim_dir}
mkdir -p ${hist_dir}
mkdir -p ${fig_dir}
#=============================
# step 1: download sra
#=============================
echo "downloading data..."
# wild-type mRNA-Seq
wget -P ${sra_dir} -N ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX386%2FSRX386982/SRR1042851/SRR1042851.sra
# dom34KO mRNA-Seq
wget -P ${sra_dir} -N ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX386%2FSRX386983/SRR1042852/SRR1042852.sra
# wild-type CHX
wget -P ${sra_dir} -N ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX386%2FSRX386984/SRR1042853/SRR1042853.sra
# dom34KO CHX
wget -P ${sra_dir} -N ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX386%2FSRX386985/SRR1042854/SRR1042854.sra 
# wild-type short footprints
wget -P ${sra_dir} -N ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX387%2FSRX387003/SRR1042876/SRR1042876.sra
# dom34KO short footprints
wget -P ${sra_dir} -N ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX387%2FSRX387004/SRR1042877/SRR1042877.sra
# wild-type disome footprints
wget -P ${sra_dir} -N ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX387%2FSRX387005/SRR1042878/SRR1042878.sra
wget -P ${sra_dir} -N ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX387%2FSRX387005/SRR1042879/SRR1042879.sra
# dom34KO disome footprints
wget -P ${sra_dir} -N ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX387%2FSRX387006/SRR1042880/SRR1042880.sra
wget -P ${sra_dir} -N ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX387%2FSRX387006/SRR1042881/SRR1042881.sra
#==============================
# step 2: convert sra to fastq
#==============================
echo "converting sra to fastq..."
for f in ${sra_dir}*.sra; do
    fastq-dump $f -O ${fasta_dir} --fasta --gzip
done
#==============================
# step 3: rename files
#==============================
echo "renaming fastq files to be more informative..."
mv ${fasta_dir}SRR1042851.fasta.gz ${fasta_dir}WT_RNAseq.fasta.gz
mv ${fasta_dir}SRR1042852.fasta.gz ${fasta_dir}Dom34KO_RNAseq.fasta.gz
mv ${fasta_dir}SRR1042853.fasta.gz ${fasta_dir}WT-Chx.fasta.gz
mv ${fasta_dir}SRR1042854.fasta.gz ${fasta_dir}Dom34KO_Chx.fasta.gz
mv ${fasta_dir}SRR1042876.fasta.gz ${fasta_dir}WT_short.fasta.gz
mv ${fasta_dir}SRR1042877.fasta.gz ${fasta_dir}Dom34KO_short.fasta.gz
mv ${fasta_dir}SRR1042878.fasta.gz ${fasta_dir}WT_disome1.fasta.gz
mv ${fasta_dir}SRR1042879.fasta.gz ${fasta_dir}WT_disome2.fasta.gz
mv ${fasta_dir}SRR1042880.fasta.gz ${fasta_dir}Dom34KO_disome1.fasta.gz
mv ${fasta_dir}SRR1042881.fasta.gz ${fasta_dir}Dom34KO_disome2.fasta.gz
#==============================
# step 4: trim adapter
#==============================
cutadapt -u -16 -o ${trim_dir}WT_RNAseq.fasta.gz ${fasta_dir}WT_RNAseq.fasta.gz  
cutadapt -u -16 -o ${trim_dir}Dom34KO_RNAseq.fasta.gz ${fasta_dir}Dom34KO_RNAseq.fasta.gz 
cutadapt -a CTGTAGGCACCATCAAT -m 1 -o ${trim_dir}WT-Chx.fasta.gz ${fasta_dir}WT-Chx.fasta.gz
cutadapt -a CTGTAGGCACCATCAAT -m 1 -o ${trim_dir}Dom34KO_Chx.fasta.gz ${fasta_dir}Dom34KO_Chx.fasta.gz
cutadapt -a CTGTAGGCACCATCAAT -m 1 -o ${trim_dir}WT_short.fasta.gz ${fasta_dir}WT_short.fasta.gz
cutadapt -a CTGTAGGCACCATCAAT -m 1 -o ${trim_dir}Dom34KO_short.fasta.gz ${fasta_dir}Dom34KO_short.fasta.gz 
cutadapt -a CTGTAGGCACCATCAAT -m 1 ${fasta_dir}WT_disome1.fasta.gz > ${trim_dir}WT_disome.fasta
cutadapt -a CTGTAGGCACCATCAAT -m 1 ${fasta_dir}WT_disome2.fasta.gz >> ${trim_dir}WT_disome.fasta
gzip -f ${trim_dir}WT_disome.fasta
cutadapt -a CTGTAGGCACCATCAAT -m 1 ${fasta_dir}Dom34KO_disome1.fasta.gz > ${trim_dir}Dom34KO_disome.fasta
cutadapt -a CTGTAGGCACCATCAAT -m 1 ${fasta_dir}Dom34KO_disome2.fasta.gz >> ${trim_dir}Dom34KO_disome.fasta
gzip -f ${trim_dir}Dom34KO_disome.fasta
#==============================
# step 5: star align
#==============================
for f in ${trim_dir}*.gz ; do
    ./read_align.sh $f ${align_dir}
done
#==============================
# step 6: read len hist
#==============================
for f in ${align_dir}*.bam; do
    fcore=`basename $f`
    fcore=${fcore%%_transcript_Aligned*}
    ${src_dir}read_len_hist $f ${hist_dir}${fcore}.hist
done
for f in ${hist_dir}*.hist; do
    python read_len_hist.py ${cds_range} $f ${fig_dir}
    python meta_profile.py ${cds_range} $f ${fig_dir}
done
