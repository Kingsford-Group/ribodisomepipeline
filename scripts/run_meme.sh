#!/bin/bash
nproc=25
ref_fn=/home/hw1/scratch/ribodisomepipeline/ref/protein_coding_100_filtered.fasta
order=6
bg_fn=/home/hw1/scratch/ribodisomepipeline/ref/transcript.${order}.freq
#fasta-get-markov -m ${order} -dna -norc ${ref_fn} ${bg_fn}

ifa=$1
odir=$2
if [ -z "${ifa}" ] || [ -z "${odir}" ]; then
    echo "usage bash run_meme.sh input.fasta output_dir"
    exit
fi
mkdir -p ${odir}
meme ${ifa} -dna -mod anr -minw 6 -maxw 30 -bfile ${bg_fn} -nmotifs 100 -p ${nproc} -oc ${odir} -evt 1e600 

# for i in `seq 6 30`;
# do
#     mkdir -p ${odir}$i
#     meme ${ifa} -dna -mod anr -w $i -bfile ${bg_fn} -nmotifs 100 -evt 1e600 -p ${nproc} -oc ${odir}$i
# done
# quick test
# meme ${ifa} -dna -mod oops -minw 9 -maxw 12 -nmotifs 10 -bfile ${bg_fn} -p ${nproc} -oc ${odir} 

###### notes
# -text only output plain text without html to stdout
