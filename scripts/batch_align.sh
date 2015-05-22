#!/bin/bash
for f in ../data/fasta/*.fastq.gz; do
    ./read_align.sh $f
done
