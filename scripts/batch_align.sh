#!/bin/bash
# after barcode collapse
for f in ../data/fasta/*_nodup.fastq.gz; do
    ./read_align.sh $f
done
