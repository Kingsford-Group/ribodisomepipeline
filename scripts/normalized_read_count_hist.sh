#!/bin/bash
ribomap_odir=`dirname $0`/../ribomap/
for f in ${ribomap_odir}*base; do
    python codon_count_hist.py ../ref/cds_range.txt ../ref/protein_coding_100_filtered.fasta $f ../figures/prelim/codon_usage/
done
