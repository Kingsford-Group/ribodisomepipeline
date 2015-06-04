#!/usr/bin/env python
from codon_peak import *
f __name__ == "__main__":
    cds_txt = "../ref/cds_range.txt"
    tfasta = "../ref/protein_coding_100_filtered.fasta"
    p_nc_nb_fn = "../ribomap/Lib-5-5-15_2_AGTTCC_R1_nodup.base"
    p_nc_wb_fn = "../ribomap/Lib-5-5-15_2_AGTTCC_R1_nonempty.base"
    p_wc_nb_fn = "../ribomap/Lib-5-5-15_6_TCCCGA_R1_nodup.base"
    p_wc_wb_fn = "../ribomap/Lib-5-5-15_6_TCCCGA_R1_nonempty.base"
    odir = "../figures/prelim/codon_usage/"
    start=20
    stop=20

    cds_range = get_cds_range(cds_txt)
    tseq = get_tseq(tfasta, cds_range)
    # no Chx collapse barcode
    p_ncnb = get_peaks_from_histfile(p_nc_nb_fn, cds_range, start, stop)
    # no Chx no collapse barcode
    p_ncwb = get_peaks_from_histfile(p_nc_wb_fn, cds_range, start, stop)
    # Chx collapse barcode
    p_wcnb = get_peaks_from_histfile(p_wc_nb_fn, cds_range, start, stop)
    # Chx no collapse barcode
    p_wcwb = get_peaks_from_histfile(p_wc_wb_fn, cds_range, start, stop)
    print "no chx barcode comparison 1: collapsed barcode 2:no collapse"
    compare_two_peak_sets(p_ncnb, p_ncwb, start, stop, tseq, 'no Chx', odir+'noChxYNbarcode')
    print "10 chx barcode comarison 1: collapsed barcode 2: no collpase"
    compare_two_peak_sets(p_wcnb, p_wcwb, start, stop, tseq, '10x Chx', odir+'Chx10YNbarcode')
    print "chx comarison (barcode collapsed) 1: no chx 2: 10 chx"
    compare_two_peak_sets(p_ncnb, p_wcnb, start, stop, tseq, 'barcode collapsed', odir+'BarcodeYNChx')
