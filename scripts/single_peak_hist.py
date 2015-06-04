#!/usr/bin/env python
from io_utils import *
from codon_peak import *
if len(sys.argv) != 5:
    print "Usage: python single_peak_hist cds_range.txt ref.fasta ribomap.base ouput_dir"
    exit(1)

cds_txt = sys.argv[1]
tfasta = sys.argv[2]
pfn = sys.argv[3]
odir = sys.argv[4]

start = 20
stop = 20
ensure_dir(odir)
cds_range = get_cds_range(cds_txt)
tseq = get_tseq(tfasta, cds_range)
peak = get_peaks_from_histfile(pfn, cds_range, start, stop)
single_peak_set_analysis(peak, start, stop, tseq, "", odir+get_file_core(pfn))

