#!/usr/bin/env python
cds_txt = "../ref/cds_range.txt"
# Joel's data
nchx_dfn = "../rlen_hist/prelim/Lib-5-5-15_3_CACGAT_R1_nodup.hist"
nchx_sfn = "../ribomap/Lib-5-5-15_2_AGTTCC_R1_nodup.base"
nchx_nbc_dfn = "../rlen_hist/prelim/Lib-5-5-15_3_CACGAT_R1_nonempty.hist"
nchx_nbc_sfn = "../ribomap/Lib-5-5-15_2_AGTTCC_R1_nonempty.base"
chx_dfn = "../rlen_hist/prelim/Lib-5-5-15_7_ATCACG_R1_nodup.hist"
chx_sfn = "../ribomap/Lib-5-5-15_6_TCCCGA_R1_nodup.base"
chx_nbc_dfn = "../rlen_hist/prelim/Lib-5-5-15_7_ATCACG_R1_nonempty.hist"
chx_nbc_sfn = "../ribomap/Lib-5-5-15_6_TCCCGA_R1_nonempty.base"
# Guydosh's data
wt_dfn = "../rlen_hist/Guydosh14/Guydosh14_wt_disome.hist"
wt_sfn = "../ribomap/Guydosh14/Guydosh14_wt_CHX.base"
dom34_dfn = "../rlen_hist/Guydosh14/Guydosh14_dom34KO_disome.hist"
dom34_sfn = "../ribomap/Guydosh14/Guydosh14_dom34KO_CHX.base"
# doublet filter settings
dlen_min = 57
dlen_max = 62
# peak calling settings
peak_width=3
peak_min_cutoff=3
