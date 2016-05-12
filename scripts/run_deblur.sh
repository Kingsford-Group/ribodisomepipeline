#!/bin/bash
source file_names.sh
source deblur_unique_reads.sh

get_rlen_hist ${nchx_dbam} ${rlen_hist_dir}
get_rlen_hist ${chx_dbam} ${rlen_hist_dir}
get_rlen_hist ${wt_dbam} ${rlen_hist_dir}
get_rlen_hist ${dom34_dbam} ${rlen_hist_dir}

deblur_reads ${nchx_sbam} ${rlen_hist_dir}
deblur_reads ${chx_sbam} ${rlen_hist_dir}
deblur_reads ${wt_sbam} ${rlen_hist_dir}
deblur_reads ${dom34_sbam} ${rlen_hist_dir}
