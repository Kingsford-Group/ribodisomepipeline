disome pipeline
------

1. download data `scripts/get_disome_data.sh`
2. separate barcode and ribo-seq portion from reads `scripts/barcode_seperation.sh`
3. collapse barcode `scripts/barcode_collapse.sh`
4. download STAR `scripts/download_star.sh`
5. download referecences `scripts/build_refs.sh`
6. build STAR indices `scripts/build_star_index.sh`
7. align reads `scripts/batch_align.sh`
8. align reads without barcode collapsing `scripts/batch_align_no_collapse.sh`
9. meta profile analysis per read length on alignments `scripts/meta_analysis.sh`

#### Multi-mapping steps (with ribomap)
1. call ribomap `scripts/ribomap_runner.sh` (for _singlet_ profiles)
2. accumulate read count (for _doublet_ profiles):
  1. under `src` directory, `make read_len_hist_all`
  2. run `src/read_len_hist_all` on doublet library

#### unique-mapping steps (with ribodeblur)
1. accumulate read count (for _doublet_ profiles):
  1. under `src` directory, `make read_len_hist_unique`
  2. call `src/read_len_hist_unique` on doublet library in step 2
2. run `scripts/run_deblur.sh` to get deblur results
3. _notes:_ location of deblur codes is specified under `scripts/deblur_unique_reads.sh`, input bams, output directory are specified under `scripts/file_names.sh`

#### Analysis
1. cdf on doublet positions: `scripts/positional_analysis.py`

#### Notes on analysis settings
* __transcriptome__ (for alignment) all transcripts included (no filter for overlapping genes)
* __read count summary__ all reads used (including multi-mapped reads)