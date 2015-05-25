disome pipeline
------

1. download data `scripts/get_disome_data.sh`
2. separate barcode and ribo-seq portion from reads `scripts/barcode_seperation.sh`
3. collapse barcode `src/collapse_barcode read_in.fq.gz barcode.fq.gz read_out.fq.gz`
4. download STAR `scripts/download_star.sh`
5. download referecences `scripts/build_refs.sh`
6. build STAR indices `scripts/build_star_index.sh`
7. align reads `scripts/batch_align.sh`