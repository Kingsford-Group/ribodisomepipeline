disome pipeline
------

1. download data `./get_disome_data.sh`
2. clip off Illumina-specific library sequence before alignment `./clip_illumina_seq.sh`
3. download STAR `./download_star.sh`
4. download referecences `./build_refs.sh`
5. build STAR indices `./build_star_index.sh`
6. align reads `./batch_align.sh`