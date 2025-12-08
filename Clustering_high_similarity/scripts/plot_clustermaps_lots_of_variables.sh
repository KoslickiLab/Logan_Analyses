#!/bin/bash
while read value; do
    python annotate_clustermap.py --linkage ../data/full_label_prop_viz_50k/expanded_linkage_matrix.npy --similarity ../data/full_label_prop_viz_50k/expanded_similarity_matrix.npz --sample-ids ../data/full_label_prop_viz_50k/expanded_sample_ids.txt  --database /scratch/shared_data_new/Logan_yacht_data/metadata/aws_sra_metadata/metadata_geo_joined.duckdb --metadata ${value}:50 -o ../data/full_label_prop_viz_50k/clustermap_${value}.png --figsize 50 37 --cmap viridis --dpi 400 --subsample 10000
done < metadata_cat_values.txt
