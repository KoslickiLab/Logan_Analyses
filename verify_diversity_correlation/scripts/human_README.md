# Quickstart
1. Run `iterate_through_all.sh`. This will iterate over all YACHT coverage values, calling 
   `hash_diversity_correlation.py` which assembles all the YACHT results and distinct hash counting (actually, 
   grabbing from the YACHT output)  
2. (optional) `run hash_diversity_sensitivity.py`: This runs a quick sampling of the same as above, but specifically 
   comapres how things change over the different coverage levels. 
3. then run `get_correlation_figures.sh` which calls `analyze_parquet_data.py` for all the per-metadata plots. This takes a while as the distance corr and MIC take a long time due to O(n^2). 
4. run `plot_custom_correlation.py` for fine-grained plotting by selected metadata variables
5. `analyze_jattr_keys.py` will create a bunch of files showing the most common jattr attributes in the SRA metadata 
   and give you the top ~20 most frequently appearing values. To find interesting metadata variables.

# Notes
For playing with things like: only samples with > N megabases, or > N number of sketches, or > N megabytes in file 
size. This should be done with filters on the `analyze_parquet_data.py` or `plot_custom_correlation.py`. Not all 
filters are implemented in those though as of 12/15/2025