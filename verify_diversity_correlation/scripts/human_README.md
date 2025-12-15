# Quickstart
1. Run `iterate_through_all.sh`. This will iterate over all YACHT coverage values, calling 
   `hash_diversity_correlation.py` which assembles all the YACHT results and distinct hash counting (actually, 
   grabbing from the YACHT output)  
2. (optional) `run hash_diversity_sensitivity.py`: This runs a quick sampling of the same as above, but specifically 
   comapres how things change over the different coverage levels. 
3. then `run analyze_parquet_data.py` for all the per-metadata plots. This takes a while as the distance corr and 
   MIC take a long time due to O(n^2). 
4. run `plot_custom_correlation.py` for fine-grained plotting by selected metadata variables