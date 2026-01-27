#!/bin/bash
#set -euo pipefail

# HPRC Y1


# Get the index
wget -q https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Assemblies/main/assembly_index/Year1_assemblies_v2_genbank.index 

# Download all assemblies (maternal and paternal)
# Skip header lines (if any) and download both columns 2 and 3
tail -n +1 Year1_assemblies_v2_genbank.index | while read line; do
    maternal=$(echo "$line" | awk '{print $2}')
    paternal=$(echo "$line" | awk '{print $3}')
    
    if [[ -n "$maternal" ]]; then
        echo "Downloading: $maternal"
        aws s3 --no-sign-request cp ${maternal} ../data &
    fi
    if [[ -n "$paternal" ]]; then
        echo "Downloading: $paternal"
        aws s3 --no-sign-request cp ${paternal} ../data &
    fi
    
    # Limit parallel downloads
    while [ $(jobs -r | wc -l) -ge 8 ]; do sleep 1; done
done

wait
echo "Done! Downloaded $(ls *.fa.gz | wc -l) assemblies"
