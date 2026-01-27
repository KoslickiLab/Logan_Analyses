#!/bin/bash
aws s3 --no-sign-request ls s3://human-pangenomics/working/HPRC/ --recursive \
    | grep release2 \
    | grep "\.fa\.gz$" \
    | awk '{print "s3://human-pangenomics/" $NF}' \
    | xargs -P 8 -I {} aws s3 --no-sign-request cp {} ../data
