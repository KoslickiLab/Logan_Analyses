#!/bin/bash
# T2T-CHM13v2.0 with Y chromosome
aws s3 --no-sign-request cp \
    s3://human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz \
    ../data
