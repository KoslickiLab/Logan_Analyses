#!/bin/bash
sourmash sketch dna -p k=31,scaled=1000 --merge "human_genomes" -o ../data/human_genome_sketches.sig.zip ../*a.gz
