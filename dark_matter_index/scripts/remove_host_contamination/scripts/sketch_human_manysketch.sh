#!/bin/bash
printf "name,genome_filename,protein_filename\n" > manysketch.csv
for file in `ls ../data/*a.gz`; do
	#printf "${file},${file},\n" >> manysketch.csv
	base=$(basename "${file}")
        abspath=$(realpath "${file}")
        printf "${base},${abspath},\n" >> manysketch.csv
done

sourmash scripts manysketch -o ../data/human_genome_sketches_all.sig.zip -p dna,k=31,scaled=1000,noabund -c 200 -f manysketch.csv
