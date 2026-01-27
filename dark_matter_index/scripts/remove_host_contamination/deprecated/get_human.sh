#!/bin/bash
#grep -F -w 9606 *.accession2taxid | cut -d' ' -f 1
grep -F -w 9606 nucl_gb.accession2taxid | sed 's/\t.*//g' > human_accessions.txt
grep -F -w 9606 nucl_wgs.accession2taxid | sed 's/\t.*//g' >> human_accessions.txt
grep -F -w 9606 nucl_wgs_EXTRA.accession2taxid | sed 's/\t.*//g' >> human_accessions.txt
grep -F -w 9606 wgs.accession2taxid | sed 's/\t.*//g' >> human_accessions.txt
