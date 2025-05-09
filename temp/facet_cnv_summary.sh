#!/bin/bash

echo -e "sample,purity,ploidy,dipLogR,est_insert_size" > facet_cnv_summary.csv

sample_count=0

find /home/zhonggr/projects/250224_DFSP_WES/data/wes/variant_calling/cnv/facets -name "*.vcf.gz" | while read vcf; do
    dir_name=$(basename "$(dirname "$vcf")")
    
    sample=$(echo "$dir_name" | sed 's/_vs_.*//')

    purity=$(zcat "$vcf" | grep "^##purity=" | cut -d= -f2)
    ploidy=$(zcat "$vcf" | grep "^##ploidy=" | cut -d= -f2)
    dipLogR=$(zcat "$vcf" | grep "^##dipLogR=" | cut -d= -f2)
    est_insert_size=$(zcat "$vcf" | grep "^##est_insert_size=" | cut -d= -f2)
    echo -e "${sample},${purity},${ploidy},${dipLogR},${est_insert_size}" >> facet_cnv_summary.csv

    sample_count=$((sample_count + 1))
done

# Count lines in CSV (minus header) to get sample count
sample_count=$(wc -l < facet_cnv_summary.csv)
sample_count=$((sample_count - 1))  # Subtract 1 for the header line

# Print the total count
echo "Total samples processed: $sample_count"

