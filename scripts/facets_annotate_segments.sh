#!/bin/bash

## Use source to activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
if ! conda activate varcall; then
    echo "ERROR: Failed to activate conda environment 'varcall'"
    exit 1
fi

## Setup the working directory
facet_dir=/home/zhonggr/projects/250224_DFSP_WES/data/wes/variant_calling/cnv/facets
# facet_dir=/home/zhonggr/projects/250224_DFSP_WES/data/wes/variant_calling/cnv/facets_test_normal

annotation=/home/zhonggr/projects/250224_DFSP_WES/data/reference/Gencode/annotation_protein_coding.bed

samples=${facet_dir}/samples.txt

find ${facet_dir} -maxdepth 1 -mindepth 1 -type d -exec basename {} \; > ${samples}

for tumour_id in $(cat ${samples}); do

    echo $(date +"%F") $(date +"%T") "Processing : ${tumour_id}"
    
    # Create the output file path
    out_dir=${facet_dir}/${tumour_id}

    ## Get the header from the segment file
    seg_file=${out_dir}/${tumour_id}.seg
    seg_header_file=${out_dir}/${tumour_id}.seg_header.txt
    head -n 1 ${seg_file} > ${seg_header_file}
    
    ## Create the header for the annotation file
    annotation_header_file=${out_dir}/${tumour_id}.annotation_header.txt
    echo -e "Chr\tStart\tEnd\tGene" > ${annotation_header_file}
    
    ## Combine the headers
    header_file=${out_dir}/${tumour_id}.header.txt

    paste ${seg_header_file} ${annotation_header_file} > ${header_file}

    ## Remove the header and convert to BED format
    # echo "Converting segment file to BED format: ${segment}"
    bed_file=${out_dir}/${tumour_id}.bed

    tail -n +2 ${seg_file} > ${bed_file}

    ## Run bedtools intersect using the converted BED file
    output_file=${out_dir}/${tumour_id}.intersect.tsv

    bedtools intersect \
        -wa \
        -wb \
        -a ${bed_file} \
        -b ${annotation} \
        > ${output_file} 
        # -nonamecheck \
    
    final_file=${out_dir}/${tumour_id}.annotated.tsv

    ## Reattach header
    cat ${header_file} ${output_file} > ${final_file}

    ## Clean up intermediate files
    rm ${seg_header_file} ${annotation_header_file} ${header_file} ${bed_file} ${output_file} 

done

rm ${samples}

echo $(date +"%F") $(date +"%T") "Facets annotation completed."