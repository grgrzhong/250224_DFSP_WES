#!/bin/bash

## Convert VCF to TSV using vcf2tsvpy
conda install -c bioconda vcf2tsvpy

conda activate vcf2tsv

## Convert a VCF file to TSV format
input_vcf=/mnt/m/WES/SARC/CNV-call/facets/SARC-022-T/SARC-022-T.vcf.gz
output_tsv=/mnt/m/WES/SARC/CNV-call/facets/SARC-022-T/SARC-022-T.tsv

zgrep "^#CHROM" ${input_vcf}

## View INFO column of VCF
# echo "Showing INFO column sample from VCF file:"
# zcat ${input_vcf} | grep -v "^#" | head -5 | cut -f 8

# ## Alternative using bcftools (if installed)
# if command -v bcftools &> /dev/null; then
#     echo "INFO fields and descriptions:"
#     bcftools view -h ${input_vcf} | grep "^##INFO"
    
#     echo "Sample INFO column values:"
#     bcftools query -f '%INFO\n' ${input_vcf} | head -5
# fi

## Check if input VCF file exists
if [[ ! -f ${input_vcf} ]]; then
    echo "Input VCF file does not exist: ${input_vcf}"
    exit 1
fi

vcf2tsvpy \
    --input_vcf ${input_vcf} \
    --out_tsv ${output_tsv}

## Remove the first lines starting with # from the TSV output
grep -v "^#" ${output_tsv} > ${output_tsv}.tmp && mv ${output_tsv}.tmp ${output_tsv}

echo "Conversion complete and comment lines removed: ${output_tsv}"
