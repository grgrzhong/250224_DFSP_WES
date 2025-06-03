#!/bin/bash

singularity shell \
    --bind /home/zhonggr/projects/250224_DFSP_WES/data/mfs:/data/mfs \
    --bind /home/zhonggr/projects/250224_DFSP_WES/data/reference:/reference \
    /home/zhonggr/projects/250224_DFSP_WES/containers/singularity/cnv_facets.sif

export dbsnp="/reference/dbSNP.vcf.gz"
export input_normal="/data/mfs/bam/SARC-022-N/SARC-022-N_recalibrated.bam"
export input_tumour="/data/mfs/bam/SARC-022-T/SARC-022-T_recalibrated.bam"
export prefix="/data/mfs/cnv/facets/SARC-022-T/SARC-022-T"

cnv_facets.R \
    --snp-normal ${input_normal} \
    --snp-tumour ${input_tumour} \
    --snp-vcf ${dbsnp} \
    --snp-nprocs 4 \
    --out ${prefix}


vcftools --gzvcf /home/zhonggr/projects/250224_DFSP_WES/data/mfs/cnv/facets/SARC-002-T/SARC-022-T.vcf.gz --out /home/zhonggr/projects/250224_DFSP_WES/data/mfs/cnv/facets/SARC-002-T/SARC-022-T --extract-FORMAT-info GT

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\n' \
/home/zhonggr/projects/250224_DFSP_WES/data/mfs/cnv/facets/SARC-002-T/SARC-022-T.vcf.gz \
> /home/zhonggr/projects/250224_DFSP_WES/data/mfs/cnv/facets/SARC-002-T/SARC-022-T.txt

zcat /home/zhonggr/projects/250224_DFSP_WES/data/mfs/cnv/facets/SARC-002-T/SARC-022-T.vcf.gz | grep "^#CHROM" | \
sed 's/^#//' > /home/zhonggr/projects/250224_DFSP_WES/data/mfs/cnv/facets/SARC-002-T/SARC-022-T.txt

zcat /home/zhonggr/projects/250224_DFSP_WES/data/mfs/cnv/facets/SARC-002-T/SARC-022-T.vcf.gz | \
grep -v "^##" > /home/zhonggr/projects/250224_DFSP_WES/data/mfs/cnv/facets/SARC-002-T/SARC-022-T.txt