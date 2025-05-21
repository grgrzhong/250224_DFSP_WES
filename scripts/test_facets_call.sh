#!/bin

singularity shell \
    --bind /home/zhonggr/projects/250224_DFSP_WES/data/:/data \
    --bind /home/zhonggr/projects/250224_DFSP_WES/data/reference:/reference \
    /home/zhonggr/projects/250224_DFSP_WES/containers/singularity/cnv_facets.sif

tumour_bam="/data/wes/preprocessing/recalibrated/DFSP-336-T/DFSP-336-T_recalibrated.bam"
tumour_bai="/data/wes/preprocessing/recalibrated/DFSP-336-T/DFSP-336-T_recalibrated.bai"

normal_bam="/data/wes/preprocessing/recalibrated/DFSP-336-N/DFSP-336-N_recalibrated.bam"
normal_bai="/data/wes/preprocessing/recalibrated/DFSP-336-N/DFSP-336-N_recalibrated.bai"

dbsnp="/reference/dbSNP.vcf.gz"
dbsnp_index="/reference/dbSNP.vcf.gz.tbi"

prefix="/data/test/cnv/facets/DFSP-336-T/DFSP-336-T"

cnv_facets.R \
    --snp-normal ${normal_bam} \
    --snp-tumour ${tumour_bam} \
    --snp-vcf ${dbsnp} \
    --snp-nprocs 4 \
    --out ${prefix}