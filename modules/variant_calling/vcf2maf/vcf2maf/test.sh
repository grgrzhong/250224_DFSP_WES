#/bin

# Open a shell for debugging
singularity shell \
    -B /home/zhonggr/projects/250224_DFSP_WES/data/SARC/annotation/annovar/SARC-004-T:/data \
    /home/zhonggr/projects/250224_DFSP_WES/containers/singularity/vcf2maf.sif

# Verify the VEP installation
singularity exec /home/zhonggr/projects/250224_DFSP_WES/containers/singularity/vcf2maf.sif \
    /opt/conda/bin/vep --help

singularity exec \
    -B /home/zhonggr/projects/250224_DFSP_WES/data/SARC/annotation/annovar/SARC-004-T:/data \
    /home/zhonggr/projects/250224_DFSP_WES/containers/singularity/vcf2maf.sif \
    perl /opt/conda/bin/vcf2maf.pl \
    --input-vcf /data/SARC-004-T.hg38.multianno.vcf \
    --output-maf /data/SARC-004-T.hg38.multianno.maf \
    --tumor-id SARC-004-T \
    --normal-id SARC-004-N \
    --vep-path /opt/conda/bin/vep \
    --vep-data /data/vep_data \
    --ncbi-build GRCh38