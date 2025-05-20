#/bin

# Open a shell for debugging
singularity shell \
    -B /home/zhonggr/projects/250224_DFSP_WES/data/sarc/annotation/annovar/SARC-004-T:/data \
    -B /home/zhonggr/projects/250224_DFSP_WES/data/reference:/reference \
    /home/zhonggr/projects/250224_DFSP_WES/containers/singularity/vcf2maf.sif

fasta="/reference/Gencode/gencode.hg38.v36.primary_assembly.fa"
vep_cache="/reference/ensembl_vep"

perl vcf2maf.pl \
        --ref-fasta ${fasta} \
        --tumor-id SARC-004-T \
        --normal-id SARC-004-N \
        --vep-data ${vep_cache} \\
        --input-vcf /data/SARC-004-T.hg38.multianno.vcf \\
        --output-maf /data/SARC-004-T.maf