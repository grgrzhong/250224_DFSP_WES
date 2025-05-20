
params.input_vcf = "/home/zhonggr/projects/250224_DFSP_WES/data/SARC/variant_calling/annovar/SARC-004-T/SARC-004-T.hg38_multianno.vcf"

process Convert_MAFtoMAF {
    
    publishDir "${params.outdir}/maf2maf", mode: 'copy', overwrite: true

    input:
    path maf_file

    output:
    path "maf2maf/${maf_file.getName()}"

    script:
    """
    maf2maf -i ${maf_file} -o maf2maf/${maf_file.getName()} --ref_genome ${params.ref_genome} --output_format maf
    """
}

