process FORMAT_VCF_TUMOUR_NORMAL {

    tag "${meta.id}"

    label 'process_low'

    input:
    tuple val(meta), path(vcf), path(vcf_index)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    pcgr_reformat_vcf_tumour_normal.py \\
        --input $vcf \\
        --output ${prefix}.reformat.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        pysam: \$(python3 -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    
    """
}