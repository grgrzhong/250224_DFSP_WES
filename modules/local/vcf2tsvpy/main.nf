process VCF2TSVPY {
    
    tag "$meta.id"
    
    label 'process_low'
        
    input:
    tuple val(meta), path(vcf)
    
    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    vcf2tsvpy \\
        --input_vcf $vcf \\
        --out_tsv ${prefix}.tsv \\
        $args
    
    ## Remove the first lines starting with # from the TSV output
    
    grep -v "^#" ${prefix}.tsv > ${prefix}.tsv.tmp
    
    mv ${prefix}.tsv.tmp ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        vcf2tsvpy: \$(vcf2tsvpy --version 2>&1 | sed 's/vcf2tsvpy //g' || echo "unknown")
    END_VERSIONS
    """
}