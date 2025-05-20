process CNVKIT_REFERENCE {
    
    tag "$fasta"
    
    label 'process_low'

    input:
    path(normal_bams)
    path(normal_bais)
    path(fasta)
    path(fai)
    path(dict)
    path(targets)

    output:
    path "*.cnn"       , emit: cnn
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def normal_args = normal_bams.join(' ')
    def prefix = "cnvkit"
    """
    # Create reference using all normal samples

    cnvkit.py \\
        batch \\
        --normal ${normal_args} \\
        --fasta $fasta \\
        --targets $targets \\
        --output-reference ${prefix}.reference.cnn \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvkit: \$(cnvkit.py version | sed -e "s/cnvkit v//g")
    END_VERSIONS
    """
}
