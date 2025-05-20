process GATK4_FUNCOTATOR {
    
    tag "$meta.id"

    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path  fasta
    path  fai
    path  dict
    path  intervals
    path  funcotator_source
    val   funcotator_ref_version

    output:
    tuple val(meta), path("*.maf.gz"), emit: maf
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val(meta), path("*.log")   , emit: log
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    gatk Funcotator \\
        --reference $fasta \\
        --variant $vcf \\
        --output ${prefix}.annotated.maf.gz \\
        -L $intervals \\
        --output-file-format MAF \\
        --data-sources-path $funcotator_source \\
        --ref-version $funcotator_ref_version \\
        --remove-filtered-variants true \\
        $args \\
        &> ${prefix}.funcotator.log

    cat ${prefix}.annotated.maf.gz | grep -v "#" > ${prefix}.annotated.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
