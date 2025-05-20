process GATK4_FUNCOTATOR_EXPORT {
    
    tag "$meta.id"

    label 'process_low'

    input:
    tuple val(meta), path(maf)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    cat ${prefix}.annotated.maf.gz | grep -v "#" > ${prefix}.annotated.tsv

    """
}
