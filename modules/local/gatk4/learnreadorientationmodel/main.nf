process GATK4_LEARNREADORIENTATIONMODEL {

    tag "$meta.id"

    label 'process_low'

    input:
    tuple val(meta), path(f1r2)

    output:
    tuple val(meta), path("*.readorientationmodel.tar.gz"), emit: orientation
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK LearnReadOrientationModel] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        LearnReadOrientationModel \\
        -I ${f1r2} \\
        --output ${prefix}.readorientationmodel.tar.gz \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
