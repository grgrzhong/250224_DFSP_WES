process GATK4_MUTECT2_NORMAL {

    tag "$meta.id"

    label 'process_medium'

    input:
    tuple val(meta), path(input), path(input_index)
    path(fasta)
    path(fai)
    path(dict)
    path(intervals)

    output:
    tuple val(meta), path("*.vcf.gz")     , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi") , emit: tbi
    tuple val(meta), path("*.log")        , emit: log  
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.tumour_id}"

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK Mutect2] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        Mutect2 \\
        --input $input \\
        --reference $fasta \\
        -max-mnp-distance 0 \\
        -L $intervals \\
        --output ${prefix}.pon.vcf.gz \\
        --native-pair-hmm-threads ${task.cpus} \\
        $args \\
        2> ${prefix}.mutect2normal.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
