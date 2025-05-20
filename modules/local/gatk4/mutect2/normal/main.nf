process GATK4_MUTECT2_NORMAL {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(normal), path(normal_index)
    path(fasta)
    path(fai)
    path(dict)
    path(intervals)

    output:
    tuple val(meta), path("*.vcf.gz")     , emit: vcf
    tuple val(meta), path("*.tbi")        , emit: tbi
    tuple val(meta), path("*.log")        , emit: log  
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.normal_id}"
    // Only use intervals for WES assay type
    def use_intervals = (params.assay_type == "wes") && intervals
    def interval_command = use_intervals ? "--intervals $intervals" : ""
    
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK Mutect2] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        Mutect2 \\
        --reference $fasta \\
        --input $normal \\
        -max-mnp-distance 0 \\
        $interval_command \\
        --output ${prefix}.pon.vcf.gz \\
        --native-pair-hmm-threads ${task.cpus} \\
        $args \\
        2> ${prefix}.mutect2.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}