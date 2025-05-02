process GATK4_MUTECT2 {

    tag "$meta.id"

    label 'process_medium'

    input:
    tuple val(meta), path(tumour), path(tumour_index)
    tuple val(meta2), path(normal), path(normal_index)
    path(fasta)
    path(fai)
    path(dict)
    path(germline_resource)
    path(germline_resource_tbi)
    path(panel_of_normals)
    path(panel_of_normals_tbi)
    path(intervals)

    output:
    tuple val(meta), path("*.vcf.gz")     , emit: vcf
    tuple val(meta), path("*.tbi")        , emit: tbi
    tuple val(meta), path("*.bam")        , emit: bam
    tuple val(meta), path("*.stats")      , emit: stats
    tuple val(meta), path("*.f1r2.tar.gz"), optional:true, emit: f1r2
    tuple val(meta), path("*.log")        , emit: log  
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.tumour_id}"
    def interval_command = intervals ? "--intervals $intervals" : ""
    def pon_command = panel_of_normals ? "--panel-of-normals $panel_of_normals" : ""
    def gr_command = germline_resource ? "--germline-resource $germline_resource" : ""

    // Handle normal sample if provided
    def normal_command = ""
    if (normal && normal.toString() != "null" && normal.name != "null") {
        log.info "Running Mutect2 in paired tumor-normal mode for ${meta.id} with normal ${meta2.id}"
        normal_command = "--input ${normal} --normal ${meta2.normal_id}"
    } else {
        log.info "Running Mutect2 in tumor-only mode for ${meta.id} (no matching normal sample)"
    }

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK Mutect2] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        Mutect2 \\
        --input $tumour \\
        ${normal_command} \\
        --reference $fasta \\
        $pon_command \\
        $gr_command \\
        $interval_command \\
        --output ${prefix}.unfiltered.vcf.gz \\
        -bamout ${prefix}.realigned.bam \\
        --f1r2-tar-gz ${prefix}.f1r2.tar.gz \\
        --native-pair-hmm-threads ${task.cpus} \\
        --tmp-dir . \\
        $args \\
        2> ${prefix}.mutect2.log

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
    touch ${prefix}.vcf.gz.stats
    touch ${prefix}.f1r2.tar.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
