process GATK4_FILTERMUTECTCALLS {
    tag "$meta.id"
    
    label 'process_low'

    input:
    tuple val(meta), path(vcf), path(vcf_tbi), path(stats), path(orientationmodel), path(contamination), path(segmentation)
    path(fasta)
    path(fai)
    path(dict)

    output:
    tuple val(meta), path("*.vcf.gz")      , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi")  , emit: tbi
    tuple val(meta), path("*.stats")       , emit: stats
    tuple val(meta), path("*.log")         , emit: log 
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.tumour_id}"

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK FilterMutectCalls] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        FilterMutectCalls \\
        --variant $vcf \\
        --reference $fasta \\
        --ob-priors $orientationmodel \\
        --contamination-table $contamination \\
        --tumor-segmentation $segmentation \\
        --min-allele-fraction 0.01 \\
        --unique-alt-read-count 1 \\
        --stats ${prefix}.unfiltered.vcf.gz.stats \\
        --output ${prefix}.filtered.vcf.gz \\
        --tmp-dir . \\
        $args \\
        2> ${prefix}.filtermutectcalls.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
