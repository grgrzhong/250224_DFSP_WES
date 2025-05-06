process GATK4_GENOMICSDBIMPORT {

    tag "$meta.id"

    label 'process_medium'

    input:
    tuple val(meta), path(vcfs)        // List of VCF files
    tuple val(meta2), path(tbis)        // List of tbi files
    path fasta
    path fai
    path dict
    path interval

    output:
    tuple val(meta), path("pon_db") , emit: pon_db
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--merge-input-intervals true'
    
    // Dynamically create input string for all VCF files
    def vcf_inputs = vcfs.collect{ "-V $it" }.join(' ')

    def avail_mem = 64
    if (!task.memory) {
        log.info '[GATK GenomicsDBImport] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        GenomicsDBImport \\
        $vcf_inputs \\
        -R $fasta \\
        -L $interval \\
        --genomicsdb-workspace-path pon_db \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

   stub:
    """
    mkdir -p pon_db
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
