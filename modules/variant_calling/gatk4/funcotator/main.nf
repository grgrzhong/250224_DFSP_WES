process GATK4_FUNCOTATOR {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::gatk4=4.4.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0':
        'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf)
    path  fasta
    path  fai
    path  dict
    path  funcotator_source
    val   genome_ver
    val   output_format

    output:
    tuple val(meta), path("*.{maf.gz,vcf.gz}"), emit: annotated_variants
    tuple val(meta), path("*.{maf.gz.tbi,vcf.gz.tbi}"), optional:true, emit: index
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = output_format == "MAF" ? "maf.gz" : "vcf.gz"
    
    """
    gatk Funcotator \\
        --variant $vcf \\
        --reference $fasta \\
        --ref-version $genome_ver \\
        --data-sources-path $funcotator_source \\
        --output ${prefix}.funcotator.$extension \\
        --output-file-format $output_format \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
