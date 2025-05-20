process VCF2MAF {
    
    tag "${meta.id}"

    label "process_medium"
    
    input:
    tuple val(meta), path(vcf)
    path fasta      // Required
    path vep_cache  // Required for VEP running. A default of /.vep is supplied.

    output:
    tuple val(meta), path("*.maf"), emit:maf
    path  "version.yml",            emit:versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    // def assembl_version = params.genome 
    def prefix = task.ext.prefix ?: "${meta.id}"
    def vep_cache_cmd = vep_cache ? "--vep-data $vep_cache" : ""     // If VEP is present, it will find it and add it to commands otherwise blank
    def VERSION = '1.6.22' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    #!/bin/bash

    if [ "$vep_cache" ]; then
        VEP_CMD="--vep-path \$(dirname \$(type -p vep))"
        VEP_VERSION=\$(echo -e "\\n    ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')")
    else
        VEP_CMD=""
        VEP_VERSION=""
    fi

    vcf2maf.pl \\
        $args \\
        \$VEP_CMD \\
        $vep_cache_cmd \\
        --offline \\
        --ref-fasta $fasta \\
        --tumor-id ${meta.tumour_id} \\
        --normal-id ${meta.normal_id} \\
        --vep-data ${vep_cache} \\
        --input-vcf ${vcf} \\
        --output-maf ${prefix}.maf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2maf: $VERSION\$VEP_VERSION
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.6.22' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    if [ "$vep_cache" ]; then
        VEP_VERSION=\$(echo -e "\\n    ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')")
    else
        VEP_VERSION=""
    fi

    touch ${prefix}.maf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2maf: $VERSION\$VEP_VERSION
    END_VERSIONS
    """
}