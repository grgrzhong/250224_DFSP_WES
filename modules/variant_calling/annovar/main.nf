process ANNOVAR {
    tag "$meta.id"
    label 'process_medium'
    
    // Annovar doesn't have a standard conda or container, usually installed locally
    // container "quay.io/biocontainers/annovar:xxxxxxxx" // Use if available
    
    input:
    tuple val(meta), path(vcf), path(tbi)
    path annovar_db
    val build_version
    val protocols
    val operations
    
    output:
    tuple val(meta), path("*_annovar.txt"), emit: annovar_txt
    tuple val(meta), path("*_multianno.txt"), emit: multianno
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    perl /path/to/annovar/table_annovar.pl $vcf \\
        $annovar_db \\
        -buildver $build_version \\
        -out $prefix \\
        -remove \\
        -protocol $protocols \\
        -operation $operations \\
        -nastring . \\
        -polish \\
        -xreffile /path/to/annovar/example/gene_fullxref.txt \\
        --otherinfo \\
        --vcfinput \\
        $args
    
    # Create a more readable output with added fields for AD, AF, DP
    less -S ${prefix}.${build_version}_multianno.txt | \\
        awk 'BEGIN {FS=OFS="\\t"} \\
        NR==1 {print \$0, "AD", "AF", "DP"}; \\
        NR >1 {split(\$NF, a, ":"); \$(NF+1)=a[2]; \$(NF+1)=a[3]; \$(NF+1)=a[4]; print}' \\
        > ${prefix}_annovar.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annovar: \$(grep "Version" \$(dirname \$(which table_annovar.pl))/README.txt | sed 's/.*Version: //; s/ .*//')
    END_VERSIONS
    """
}