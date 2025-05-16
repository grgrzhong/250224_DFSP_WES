process ANNOVAR_EXPORT {
    
    tag "$meta.id"
    
    label 'process_low'
    
    input:
    tuple val(meta), path(maf)
    val annovar_buildver
    
    output:
    tuple val(meta), path("*annovar.txt"), emit: txt

    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    # Create a more readable output with added fields for AD, AF, DP
    cat ${prefix}.${annovar_buildver}_multianno.txt | \\
        awk 'BEGIN {FS=OFS="\\t"} \\
        NR==1 {print \$0, "AD", "AF", "DP"}; \\
        NR >1 {split(\$NF, a, ":"); \$(NF+1)=a[2]; \$(NF+1)=a[3]; \$(NF+1)=a[4]; print}' \\
        > ${prefix}.annovar.txt

    """
}