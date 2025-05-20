process ANNOVAR {
    
    tag "$meta.id"
    
    label 'process_low'
    
    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path annovar_db
    val annovar_buildver
    val annovar_protocol
    val annovar_operation
    path annovar_xreffile
    
    output:
    tuple val(meta), path("*multianno.txt") , emit: multianno
    tuple val(meta), path("*.annovar.txt") , emit: annovar
    tuple val(meta), path("*.log")  , emit: log
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    perl ${annovar_db}/table_annovar.pl $vcf \\
        $annovar_db/humandb/ \\
        -buildver $annovar_buildver \\
        -out $prefix \\
        -remove \\
        -protocol $annovar_protocol \\
        -operation $annovar_operation \\
        -nastring . \\
        -polish \\
        -xreffile $annovar_xreffile \\
        --otherinfo \\
        --vcfinput \\
        $args \\
        &> ${prefix}.annovar.log
        # Create a more readable output with added fields for AD, AF, DP
    
    cat ${prefix}.${annovar_buildver}_multianno.txt | \\
        awk 'BEGIN {FS=OFS="\\t"} \\
        NR==1 {print \$0, "AD", "AF", "DP"}; \\
        NR >1 {split(\$NF, a, ":"); \$(NF+1)=a[2]; \$(NF+1)=a[3]; \$(NF+1)=a[4]; print}' \\
        > ${prefix}.annovar.txt

    """
}