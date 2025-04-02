process ANNOVAR_SOMATIC {
    tag "${meta.id}"
    
    publishDir "${params.outdir}/annovar/${meta.patient}", mode: 'copy'
    
    input:
    // Input is a tuple with metadata and VCF file
    tuple val(meta), path(vcf)
    path(annovar_db)
    
    output:
    tuple val(meta), path("${meta.id}.annovar.txt"), emit: annovar_txt
    tuple val(meta), path("${meta.id}.annovar.vcf"), emit: annotated_vcf
    tuple val(meta), path("${meta.id}.annovar.*.txt"), emit: all_annovar_files
    
    script:
    def avdbdir = annovar_db ? "--dbtype ${params.annovar_dbtype} --buildver ${params.annovar_buildver} --database_dir ${annovar_db}" : ""
    def protocols = params.annovar_protocols.join(',')
    def operations = params.annovar_operation.join(',')
    
    """
    # Convert VCF to ANNOVAR format
    convert2annovar.pl \\
        -format vcf4 \\
        ${vcf} \\
        -outfile ${meta.id}.avinput \\
        -includeinfo
    
    # Run ANNOVAR annotation
    table_annovar.pl \\
        ${meta.id}.avinput \\
        ${annovar_db} \\
        -buildver ${params.annovar_buildver} \\
        -out ${meta.id}.annovar \\
        -remove \\
        -protocol ${protocols} \\
        -operation ${operations} \\
        -nastring . \\
        -vcfinput \\
        -polish
        
    # If this is a somatic mutation file, add cancer-specific annotations if available
    if [ "${meta.status}" == "somatic" ] || [ "${meta.status}" == "tumor-normal" ]; then
        # Add cancer-specific annotations if CancerVar DB is available
        if [ -d "${annovar_db}/CancerVar" ]; then
            CancerVar.pl \\
                -f ${meta.id}.annovar.txt \\
                -o ${meta.id}.cancervar \\
                --db ${annovar_db}/CancerVar \\
                --buildver ${params.annovar_buildver}
            
            # Merge annotations
            paste ${meta.id}.annovar.txt ${meta.id}.cancervar.txt > ${meta.id}.combined_annotations.txt
        fi
    fi
    """
}

// Example of a workflow that uses this process
workflow {
    // Define parameters with defaults
    params.vcf = null
    params.outdir = "results/annotation"
    params.annovar_db = "/path/to/annovar/humandb"
    params.annovar_buildver = "hg38"
    params.annovar_dbtype = "refGene"
    params.annovar_protocols = ["refGene", "knownGene", "ensGene", "gnomad211_exome", 
                               "clinvar_20220320", "cosmic70", "dbnsfp42c"]
    params.annovar_operation = ["g", "g", "g", "f", "f", "f", "f"]
    
    // Input channel
    if (params.vcf) {
        vcf_ch = Channel.fromPath(params.vcf)
            .map { vcf -> 
                def meta = [
                    id: vcf.simpleName,
                    patient: vcf.simpleName.split('_')[0],
                    status: params.vcf_type ?: "somatic" // default to somatic
                ]
                return [meta, vcf]
            }
    }
    
    // Database
    annovar_db_ch = Channel.value(file(params.annovar_db))
    
    // Run ANNOVAR
    ANNOVAR_SOMATIC(vcf_ch, annovar_db_ch)
}