process CNV_FACETS_PAIRED {
    
    tag "${meta.id}"

    label "process_medium"
    
    input:
    tuple val(meta), path(input_tumour), path(input_tumour_index), path(input_normal), path(input_normal_index)
    path dbsnp
    path dbsnp_tbi
    
    output:
    tuple val(meta), path("*.cnv.png"),     emit: cnv
    tuple val(meta), path("*.cov.pdf"),     emit: cov
    tuple val(meta), path("*.csv.gz"),      emit: csv
    tuple val(meta), path("*.spider.pdf"),  emit: spider
    tuple val(meta), path("*.vcf.gz"),      emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"),  emit: vcf_tbi
    
    script:
    def prefix = meta.tumour_id

    """
    # Run CNV-FACETS for tumor-normal pair
    cnv_facets.R \\
        --snp-normal ${input_normal} \\
        --snp-tumour ${input_tumour} \\
        --snp-vcf ${dbsnp} \\
        --snp-nprocs ${task.cpus} \\
        --out ${prefix}
    """
}
