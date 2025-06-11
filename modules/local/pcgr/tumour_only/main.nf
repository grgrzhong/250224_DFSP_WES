process PCGR_TUMOUR_ONLY {
    tag "$meta.id"
    label 'process_high'
    
    conda "bioconda::pcgr=1.4.1"
    
    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path(vep_dir)
    path(refdata_dir)
    
    output:
    tuple val(meta), path("${meta.id}_tumour_only/*"), emit: results
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_tumour_only"
    """
    mkdir -p ${prefix}
    
    pcgr \\
        --input_vcf $vcf \\
        --vep_dir $vep_dir \\
        --refdata_dir $refdata_dir \\
        --output_dir ${prefix} \\
        --genome_assembly grch38 \\
        --sample_id ${meta.id} \\
        --assay WES \\
        --effective_target_size_mb 34 \\
        --tumor_only \\
        --tumor_dp_tag TDP \\
        --tumor_af_tag TAF \\
        --tumor_dp_min 20 \\
        --tumor_af_min 0.05 \\
        --estimate_tmb \\
        --tmb_dp_min 20 \\
        --tmb_af_min 0.05 \\
        --estimate_msi \\
        --estimate_signatures \\
        --vcf2maf \\
        --ignore_noncoding \\
        --force_overwrite \\
        $args
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pcgr: \$(pcgr --version 2>&1 | grep -oP 'PCGR v\\K[0-9.]+')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_tumour_only"
    """
    mkdir -p ${prefix}
    touch ${prefix}/${meta.id}.pcgr_acmg.grch38.html
    touch ${prefix}/${meta.id}.pcgr_acmg.grch38.json.gz
    touch ${prefix}/${meta.id}.pcgr_acmg.grch38.maf
    touch ${prefix}/${meta.id}.pcgr_acmg.grch38.snvs_indels.tiers.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pcgr: \$(echo "1.4.1")
    END_VERSIONS
    """
}