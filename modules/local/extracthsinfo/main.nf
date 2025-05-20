process EXTRACT_HSINFO {
    tag "Combine HS metrics"
    label 'process_low'

    publishDir(
        path: { "${params.outdir}/reports/combined" },
        mode: params.publish_dir_mode,
        pattern: "*.txt"
    )

    input:
    path(hs_metrics_files)  // All HS metrics files as a list

    output:
    path("HS_combined.txt"), emit: combined_metrics
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Extract the header line (line 7) from the first file
    head -n 7 `ls *.txt | head -n 1` | tail -n 1 > heading.txt

    # Extract line 8 from each file with the sample ID
    for file in *.txt; do
        sample_id=\$(basename \$file | sed 's/_hs_metrics.txt//')
        awk -v sample="\$sample_id" 'NR==8 {print sample "\\t" \$0}' \$file >> metrics.txt
    done

    # Create final file with header
    cat heading.txt > HS_combined.txt
    sed -i '1s/^/Sample_ID\\t/' HS_combined.txt
    cat metrics.txt >> HS_combined.txt

    # Clean up
    rm heading.txt metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n 1 | sed 's/.*version //; s/ .*//')
    END_VERSIONS
    """
}