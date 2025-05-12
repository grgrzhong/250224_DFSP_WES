workflow PREPARE_SAMPLE {
    take:
        input_csv  // Path to the input CSV file

    main:
        // Create unified channel from the CSV with support for both input types
    samples = Channel.fromPath(input_csv)
        .ifEmpty { exit 1, "Samplesheet not found: ${input_csv}" }
        .splitCsv(header: true)
        .map { row ->
            // Extract sample info
            def patient_id = row.patient
            def sample_id = row.sample
            def status = row.status
            
            // Validate required fields
            if (!patient_id) error "Missing or empty 'patient' field in row: ${row}"
            if (!sample_id) error "Missing or empty 'sample' field in row: ${row}"
            if (!status) error "Missing or empty 'status' field in row: ${row}"
            
            // Process BAM and BAI paths
            def bam = file(row.bam)
            def bai = file(row.bai)
            
            // Check if files exist
            // if (!bam.exists()) error "BAM file not found: ${bam}"
            // if (!bai.exists()) error "BAI file not found: ${bai}"
            
            // Return a tuple with patient_id, sample_id, status, bam, and bai
            [
                patient_id: patient_id, 
                sample_id: sample_id, 
                status: status.toInteger(), 
                bam: bam, 
                bai:bai
            ]
        }
        // .view()
    
    samples.count().view { count -> "Total samples: ${count}" }

    // Split samples into tumour and normal
    tumour_samples = samples.filter { it.status == 1 }
    normal_samples = samples.filter { it.status == 0 }
    
    tumour_samples.count().view { count -> "Total tumour samples: ${count}" }
    normal_samples.count().view { count -> "Total normal samples: ${count}" }

    // tumour_samples
    //     .filter { it.patient_id == "DFSP-332" }
    //     .view { sample -> "Tumour sample: ${sample}" }

    // Match tumour samples with normal samples by patient ID
    tumour_normal_pairs = tumour_samples
        .map { tumour -> [tumour.patient_id, tumour] }
        .combine(
            normal_samples.map { normal -> [normal.patient_id, normal] },
            by: 0 // Match by patient_id
        )
        .map { 
            patient_id, tumour, normal ->
            def tumour_id = tumour.sample_id
            def normal_id = normal.sample_id
            def tumour_bam = tumour.bam
            def tumour_bai = tumour.bai
            def normal_bam = normal.bam
            def normal_bai = normal.bai
            
            // Create meta map required by the process
            def meta = [
                id: "${tumour_id}_vs_${normal_id}",
                patient_id: patient_id,
                tumour_id: tumour_id,
                normal_id: normal_id
            ]
            
            // Return in format expected by CNV_FACETS process
            [meta, tumour_bam, tumour_bai, normal_bam, normal_bai]
        }
    
    tumour_normal_pairs
        .count()
        .view{count -> "Total tumour-normal pairs: ${count}"}
        
    emit:
    input_samples       = samples
    tumour_normal_pairs = tumour_normal_pairs
}