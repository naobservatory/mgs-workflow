/***********
| WORKFLOW |
***********/

workflow LOAD_SAMPLESHEET {
    take:
        sample_sheet
    main:
        // Start time
        start_time = new Date()
        start_time_str = start_time.format("YYYY-MM-dd HH:mm:ss z (Z)")

        // Check pairing and validate headers
        def headers = file(sample_sheet).readLines().first().tokenize(',')*.trim()
        def expected_headers_se = ['sample', 'fastq']
        def expected_headers_pe = ['sample', 'fastq_1', 'fastq_2']
        if (headers.take(2) == expected_headers_se) {
            single_end = true
        } else if (headers.take(3) == expected_headers_pe) {
            single_end = false
        } else {
            throw new Exception("""Invalid samplesheet header. 
                Expected ${expected_headers_se.join(', ')} or ${expected_headers_pe.join(', ')}
                Found: ${headers.join(', ')}
                Please ensure the samplesheet has the correct columns in the specified order.""".stripIndent())
        }
        // Construct samplesheet channel
        if (single_end) {
            samplesheet = Channel
                .fromPath(sample_sheet)
                .splitCsv(header: true)
                .map { row -> tuple(row.sample, file(row.fastq)) }
            samplesheet_ch = samplesheet.map { sample, read -> tuple(sample, [read]) }
        } else {
            samplesheet = Channel
                .fromPath(sample_sheet)
                .splitCsv(header: true)
                .map { row -> tuple(row.sample, file(row.fastq_1), file(row.fastq_2)) }
            samplesheet_ch = samplesheet.map { sample, read1, read2 -> tuple(sample, [read1, read2]) }
        }

    emit:
        single_end = single_end
        samplesheet = samplesheet_ch
        start_time_str = start_time_str
        test_input = sample_sheet
}
