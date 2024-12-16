/***********
| WORKFLOW |
***********/

workflow LOAD_SAMPLESHEET {
    // Start time
    start_time = new Date()
    start_time_str = start_time.format("YYYY-MM-dd HH:mm:ss z (Z)")

    // Check if grouping column exists in samplesheet
    check_grouping = new File(params.sample_sheet).text.readLines()[0].contains('group') ? true : false
    if (params.grouping != check_grouping) {
        if (params.grouping && !check_grouping) {
            throw new Exception("Grouping enabled in config file, but group column absent from samplesheet.")
        } else if (!params.grouping && check_grouping) {
            throw new Exception("Grouping is not enabled in config file, but group column is present in the samplesheet.")
        }
    }
    take:
        sample_sheet
    main:
        if (params.single_end) {
            if (params.grouping) {
                samplesheet = Channel
                    .fromPath(params.sample_sheet)
                    .splitCsv(header: true)
                    .map { row -> tuple(row.sample, file(row.fastq), row.group) }
                samplesheet_ch = samplesheet.map { sample, read, group -> tuple(sample, [read]) }
                group_ch = samplesheet.map { sample, read, group -> tuple(sample, group) }
            } else {
                samplesheet = Channel
                    .fromPath(params.sample_sheet)
                    .splitCsv(header: true)
                    .map { row -> tuple(row.sample, file(row.fastq)) }
                samplesheet_ch = samplesheet.map { sample, read -> tuple(sample, [read]) }
                group_ch = Channel.empty()
            }
        } else {
            if (params.grouping) {
                samplesheet = Channel
                    .fromPath(params.sample_sheet)
                    .splitCsv(header: true)
                    .map { row -> tuple(row.sample, file(row.fastq_1), file(row.fastq_2), row.group) }
                samplesheet_ch = samplesheet.map { sample, read1, read2, group -> tuple(sample, [read1, read2]) }
                group_ch = samplesheet.map { sample, read1, read2, group -> tuple(sample, group) }
            } else {
                samplesheet = Channel
                    .fromPath(params.sample_sheet)
                    .splitCsv(header: true)
                    .map { row -> tuple(row.sample, file(row.fastq_1), file(row.fastq_2)) }
                samplesheet_ch = samplesheet.map { sample, read1, read2 -> tuple(sample, [read1, read2]) }
                group_ch = Channel.empty()
                }
            }
    emit:
        samplesheet = samplesheet_ch
        group = group_ch
}
