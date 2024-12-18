/***********
| WORKFLOW |
***********/

workflow LOAD_SAMPLESHET {
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
