process BATCH_POD_5 {
    label "pandas"
    label "batch_pod_5"

    input:
        path(pod_5_dir)
        val batch_size
    output:
        path('pod_5_batch/*')

    shell:
        '''
            batch_pod5.py --batch_size !{batch_size} --pod5_dir !{pod_5_dir} --output_dir pod_5_batch/
        '''
}