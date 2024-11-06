process BATCH_POD_5 {
    label "pandas"
    label "batch_pod_5"

    input:
        path(pod_5_dir)
        val batch_size
    output:
        path('pod_5_batch/*')
    //  ## batch_pod5.py --batch_size !{batch_size} --pod5_dir !{pod_5_dir} --output_dir pod_5_batch/
    shell:
        '''
        batch_num=0
        current_batch_size=0
        current_batch=()

        for fname in $(ls !{pod_5_dir}); do
            path="!{pod_5_dir}/$fname"
            size=$(stat -c %s "$path")

            if [ $((current_batch_size + size)) -gt !{batch_size} ]; then
                batch_num=$((batch_num + 1))
                batch_dir="pod_5_batch/batch_$(printf "%03d" $batch_num)"
                mkdir -p "$batch_dir"
                for pod5_file in "${current_batch[@]}"; do
                    ln -s "$pod5_file" "$batch_dir/$(basename $pod5_file)"
                done
                current_batch=()
                current_batch_size=0
            fi

            current_batch+=("$path")
            current_batch_size=$((current_batch_size + size))
        done

        # Handle remaining files in last batch
        if [ ${#current_batch[@]} -gt 0 ]; then
            batch_num=$((batch_num + 1))
            batch_dir="pod_5_batch/batch_$(printf "%03d" $batch_num)"
            mkdir -p "$batch_dir"
            for pod5_file in "${current_batch[@]}"; do
                ln -s "$pod5_file" "$batch_dir/$(basename $pod5_file)"
            done
        fi
        '''
}