#!/usr/bin/env bash

# Arguments
run_id=${1}
sample_process=${2}
sample_pattern=${3}
job_process=${4}
s3_prefix=${5}

# Break if any variable missing
if [[ -z ${run_id} ]]; then
    echo "Error: no run ID given."
    exit 1
fi
if [[ -z ${sample_process} ]]; then
    echo "Error: no sample process given."
    exit 1
fi
if [[ -z ${sample_pattern} ]]; then
    echo "Error: no sample pattern given."
    exit 1
fi
if [[ -z ${job_process} ]]; then
    echo "Error: no job process given."
    exit 1
fi

# Define auxiliary commands
run="nextflow log ${run_id} -fields hash,name,exit,status,script"

# Get sample IDs
sample_ids=$(${run} -filter "name=~ /.*${sample_process}.*/ && status == 'FAILED'" | grep -o "${sample_pattern}" | sort | uniq)
echo "Sample IDs: $(echo ${sample_ids})"

# Get job IDs
job_ids=$(for id in ${sample_ids}; do ${run} -filter "name=~ /.*${job_process}.*/ && script.contains('${id}')" | grep ${job_process} | cut -f 1; done)
echo "Job IDs: $(echo ${job_ids})"

# Get directories
job_dirs=$(for id in ${job_ids}; do dir1=$(dirname ${id}); dir2=$(aws s3 ls ${s3_prefix}/${id} | awk '{print $2}'); if [[ -z ${dir2} ]]; then continue; fi; job_dir="${s3_prefix}/${dir1}/${dir2}"; echo ${job_dir}; done)
echo "Job directories to delete: $(echo ${job_dirs})"
echo

# Delete directories
for dir in ${job_dirs}; do
    echo ${dir}
    aws s3 rm --recursive ${dir}
    echo
done
