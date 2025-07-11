process COPY_REF_DIRECTORY {
    label "single"
    label "coreutils"
    maxForks 1
    input:
        val(dir)
        val(outname)
    output:
        path("${outname}")
    script:
        '''
        if [ ! -d "!{outname}" ]; then
            if [[ "!{dir}" == s3://* ]]; then
                aws s3 cp --recursive !{dir} /scratch/tmp_krakendb
            else
                cp -r !{dir} /scratch/tmp_krakendb
            fi
            mv /scratch/tmp_krakendb !{outname}
        fi
        '''
}
