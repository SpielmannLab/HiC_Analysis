// Here we define modules related to getting TAD boundaries
process write_params{ 
    // write input parameters to file
    label 'tiny'
    output:
        path "*.txt"
    shell:
    '''
    suffix=$(date "+%Y_%m_%d_%H_%M")
    echo "Parameters of the job" >> parameters_sc_p_sample_${suffix}.txt
    echo "!{params}" | tr , '\n' >> parameters_sc_p_sample_${suffix}.txt
    '''
}

process perform_arrowhead {
    // perform arrowhead tad calling
    tag "${meta.samplename}"
    input:
        tuple val(meta), path(hic_file)
        tuple val(juicer_name), path(juicer_path)
    output:
        tuple val(meta), path("**.bedpe"), emit: arrowhead_out, optional: true
    shell:
    '''
    if [[ !{params.chromosome} = "null" ]]; then
        java -jar !{juicer_path} arrowhead \
             -k !{params.normalization} \
             -r !{params.hic_resolution} \
             -k !{hic_file} \
             --threads 0 \
             ./
    else
        java -jar !{juicer_path} arrowhead \
             -k !{params.normalization} \
             -r !{params.hic_resolution} \
             -c !{params.chromosome} \
             -k !{hic_file} \
             --threads 0 \
             ./
    fi
    '''
}
