// Here are modules related to the Aggregate Peak Analysis 
process do_juicer_apa {
    tag "${meta.id}"
    input:
        tuple val(meta), path(hic_file), path("unique.loops")
        tuple val(juicer_name), path(juicer_path)
    output:
        tuple val(meta), path("results")
    shell:
    '''
    # First convert unique.loops to bedpe format
    convert_uniqueloops_to_bedpe.awk unique.loops > unique.loops.bedpe
    
    # View the first few lines for troubleshooting...
    cat -T unique.loops.bedpe | head

    # Now do APA
    java -jar !{juicer_path} apa -k !{params.normalization} -r !{params.hic_resolution} -u !{hic_file} unique.loops.bedpe results 
    '''
}

// Write the parameters
process write_params { 
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
