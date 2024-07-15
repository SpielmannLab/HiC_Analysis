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
        tuple val(meta), path(hic_file), path("arrowhead_tads/*.bedpe"), emit: arrowhead_out
        path("output.log"), optional: true
    shell:
    '''
    if [[ !{params.chromosome} = "null" ]]; then
        java -jar !{juicer_path} arrowhead -m !{params.sliding_window_size} -k !{params.normalization} -r !{params.hic_resolution} !{hic_file} arrowhead_tads --threads !{task.cpus} --ignore_sparsity | tee output.log
    else
        java -jar !{juicer_path} arrowhead -m !{params.sliding_window_size} -k !{params.normalization} -r !{params.hic_resolution} -c !{params.chromosome} !{hic_file} arrowhead_tads --threads !{task.cpus} --ignore_sparsity | tee output.log
    fi
    # see the output
    tree 
    '''
}

process perform_diffDomain {
    tag "${meta1.samplename}_vs_${meta2.samplename}"
    input:
        tuple val(meta1), path("hic_file1.hic"), path("tad_file1.bedpe")
        tuple val(meta2), path("hic_file2.hic"), path("tad_file2.bedpe")
        tuple val(diffDomain_name), path(diffDomain_path)
    output:
        tuple val(meta1), val(meta2), path("hic_file1.hic"), path("hic_file2.hic"), path("dvsd.fdr_bh.tsv"), path("dvsd.fdr_bh.tsv_types.txt")
    shell:
    '''
    # prepare TAD.bedpe files for diffDomain
    # Remove # from the 1st header
    # Remove all other headers
    sed -i '1 s/^#//1' tad_file1.bedpe
    grep -v "#" tad_file1.bedpe > tad_file1_for_diffDomain.bedpe
    # Same for the second one
    sed -i '1 s/^#//1' tad_file2.bedpe
    grep -v "#" tad_file2.bedpe > tad_file2_for_diffDomain.bedpe
    
    # run diffdomains python script
    python !{diffDomain_path}/diffdomains.py dvsd multiple hic_file1.hic hic_file2.hic tad_file1_for_diffDomain.bedpe --reso !{params.hic_resolution} --hicnorm !{params.normalization} --ncore !{task.cpus} --ofile "dvsd.tsv"

    # adjust for fdr using Benjamini Hochburg
    python !{diffDomain_path}/diffdomains.py adjustment fdr_bh dvsd.tsv dvsd.fdr_bh.tsv --filter false --ncore !{task.cpus}

    # Classify the differential tads
    python !{diffDomain_path}/classification.py -t tad_file2_for_diffDomain.bedpe -d dvsd.fdr_bh.tsv -o "./"
    '''
}

process visualize_diffDomain {
    tag "${meta1.samplename}_vs_${meta2.samplename}"
    input:
        tuple val(meta1), val(meta2), path("hic_file1.hic"), path("hic_file2.hic"), path("dvsd.fdr_bh.tsv"), path("dvsd.fdr_bh.tsv_types.txt")
        tuple val(diffDomain_name), path(diffDomain_path)
    output:
        path "to_plot.txt"
        path "*.pdf", optional: true
    shell:
    '''
    # Plot all the statistically relevant ones
    awk 'NR==1 {} $9 == "1" {print "chr"$1,$2,$3}' dvsd.fdr_bh.tsv_types.txt > to_plot.txt
    cat to_plot.txt
    if [[ -s to_plot.txt ]]; then
        while IFS= read -r line; do
            python !{diffDomain_path}/diffdomains.py visualization ${line} hic_file1.hic hic_file2.hic --ncore !{task.cpus} --hicnorm !{params.normalization} --reso !{params.hic_resolution} --ofile ${line// /_}
        done < to_plot.txt
    fi
    '''
}

process perform_insulation_score {
    // perform arrowhead tad calling
    tag "${meta.samplename}"
    input:
        tuple val(meta), path(hic_file)
        tuple val(juicer_name), path(juicer_path)
    output:
        tuple val(meta), path("**.bedpe"), emit: arrowhead_out, optional: true
        path("output.log"), optional: true
    shell:
    '''

    '''
}
