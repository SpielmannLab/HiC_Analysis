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

process hic_to_cool_raw {
  // Convert *.hic file to cool file with raw counts
  tag "${meta.samplename}"
  input:
    tuple val(meta), path("file.hic")
	output:
    path("*.cool"), emit: cool_files
    path("*.out")
	shell:
    '''
    hicConvertFormat -m file.hic --inputFormat hic --outputFormat cool -o !{meta.samplename}_!{meta.resolution}.cool --resolution !{meta.resolution} |& tee !{meta.samplename}_!{meta.resolution}_hicexplorer.out
    '''
}

process cool_normalizeNcorrect {
  input:
    path cool_files
	output:
    path "*.cool"
    path "*.png"
    path "*.out"
	shell:
    '''
    in_cool_files=!{cool_files}
    norm_cool_files=${in_cool_files/.cool/_norm.cool}
    diagnostic_plot_files=${norm_cool_files/_norm.cool/_diagnostics.png}
    corr_cool_files=${norm_cool_files/_norm.cool/_corr.cool}

    hicNormalize --matrices ${in_cool_files} --normalize !{params.hicexplorer_normalization} -o ${norm_cool_files} |& tee hicNormalize.out
    hicCorrectMatrix diagnostic_plot -m ${norm_cool_files} -o ${diagnostic_plot_files} |& tee hicDiagnostics.out
    hicCorrectMatrix correct -m ${norm_cool_files} -o ${corr_cool_files} --correctionMethod !{params.hicexplorer_correction_method} --filterThreshold !{params.hicexplorer_threshold_low} !{params.hicexplorer_threshold_high} |& tee hicCorrect.out
    '''
}
