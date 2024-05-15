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
    tuple val(meta), path("*.cool"), emit: samples
    path("*.out")
	shell:
    '''
    hicConvertFormat -m file.hic --inputFormat hic --outputFormat cool -o !{meta.samplename}_!{meta.resolution}.cool --resolution !{meta.resolution} --load_raw_values |& tee !{meta.samplename}_!{meta.resolution}_hicexplorer.out
    '''
}

process cool_normalizeNcorrect {
  tag "${meta.samplename}"
  input:
    tuple val(meta), path("file.cool")
	output:
    tuple val(meta), path("*.cool"), emit: samples
    path("*.png")
    path("*.out")
	shell:
    '''
  hicNormalize --matrices file.cool --normalize !{meta.hicexplorer_normalization} -o !{meta.samplename}_!{meta.resolution}_norm.cool |& tee !{meta.samplename}_!{meta.resolution}_hicNormalize.out
  hicCorrectMatrix diagnostic_plot -m !{meta.samplename}_!{meta.resolution}_norm.cool -o !{meta.samplename}_!{meta.resolution}_diagnostics.png |& tee !{meta.samplename}_!{meta.resolution}_hicDiagnostics.out
  hicCorrectMatrix correct --m !{meta.samplename}_!{meta.resolution}_norm.cool -o !{meta.samplename}_!{meta.resolution}_corr.cool --correctionMethod !{meta.hicexplorer_correction_method} --filterThreshold !{meta.hicexplorer_threshold_low} !{meta.hicexplorer_threshold_high} |& tee !{meta.samplename}_!{meta.resolution}_hicCorrect.out
    '''
}
