// write parameters to file
process write_params{ // write input parameters to file
    output:
        path "*.txt"
    shell:
    '''
    suffix=$(date "+%Y_%m_%d_%H_%M")
    echo "Parameters of the job" >> parameters_compare_hic_matrices_${suffix}.txt
    echo "!{params}" | tr , '\n' >> parameters_compare_hic_matrices_${suffix}.txt
    '''
}


// Convert to fanc file
process juicer_to_fanc {
    label "publish"
    label "fanc"
    label "multicore"
    tag "${samplename}"
    input:
        tuple val(samplename), path(hic_file)
    output:
        tuple val(samplename), path("fanc.hic") 
    script:
    """
    TMPDIR="\$PWD"
    fanc hic --threads=${task.cpus} --deepcopy "${hic_file}@${params.hic_resolution}@NONE" "fanc.hic"
    """
}

// Convert to fanc file
process fanc_downsample {
    label "publish"
    label "fanc"
    label "multicore"
    tag "${samplename}/${refname}"
    input:
        tuple val(samplename), path("query.hic"), val(refname), path("ref.hic")
    output:
        tuple val(samplename), path("fanc.hic")
    script:
    """
    TMPDIR="\$PWD"
    fanc hic --threads=${task.cpus} --downsample=ref.hic query.hic "fanc.hic"
    """
}

// Convert to fanc file
process merge_fanc {
    label "publish"
    label "fanc"
    label "multicore"
    tag "${meta.merge_name}"
    input:
        tuple val(meta), path("fanc?.hic")
    output:
        tuple val(meta), path("fanc_merged.hic") 
    script:
    """
    TMPDIR="\$PWD"
    fanc hic --threads=${task.cpus} fanc*.hic "fanc_merged.hic"
    """
}

// Convert to fanc file
process fanc_normalize {
    label "fanc"
    tag "${samplename}:${norm}"
    input:
        tuple val(samplename), path(fanc_file)
        each norm
    output:
        tuple val(samplename), val(norm), path("fanc_*_normalized.hic") 
    script:
    """
    TMPDIR="\$PWD"
    fanc hic --normalise --norm-method=${norm} $fanc_file "fanc_${norm}_normalized.hic"
    """
}

// compare the two hic_matrices using fanc compare
process diff_hic_matrices {
    label "publish"
    label "fanc"
    tag "${meta.id}:${norm}"
    input:
        tuple val(meta), val(norm), path("sample1.fanc_normalized.hic"), path("sample2.fanc_normalized.hic") 
    output:
        tuple val(meta), val(norm), path("*.hic") 
    script:
    """
    TMPDIR="\$PWD"
    # note: the -u option overrides the @norm notation. So, including it will only ever work on the uncorrected matrix
    fanc compare -c difference sample1.fanc_normalized.hic sample2.fanc_normalized.hic ${meta.sample1}_minus_${meta.sample2}_${params.hic_resolution}_${norm}.hic
    """
}

// compare the two hic_matrices using fanc compare
process foldchange_hic_matrices {
    label "publish"
    label "fanc"
    tag "${meta.id}:${norm}"
    input:
        tuple val(meta), val(norm), path("sample1.fanc_normalized.hic"), path("sample2.fanc_normalized.hic") 
    output:
        tuple val(meta), val(norm), path("*.hic") 
    script:
    """
    TMPDIR="\$PWD"
    # note: the -u option overrides the @norm notation. So, including it will only ever work on the uncorrected matrix
    fanc compare --comparison "fold-change" --log --ignore-zero --ignore-infinite sample1.fanc_normalized.hic sample2.fanc_normalized.hic ${meta.sample1}_div_${meta.sample2}_${params.hic_resolution}_${norm}.hic
    """
}

// Plot the difference matrix
process plot_difference_matrix {
    label "publish"
    label "fanc"
    tag "${meta.id}:${norm}:${plotting_region}"
    input:
        tuple val(meta), val(norm), path(diff_mtx) 
        each plotting_region
    output:
        tuple val(meta), path("*.png") 
    script:
    """
    fancplot "${plotting_region}" -o "${meta.sample1}_minus_${meta.sample2}_${params.hic_resolution}_${norm}_${plotting_region}.png" -p square -s=0 -c="coolwarm" "${diff_mtx}"
    """
}

// plot the fold change matrix
process plot_foldchange_matrix {
    label "publish"
    label "fanc"
    tag "${meta.id}:${norm}:${plotting_region}"
    input:
        tuple val(meta), val(norm), path(foldchange_mtx) 
        each plotting_region
    output:
        tuple val(meta), path("*.png") 
    script:
    """
    fancplot "${plotting_region}" -o "${meta.sample1}_div_${meta.sample2}_${params.hic_resolution}_${norm}_${plotting_region}.png" -p square -s=0 -vmin=-2 -vmax=2 -c="coolwarm" "${foldchange_mtx}"
    """
}
