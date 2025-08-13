// write parameters to file
process write_params {
    output:
    path "*.txt"

    script:
    """
    suffix=\$(date "+%Y_%m_%d_%H_%M")
    echo "Parameters of the job" >> parameters_ab_compartments_\${suffix}.txt
    echo "${params}" | tr , '\n' >> parameters_ab_compartments_\${suffix}.txt
    """
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

// Downsample fanc hic matrices
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

// Merge multiple fanc hic matrices
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

// Normalize hic matrices as per the method 
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
    fanc hic --normalise --norm-method=${norm} ${fanc_file} "fanc_${norm}_normalized.hic"
    """
}

// Get ab compartments, calculate eigen vector, plot enrichment score etc.
process get_compartments {
    label "publish"
    label "fanc"
    tag "${meta.samplename}:${norm}"

    input:
    tuple val(meta), val(norm), path("fanc.hic"), val(eigen_vector_index)

    output:
    tuple val(meta), val(norm), path("fanc.ab"), path("*.bed"), emit: compartments
    path "*.png"
    path "*.tsv"

    script:
    """
    TMPDIR="\$PWD"
    fanc compartments -d="domains.bed" -v="eigen_vectors.bed" -e="enrichments.png" -m="enrichments.tsv" -g="${params.genome_fa}" -i="${eigen_vector_index}" -s=0 "fanc.hic" "fanc.ab"
    """
}

// Plot the compartments and eigen vectors
process plot_compartments {
    label "publish"
    label "fanc"
    tag "${meta.samplename}:${norm}:${plotting_region}"

    input:
    tuple val(meta), val(norm), path("fanc.ab"), path(bedfiles)
    each plotting_region

    output:
    tuple val(meta), path("*.png")

    script:
    """
    fancplot "${plotting_region}" -o="compartments.png" -p square fanc.ab -p line eigen_vectors.bed
    """
}

process call_peaks_on_chipseq_data_with_cutoff_analysis {
    label "macs3"

    input:
    path "chipseq.bedGraph"

    output:
    path "chippeaks_cutoff_analysis.png"
    path "chippeaks_cutoff_analysis.txt"

    script:
    """
    macs3 bdgpeakcall -i chipseq.bedgraph --cutoff-analysis --ofile chippeaks_cutoff_analysis.txt
    plot_macs_cutoff_analysis.py
    """
}

process call_peaks_on_chipseq_data {
    label "macs3"
    tag "macs3 cutoff: ${params.macs3_cutoff}"

    input:
    path "chipseq.bedGraph"

    output:
    path "chippeaks.bed"

    script:
    """
    macs3 bdgpeakcall -i chipseq.bedGraph --cutoff ${params.macs3_cutoff} --ofile chippeaks.bed
    """
}

process convert_bigWig_to_bedGraph {
    label "genomics"

    input:
    path "chipseq.bigWig"

    output:
    path "chipseq.bedGraph"

    script:
    """
    bigWigToBedGraph chipseq.bigWig chipseq.bedGraph
    """
}

// Align the reference eigen vector using functional data
process align_reference_ev_to_chip_peaks {
    label "publish"
    label "genomics"

    input:
    path "ref_ev.bed"
    path "chippeaks.bed"

    output:
    path "ref_ev_corrected.bed", emit: bedfile
    path "*_stats.bed"
    script:
    """
    ls
    # Do first pass on uncorrected ref_ev.bed
    sed 's/^chr//' chippeaks.bed \
      | tail -n +2 \
      | bedtools intersect -c -a ref_ev.bed -b - \
      > peaks_in_bins.txt
    echo "---------------------" 1> ref_ev_correction_stats.bed
    echo "Stats on original eigen vector" 1>> ref_ev_correction_stats.bed 
    echo "---------------------" 1>> ref_ev_correction_stats.bed 
    align_reference_ev_to_chip_peaks.awk pass=1 peaks_in_bins.txt pass=2 ref_ev.bed 1> ref_ev_corrected.bed 2>> ref_ev_correction_stats.bed

    # Then do second pass on the ref_ev_corrected.bed
    sed 's/^chr//' chippeaks.bed \
      | tail -n +2 \
      | bedtools intersect -c -a ref_ev_corrected.bed -b - \
      > peaks_in_corrected_bins.txt
    echo "---------------------" 1> ref_ev_corrected_stats.bed
    echo "Stats on corrected eigen vector" 1>> ref_ev_corrected_stats.bed
    echo "---------------------" 1>> ref_ev_corrected_stats.bed
    align_reference_ev_to_chip_peaks.awk pass=1 peaks_in_corrected_bins.txt pass=2 ref_ev_corrected.bed 1> deleteme.bed 2>> ref_ev_corrected_stats.bed
    rm deleteme.bed
    """
}

//Align all eigen vectors and domains to a reference eigen vector that has been oriented using functional data.
process align_query_ev_and_domains_to_ref_ev {
    label "publish"
    label "genomics"
    tag "${meta.samplename}"

    input:
    tuple val(meta), path(hic_file), path("query_domains.bed"), path("query_ev.bed"), path(corr_mat), path("ref_ev.bed")

    output:
    tuple val(meta), path(hic_file), path("query_domains_corrected.bed"), path("query_ev_corrected.bed"), path(corr_mat), emit: corrected_ch
    path "*.txt"
    path "uncorrected*"

    script:
    """
    bedtools intersect -wa -wb -a ref_ev.bed -b query_ev.bed \
    | cut -f 1-5,10,11 > ./ref_query.tmp

    # First pass to check query_ev alignment to ref_ev
    # Second pass to flip the query_ev as needed
    # Third pass to flip the domain assignments in the query as needed
    align_query_ev_and_domains_to_reference_ev.awk pass=1 ref_query.tmp pass=2 ref_query.tmp pass=3 query_domains.bed >> "correction_stats.txt"

    # Check if all domain entries were corrected
    bedtools intersect -v -a query_domains.bed -b query_domains_corrected.bed > "uncorrected_domains.bed"
    # If any of the domains were not aligned to the reference eigen vector, then raise error and quit
    # Exit with error if non-interactive shell
    if [ -s "uncorrected_domains.bed" ]; then
      printf "Warning: There are %d uncorrected domains.\n" \$( cat "uncorrected_domains.bed" | wc -l ) >> "correction_stats.txt"
    fi

    # Check if all eigen vector entries were corrected
    bedtools intersect -v -a query_ev.bed -b query_ev_corrected.bed > "uncorrected_ev.bed"
    # If any of the eigen vectors were not aligned to the reference eigen vector, then raise error and quit
    # Exit with error if non-interactive shell
    if [ -s "uncorrected_ev.bed" ]; then
      printf "Warning: There are %d uncorrected eigen vector values in the query eigen vector.\n" \$( cat "uncorrected_ev.bed" | wc -l ) >> "correction_stats.txt"
    fi
    """
}

process enrichmentplot_with_corrected_ev {
    label "publish"
    label "fanc"
    tag "${meta.samplename}"

    input:
    tuple val(meta), path("fanc.hic"), path(domains_corrected), path("ev.bed"), path("fanc.ab")

    output:
    path "*.png"
    path "*.txt"

    script:
    """
    plot_enrichments.py
    """
}

// Plot the compartments and eigen vectors
process plot_compartments_with_corrected_ev {
    label "publish"
    label "fanc"
    tag "${meta.samplename}:${plotting_region}"

    input:
    tuple val(meta), path("fanc.hic"), path(domains_corrected), path("ev.bed"), path("fanc.ab")
    each plotting_region

    output:
    path "*.png"

    script:
    """
    fancplot "${plotting_region}" -o="compartments.png" -p square fanc.ab -p line ev.bed
    """
}
