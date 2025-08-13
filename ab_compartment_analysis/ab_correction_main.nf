// launch an srun (slurm) job and go to $SCRATCH
// module load nextflow/v22.04.01
// submit example: nextflow run ab_correction_main.nf -params-file ab_correction_params.yaml -resume -ansi-log false

// Use dsl 2 to separate process declaration and workflow
nextflow.enable.dsl = 2

// **** Load modules ****
include { write_params } from './modules.nf'
include { convert_bigWig_to_bedGraph ; call_peaks_on_chipseq_data_with_cutoff_analysis ; call_peaks_on_chipseq_data } from './modules.nf'
include { align_reference_ev_to_chip_peaks; align_query_ev_and_domains_to_ref_ev } from './modules.nf'
include { enrichmentplot_with_corrected_ev; plot_compartments_with_corrected_ev } from './modules.nf'

workflow {

    write_params()

    /* Collect all the samples and their corresponding files */
    samples = channel.fromList(params.samples)
        | map { it ->
            def meta = it.subMap("samplename", "prefix")
            [meta, params.hic_files_glob, params.domain_files_glob, params.ev_files_glob, params.corr_mat_files_glob]
        }
        | map { it ->
            def all_hic_files = file(it[1])
            def select_hic_file = all_hic_files.find { s -> "${s.name}" =~ it[0].prefix }
            def all_domain_files = file(it[2])
            def select_domain_file = all_domain_files.find { s -> "${s.name}" =~ it[0].prefix }
            def all_ev_files = file(it[3])
            def select_ev_file = all_ev_files.find { s -> "${s.name}" =~ it[0].prefix }
            def all_corr_mat_files = file(it[4])
            def select_corr_mat_file = all_corr_mat_files.find { s -> "${s.name}" =~ it[0].prefix }
            [it[0], select_hic_file, select_domain_file, select_ev_file, select_corr_mat_file]
        }
        | view

    /* Align eigen vectors so that A compartments are the ones with highers activity in functional data */
    // Get one eigen vector file to be used as reference. This will be directly aligned to functional data.
    ref_ev = samples
        | filter { it[0].samplename == params.ref_ev_sample }
        | map { it ->
            println("Will use ${it[0].samplename} and its eigen vector file, ${it[3].name}, as reference to align all other eigen vectors")
            it[3]
        }
    // Convert chipseq data file format if needed
    if (file(params.chip_data).extension ==~ /bigWig/) {
        println("The ChIP-seq data will be converted to bedGraph format for peak-calling with MACS3")
        chipseq_bedgraph = convert_bigWig_to_bedGraph(params.chip_data)
    } else {
        chipseq_bedgraph = channel.of(params.chip_data)
    }
    // Get threshold to call peak on chipseq data
    if (params.macs3_cutoff == null) {
        call_peaks_on_chipseq_data_with_cutoff_analysis(chipseq_bedgraph)
    }
    // Call peaks on the chipseq data, and then align the ref ev first to it. Then align the evs of all samples to this aligned reference ev
    else {
        chip_peaks = call_peaks_on_chipseq_data(chipseq_bedgraph)
        align_reference_ev_to_chip_peaks(ref_ev, chip_peaks)
        to_align_query_and_aligned_ref = samples
        | combine(align_reference_ev_to_chip_peaks.out.bedfile)
        align_query_ev_and_domains_to_ref_ev(to_align_query_and_aligned_ref)
        // Now make plots
        plotting_regions = channel.fromList(params.plotting_regions)
        enrichmentplot_with_corrected_ev(align_query_ev_and_domains_to_ref_ev.out.corrected_ch)
        plot_compartments_with_corrected_ev(align_query_ev_and_domains_to_ref_ev.out.corrected_ch, plotting_regions)
    }
}
