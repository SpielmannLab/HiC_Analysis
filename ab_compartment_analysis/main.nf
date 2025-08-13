// launch an srun (slurm) job and go to $SCRATCH
// module load nextflow/v22.04.01
// submit example: nextflow run main.nf -params-file params.yaml -resume -ansi-log false

// Use dsl 2 to separate process declaration and workflow
nextflow.enable.dsl = 2

// **** Load modules ****
include { write_params } from './modules.nf'
include { juicer_to_fanc ; fanc_downsample ; fanc_normalize ; merge_fanc } from './modules.nf'
include { get_compartments } from './modules.nf'
include { plot_compartments } from './modules.nf'

/* ------------ PEACKACHU WORKFLOW --------------------- */
workflow {

    write_params()

    normalizations = channel.fromList(params.normalizations)

    // get a list of samples, convert to fanc and normalize them
    samples = channel.fromList(params.samples)

    fanc = juicer_to_fanc(samples)
        | map { it ->
            ["samplename": it[0], "fanc_file": it[1]]
        }

    if (params.do_downsampling) {
        query_fanc = fanc
            | filter { it -> it.samplename != params.downsample_with }
        ref_fanc = fanc
            | filter { it -> it.samplename == params.downsample_with }

        downsampled = query_fanc
            | combine(ref_fanc)
            | map { q, r -> [q.samplename, q.fanc_file, r.samplename, r.fanc_file] }
            | fanc_downsample
            | map { it ->
                ["samplename": it[0], "fanc_file": it[1]]
            }
            | concat(ref_fanc)

        fanc = downsampled
    }

    if (params.merge_samples_before_comparison) {

        merges = channel.fromList(params.merges)

        fanc = merge_samples(fanc, merges)
            | concat(fanc)
    }

    fanc_normalized = fanc_normalize(fanc, normalizations)
        | map { it ->
            [["samplename": it[0]], it[1], it[2]]
        }

    eigenvectors = channel.fromList(params.eigenvectors)
        | map { it ->
            def meta = it.subMap("samplename")
            [meta, it.eigenvector_index_to_use]
        }

    ab_channel = fanc_normalized
        | combine(eigenvectors, by: 0)
        | view

    get_compartments(ab_channel)

    if (params.plot_correlation_matrix_and_ab_compartments) {
        plotting_regions = channel.fromList(params.plotting_regions)
        plot_compartments(get_compartments.out.compartments, plotting_regions)
    }
}

workflow merge_samples {
    take:
    fanc_samples_to_merge
    merges

    main:
    // Takes all the sampes into a single list
    fanc_list = fanc_samples_to_merge
        | toList

    merged_fanc_samples = merges
        | combine(fanc_list)
        | map { it ->
            def mrg = it[0]
            def smpls = it.subList(1, it.size())
            def select = smpls.findAll { s -> s.samplename in mrg.samples }
            def a = select.collect { k -> k.fanc_file }
            // [mrg, select]
            [mrg, a]
        }
        | view
        | merge_fanc
        | map { meta, merg ->
            ["samplename": meta.merge_name, "fanc_file": merg]
        }

    emit:
    merged_fanc_samples
}
