// launch an srun (slurm) job and go to $SCRATCH
// module load nextflow/v22.04.01
// submit example: nextflow run main.nf -params-file params.yaml -resume -ansi-log false

// Use dsl 2 to separate process declaration and workflow
nextflow.enable.dsl = 2

// **** Load modules ****
include { write_params } from './modules.nf'
include { juicer_to_fanc ; merge_fanc ; fanc_normalize ; fanc_downsample } from './modules.nf'
include { diff_hic_matrices ; foldchange_hic_matrices } from './modules.nf'
include { plot_difference_matrix ; plot_foldchange_matrix } from './modules.nf'

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
        | view
    }

    fanc_normalized = fanc_normalize(fanc, normalizations)
        | map { it ->
            ["samplename": it[0], "normalization": it[1], "fanc_file": it[2]]
        }

    // Get the comparisons as desired.
    comparisons = channel.fromList(params.comparisons)
        | combine(fanc_normalized)
        | filter { comp, fanc1 ->
            fanc1.samplename == comp.sample1
        }
        | combine(fanc_normalized)
        | filter { comp, _fanc1, fanc2 ->
            fanc2.samplename == comp.sample2
        }
        | filter { _comp, fanc1, fanc2 ->
            fanc1.normalization == fanc2.normalization
        }
        | map { comp, fanc1, fanc2 ->
            def meta = ["id": comp.id, "sample1": comp.sample1, "sample2": comp.sample2]
            [meta, fanc1.normalization, fanc1.fanc_file, fanc2.fanc_file]
        }

    foldchange_hic_matrices(comparisons)

    if (params.do_plotting_of_comparison) {
        plotting_regions = channel.fromList(params.plotting_regions)
        plot_foldchange_matrix(foldchange_hic_matrices.out, plotting_regions)
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
