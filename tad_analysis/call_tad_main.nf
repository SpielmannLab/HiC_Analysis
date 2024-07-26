// **** Load modules ****
include { write_params } from "./call_tad_modules.nf"
include { perform_arrowhead } from "./call_tad_modules.nf"
include { perform_diffDomain as perform_diffDomain_1vs2; perform_diffDomain as perform_diffDomain_2vs1 } from "./call_tad_modules.nf"
include { visualize_diffDomain } from "./call_tad_modules.nf"

/* ------------ WORKFLOW --------------------- */
workflow {
    write_params()

    samples = Channel.fromList(params.samples)
        | map { it -> 
            meta = it.subMap('samplename', 'group')
                [ meta, it.hic_file ]
        }

    // Peform arrowhead if asked for
    perform_arrowhead(samples, params.juicer_path)

    if (params.do_differential_analysis) {

        // Split the output into two groups to peform diffDomain
        group1_channel = perform_arrowhead.out.arrowhead_out
            | filter { it[0].group == "group1" }
            | view { it -> "Control sample (group1) is ${it[0].samplename}" }
        group2_channel = perform_arrowhead.out.arrowhead_out
            | filter { it[0].group == "group2" }
            | view { it -> "Test sample (group2) is ${it[0].samplename}" }

        // perfrom diffDomain
        // perform comparision Group1 vs Group2
        perform_diffDomain_1vs2(group1_channel, group2_channel, params.diffDomain_path)
        // perform comparision Group2 vs Group1
        perform_diffDomain_2vs1(group2_channel, group1_channel, params.diffDomain_path)

        // visualize diffDomains from both pairwise comparisons. 
        // Join them into one channel and then run visualization process
        diffdomains_channel = perform_diffDomain_1vs2.out
            | mix(perform_diffDomain_2vs1.out)
            | view
        visualize_diffDomain(diffdomains_channel, params.diffDomain_path)
    }
}
