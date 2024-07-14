// **** Load modules ****
include { write_params } from "./call_tad_modules.nf"
include { perform_arrowhead } from "./call_tad_modules.nf"
include { perform_diffDomain } from "./call_tad_modules.nf"
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

    // Split the output into two groups to peform diffDomain
    group1_channel = perform_arrowhead.out.arrowhead_out
        | filter { it[0].group == "group1" }
        | view { it -> "Control sample (group1) is ${it[0].samplename}" }
    group2_channel = perform_arrowhead.out.arrowhead_out
        | filter { it[0].group == "group2" }
        | view { it -> "Test sample (group2) is ${it[0].samplename}" }

    // pefrom diffDomain
    perform_diffDomain(group1_channel, group2_channel, params.diffDomain_path)
    
    // visualize diffDomains
    visualize_diffDomain(perform_diffDomain.out, params.diffDomain_path)
}
