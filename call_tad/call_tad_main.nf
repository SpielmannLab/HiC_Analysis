// **** Load modules ****
include { write_params } from "./call_tad_modules.nf"
include { perform_arrowhead } from "./call_tad_modules.nf"

/* ------------ WORKFLOW --------------------- */
workflow {
    samples = Channel.fromList(params.samples)
        | map { it -> 
            meta = it.subMap('samplename')
                [ meta, it.hic_file ]
        }

    write_params()

    // Peform arrowhead if asked for
    if (params.do_arrowhead) {
        perform_arrowhead(samples, params.juicer_path)
    }
}
