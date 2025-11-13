// **** Load modules ****
include { write_params } from "./apa_modules.nf"
include { do_juicer_apa } from "./apa_modules.nf"

/* ------------ WORKFLOW --------------------- */
workflow {
    apa_comparisons = Channel.fromList(params.apa_groups)
        | map { it -> 
            meta = it.subMap('id')
                [ meta, it.hic_file, it.unique_peaks ]
        }

    write_params()

    // Peform aggregage peak analysis
    do_juicer_apa(apa_comparisons, params.juicer_path)
}
