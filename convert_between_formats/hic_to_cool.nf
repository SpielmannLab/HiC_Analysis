// Use dsl 2 to separate process declaration and workflow
nextflow.enable.dsl=2

// **** Load modules ****
include { write_params } from "./hic_to_cool_modules.nf"
include { hic_to_cool_raw } from "./hic_to_cool_modules.nf"
include { cool_normalizeNcorrect } from "./hic_to_cool_modules.nf"

/* ------------ Convert of Juicer HiC to *.cool WORKFLOW --------------------- */
workflow {

    // write params
    write_params()

    // add normalization to hic
    Channel.fromList(params.samples)
    | map { it -> 
            meta = it.subMap('samplename', 'resolution', 'hicexplorer_normalization', 'hicexplorer_correction_method', 'hicexplorer_threshold_low', 'hicexplorer_threshold_high')
            [ meta, it.hic_file ] }
    | set { samples }

    // Convert to cool format 
    hic_to_cool_raw(samples)
    // Normalize (0-1) and correct (KR, ICE etc).
    cool_normalizeNcorrect(hic_to_cool_raw.out.samples)

}
