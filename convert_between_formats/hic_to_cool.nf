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
            meta = it.subMap('samplename', 'resolution' )
            [ meta, it.hic_file ] }
    | set { samples }

    // Convert to cool format 
    hic_to_cool_raw(samples)
    // Normalize (all samples together) and correct (KR, ICE etc).
    cool_normalizeNcorrect(hic_to_cool_raw.out.cool_files)

}
