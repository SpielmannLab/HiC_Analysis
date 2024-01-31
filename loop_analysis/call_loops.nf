// launch an srun (slurm) job and go to $SCRATCH
// module load nextflow/v22.04.01
// submit example: nextflow run -entry callPeaks peakachu.nf -params-file peakachu.yml --id ${SLURM_JOB_ID}
// submit example: nextflow run -entry diffPeaks peakachu.nf -params-file peakachu.yml --id ${SLURM_JOB_ID}

// Use dsl 2 to separate process declaration and workflow
nextflow.enable.dsl=2

// **** Load modules ****
// for peakachu analysis
include { peakachu_score_genome } from "./call_loops_modules.nf"
include { peakachu_pool } from "./call_loops_modules.nf"

/* ------------ PEACKACHU WORKFLOW --------------------- */
workflow {

    // add normalization to hic
    hic_channel = Channel.fromList(params.samples)

    // call loops using score_genome on the *hic file and also pool them (i.e., cluster the neighbouring loop pixels)
    peakachu_score_genome(hic_channel)
    | peakachu_pool

}

/* ------------ HICCUPS WORKFLOW --------------------- */

