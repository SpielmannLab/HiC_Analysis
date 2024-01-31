// launch an srun (slurm) job and go to $SCRATCH
// module load nextflow/v22.04.01
// submit example: nextflow run -entry callPeaks peakachu.nf -params-file peakachu.yml --id ${SLURM_JOB_ID}
// submit example: nextflow run -entry diffPeaks peakachu.nf -params-file peakachu.yml --id ${SLURM_JOB_ID}

// Use dsl 2 to separate process declaration and workflow
nextflow.enable.dsl=2

// **** Load modules ****
// for diffpeakachu analysis
include { get_union_of_loops as get_union_of_group1_loops } from "./diff_loops_modules.nf"
include { get_union_of_loops as get_union_of_group2_loops } from "./diff_loops_modules.nf"
include { do_diffpeakachu } from "./diff_loops_modules.nf"
include { annotate_loops as annotate_unique_loops1 } from "./diff_loops_modules.nf"
include { annotate_loops as annotate_unique_loops2 } from "./diff_loops_modules.nf"
include { annotate_loops as annotate_merged_loops1 } from "./diff_loops_modules.nf"
include { annotate_loops as annotate_merged_loops2 } from "./diff_loops_modules.nf"

/* ------------ PEACKACHU WORKFLOW --------------------- */
workflow {
  // Get union of filtered loops across replicates for the two groups
  filtered1 = Channel.fromList(params.group1)
    | map { it -> 
      meta = it.subMap('groupname')
        [ meta, it.filtered_loops ]
    }
    | groupTuple

  filtered2 = Channel.fromList(params.group2)
    | map { it -> 
      meta = it.subMap('groupname')
        [ meta, it.filtered_loops ]
    }
    | groupTuple

  // Merge the filtered loops from the replicates per group - ignoring the probabilites
  get_union_of_group1_loops(filtered1, params.hic_resolution*2.5) // Setting the merge resolution to 2.5x does a pretty good job
  get_union_of_group2_loops(filtered2, params.hic_resolution*2.5)

  // ---- Do diffPeakachu by pair probabilities using modified python scripts
  unfiltered1 = Channel.fromList(params.group1)
    | map { it -> 
      meta = it.subMap('groupname')
        [ meta, it.unfiltered_loops ]
    }
    | groupTuple

  unfiltered2 = Channel.fromList(params.group2)
    | map { it -> 
      meta = it.subMap('groupname')
        [ meta, it.unfiltered_loops ]
    }
    | groupTuple

  do_diffpeakachu(get_union_of_group1_loops.out.loop_format,
     get_union_of_group2_loops.out.loop_format,
     unfiltered1,
     unfiltered2)

  
  // Annotate the loops and extract the unique genes in each of the two groups
  annotate_unique_loops1(do_diffpeakachu.out.unique_loops1, params.gene_annotation_bedfile)
  annotate_unique_loops2(do_diffpeakachu.out.unique_loops2, params.gene_annotation_bedfile)
  annotate_merged_loops1(get_union_of_group1_loops.out.loop_format, params.gene_annotation_bedfile)
  annotate_merged_loops2(get_union_of_group2_loops.out.loop_format, params.gene_annotation_bedfile)
  /*
   */
}

/* ------------ HICCUPS WORKFLOW --------------------- */

