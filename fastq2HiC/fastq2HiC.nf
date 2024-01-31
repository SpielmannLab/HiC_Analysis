// Try to use slurm executor mode built into nextflow. So, spawn jobs as needed

// launch an srun (slurm) job and go to $SCRATCH
// module load nextflow/v22.04.01
// submit example: nextflow run -entry callPeaks peakachu.nf -params-file peakachu.yml --id ${SLURM_JOB_ID}
// submit example: nextflow run -entry diffPeaks peakachu.nf -params-file peakachu.yml --id ${SLURM_JOB_ID}

// **** Load modules ****
include { create_restriction_fragments; split_fastqgz; downsample_fastqgz; edit_config_file } from "./fastq2HiC_modules.nf"
include { rename_n_stage_bowtie2index } from "./fastq2HiC_modules.nf"
include { create_genome_sizes; hic_pro_in_parallel_mode; allValidPairs_to_juicertools_hic } from "./fastq2HiC_modules.nf"

// **** Real stuff happens here ****
workflow {
    // Get all the reference files in order
    create_restriction_fragments(params.hicpro_path, params.genome_reference.fa_file)
    create_genome_sizes(params.genome_reference.fa_fai_file)
    rename_n_stage_bowtie2index(params.genome_reference.bowtie2index)

    // Get the input fastq files from the provided details in the params file
    input_fastq_ch = Channel.from(params.input_fastqfiles)
    | map { it ->
        files = files(it.fastq_files)
        [it.samplename, files]
        }
    | transpose() // This operator splits the tuple into individual channel, while preserving the samplename metadata. Enables parallelization
    | view()
        
    // based on the request to downsample, either downsample or split to chunks for parallel processing
    if (params.downsample) {
        // downsample the fastqfile
        input_fastq_ch = downsample_fastqgz(input_fastq_ch)
    } 

    // split the fastq files up
    input_fastq_ch = split_fastqgz(input_fastq_ch, params.hicpro_path)
    | groupTuple() // This works as the opposite of transpose, collecting all fastq files belonging to a samplename
    | map { id, list -> 
        [id, list.flatten()] // This is needed, because the spliitting results in a list of files that need to be flattened
    }


    // Edit the config-hicpro.txt file with parameters. The only thing not implemented is the bowtie index names
    // While, the edit_config_file process does not really need the outs from the three process, it is added to make sure they have run
    config_ch = edit_config_file(create_restriction_fragments.out,
        create_genome_sizes.out,
        rename_n_stage_bowtie2index.out,
        projectDir+"/config-hicpro.txt")

    // run hicpro in parallel mode using split fastq files
    hic_pro_in_parallel_mode(params.hicpro_path,
        config_ch,
        input_fastq_ch)

    allValidPairs_ch = hic_pro_in_parallel_mode.out.allValidPairs
    | map { samplename, allvalidpairs_file ->
        [samplename:samplename, allvalidpairs_file:allvalidpairs_file]}

    juicer_ch = Channel.fromList(params.juicer_path)

    allValidPairs_juicer_ch = allValidPairs_ch
    | combine(juicer_ch)
    | map { it ->
        [samplename:it[0].samplename, allvalidpairs_file:it[0].allvalidpairs_file, juicername:it[1].name, juicerpath:it[1].path]
    }
    | view

    allValidPairs_to_juicertools_hic(allValidPairs_juicer_ch,
        create_genome_sizes.out,
        params.hicpro_path)
}

// ################################################
// **************** SUB WORKFLOWS *****************
// ################################################

workflow from_allValidPairs { // Run from allValidPairs to create a juicebox *.hic file

    // Get the necessary reference files
    create_genome_sizes(params.genome_reference.fa_fai_file)

    juicer_ch = Channel.fromList(params.juicer_path)
    | view

    allValidPairs_ch = Channel.from(params.input_allvalidpairsfile)
    | view()
    
    allValidPairs_juicer_ch = allValidPairs_ch
    | combine(juicer_ch)
    | map { it ->
        [samplename:it[0].samplename, allvalidpairs_file:it[0].allvalidpairs_file, juicername:it[1].name, juicerpath:it[1].path]
    }
    | view

    allValidPairs_to_juicertools_hic(allValidPairs_juicer_ch,
        create_genome_sizes.out,
        params.hicpro_path)
}

// If a reference is not available from iGenomes, use this worfklow to build samtools and bowtie reference indices 
// The output is a tuple with a name, fa, fai.fai and bowtie2index folder
workflow build_reference {
    take: genome_ch
    main:
        genome_ch
        | map {it -> it.subMap('name', 'fa_file')}
        | build_reference_indices
        | set {built_reference}
    emit: built_reference
}

