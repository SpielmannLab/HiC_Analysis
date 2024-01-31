process build_reference_indices {
    // This process can help create genome.fa.fai file and bowtie indices for a genome not downloadable from iGenome
    tag "${reference_name}"
    storeDir "$WORK/hic_built_references_nextflow"
    input:
        tuple val(reference_name), path(fa_file)
    output:
        tuple val(reference_name), path(fa_file), path("${fa_file}.fai"), path("${reference_name}_bowtie2_index")
    shell:
        """
        # create samtools index
        samtools faidx ${fa_file}
        
        # Make a directory to make bowtie2 index
        mkdir ${reference_name}_bowtie2_index
        bowtie2-build $fa_file ${reference_name}_bowtie2_index/${reference_name}
        """
}

process create_restriction_fragments {
    // HiC pair filtering requies the knoweldge of fragments possible by enzymatic digestion
    // Place the output in the fixed temporary_path accessible by all compute nodes
    publishDir "${params.temporary_path}", mode: 'copy'
    input:
        path "HiC-Pro"
        path "genome.fa"
	output:
        path "restriction_fragments.bed"
	shell:
        """
        HiC-Pro/bin/utils/digest_genome.py -r ${params.restriction_site_or_enzyme} -o restriction_fragments.bed genome.fa
        """
}

process rename_n_stage_bowtie2index {
    // Place the bowtie index file in a location accessible by all compute nodes. 
    // Also need a fixed full path to insert to the config file.
    publishDir "${params.temporary_path}", mode: 'copy'
    input:
        path "Bowtie2Index"
	output:
        path "Bowtie2Index/${params.genome_reference.name}*"
	shell:
        '''
        cd Bowtie2Index
        for file in *
            do
            mv ${file} !{params.genome_reference.name}.${file#*.}
        done
        '''
}

process create_genome_sizes {
    // Create a genome.sizes file from the genome.fa.fai file
    // Place the sizes file in a fixed path accessible by all compute nodes.
    publishDir "${params.temporary_path}", mode: 'copy'
    input:
        path "genome.fa.fai"
	output:
        path "genome.sizes"
	shell:
        """
        # Collect the sizes of all chromosomal scaffolds. Get rid of others and EBV (decoy)
        cut -f1-2 genome.fa.fai | grep -v "_" | grep -v EBV > genome.sizes
        """
}

process downsample_fastqgz {
    // downsample the fastqfile for testing
    // this can take nearly an hour for 70GB fastq.gz files
    tag "${samplename}"
    input:
        tuple val(samplename), path(fastqfiles)
    output:
        tuple val(samplename), path('*.fastq.gz')
    shell:
        '''
        for f in *.fastq.gz
		do
			zcat $f | head -n !{params.downsample_to} > downsampled_${f/.fastq.gz/.fastq}
            rm $f
		done
        gzip *
        '''
}

process split_fastqgz {
    tag "${samplename}"
    // splitting the fastqgz file can speed up the HiC-Pro mapping process by parallelization
    input:
        tuple val(samplename), path(fastqfiles)
        path "HiC-Pro"
    output:
        tuple val(samplename), path('*.fastq')
    shell:
        '''
        # This python script from HiC-Pro creates split fastq files
        for f in *.fastq.gz
		do
			python HiC-Pro/bin/utils/split_reads.py $f
		done
        '''
}

process edit_config_file {
    // HiC-Pro requires its own configuration file, config-hicpro.txt to create batch script. 
    // This file will be updated with certain parameters from the *_params.yaml file, so the user does not have to
    // The path to bowtie2indices, genome fragments and genome size files need to be absolute paths in the config-hicpro.txt file
    // Hence, the temporary path is chosen.
    input:
        path "restriction_fragments.bed" // added as dummy to ensure process order
        path "genome.sizes" // added as dummy to ensure process order
        path bowtie2index
        path "config-hicpro.txt"
    output:
        path "config-hicpro.txt", includeInputs: true
	shell:
        '''
        sed --in-place "s/LIGATION_SITE = .*/LIGATION_SITE = !{params.ligation_site}/" config-hicpro.txt
        sed --in-place "s@BOWTIE2_IDX_PATH = .*@BOWTIE2_IDX_PATH = !{params.temporary_path}Bowtie2Index@" config-hicpro.txt
        sed --in-place "s@GENOME_FRAGMENT =.*@GENOME_FRAGMENT = !{params.temporary_path}restriction_fragments.bed@" config-hicpro.txt
        sed --in-place "s@REFERENCE_GENOME =.*@REFERENCE_GENOME = !{params.genome_reference.name}@" config-hicpro.txt
        sed --in-place "s@GENOME_SIZE =.*@GENOME_SIZE = !{params.temporary_path}genome.sizes@" config-hicpro.txt
        '''
}

process hic_pro_in_parallel_mode {
    tag "${samplename}"
    // Start the HiC-Pro analysis to create ValidPairs and ICE-normalized matrices from fastq files
    publishDir "${params.outdir}/${samplename}", mode: 'copy'
    input:
        path "HiC-Pro"
        path "config-hicpro.txt"
        tuple val(samplename), path(fastqfiles, stageAs: "fastqfiles/sample/*")
    output:
        path "results/hic_results/pic/${samplename}/*"
        path "results/hic_results/stats/${samplename}/*"
        path "results/hic_results/pic"
        path "results/hic_results/matrix"
        path "results/logs"
        path "results/config-hicpro.txt"
        tuple val(samplename), path("results/hic_results/data/${samplename}/${samplename}.allValidPairs"), emit: allValidPairs
	shell:
        '''
        # rename the same dir to get correct output file names
        mv fastqfiles/sample fastqfiles/!{samplename}

        HiC-Pro/bin/HiC-Pro -i fastqfiles -o results -c config-hicpro.txt -p
        
        # Start the step1 script created by the parallel processing above note the "-p" in the above line.
        cd results
        # This first only needs 7GB of memorywhen fastq files are split with 20M reads
        sbatch -J !{samplename}_step1 --wait --mem-per-cpu=7GB HiCPro_step1_HiC.sh

        # This second step has a higher memory requirement. For 25 Billion read pairs 100GB was needed.
        sbatch -J !{samplename}_step2 --wait --mem-per-cpu=100GB HiCPro_step2_HiC.sh
        '''
}

process allValidPairs_to_juicertools_hic {
    tag "${samplename}_${juicername}"
    debug true
    // create hic file from *.allValidPairs file created by HiC-Pro
    publishDir "${params.outdir}/${samplename}/${juicername}", mode: 'copy'
    input:
        tuple val(samplename), path(allValidPairs), val(juicername), path(juicerpath)
        path "genome.sizes"
        path "HiC-Pro"
	output:
        path "*.hic"
	shell:
        '''
        HiC-Pro/bin/utils/hicpro2juicebox.sh -i !{allValidPairs} -g genome.sizes -j !{juicerpath} -o "./" 
        '''
}
