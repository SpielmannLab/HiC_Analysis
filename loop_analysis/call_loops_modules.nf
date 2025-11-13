// Modules for peakachu and diffpeakachu.nf

process addnorm {
/* Add normalize to juicebox hic file, if needed. It should already be there after the fastq2HiC pipeline.
Just in case you want a particular normalization, add it.
This process is not in the main nextflow pipeline. Needs to be custom integrated as needed
*/
    input:
        tuple val(sample_name), path(path_hic), path(peakachu_model), val(threshold)
        path juicerpath
    output:
        tuple val(sample_name), path("*.hic", includeInputs: true), path(peakachu_model), val(threshold)
    shell:
        '''
        java -jar !{juicerpath} addNorm -r !{params.hic_resolution} -k !{params.normalization} -j !{task.cpus} !{path_hic}
        '''
}

// Run peakachu2 on all the provided HiC files
// First score peaks across the genome
process peakachu_score_genome {
    tag "${sample_name}"
    input:
        tuple val(sample_name), path(path_hic), path(peakachu_model), val(threshold)
	output:
        tuple val(sample_name), path("*_unfiltered.loops"), val(threshold)
	shell:
        '''
        tree
        hicfile=!{path_hic} # conversion of nextflow variable to bash variable to enable substitution
        # the balance option uses ICE/KR-balanced matrix. But for whatever reason, this does not work 
        # peakachu score_genome --balance --minimum-prob 0 -p ${hicfile} -r !{params.hic_resolution} -O !{sample_name}_unfiltered.loops -m !{peakachu_model}
        peakachu score_genome --minimum-prob 0 -p ${hicfile} -r !{params.hic_resolution} -l !{params.min_anchor_distance} -u !{params.max_anchor_distance} -O !{sample_name}_unfiltered.loops -m !{peakachu_model}
        '''
}

// Pool the peaks
process peakachu_pool {
    tag "${sample_name}"
    input:
        tuple val(sample_name), path(unfiltered_loops), val(threshold)
	output:
		path "*_filtered.loops", emit: filtered_loops
	shell:
        '''
        peakachu pool -i !{sample_name}_unfiltered.loops -t !{threshold}  -o !{sample_name}_filtered.loops -r !{params.hic_resolution}
        '''
}


/*
process hiccups_loop_enrichment {
    shell:
    '''
    /data/humangen_external/test_area/varun/trying_hiccups
    java -jar /data/humangen_external/HiC/installation/Juicer/juicer_tools_1.22.01.jar hiccups -m 1024 -r 10000 -c 1 -k NONE /data/humangen_external/HiC/results/hicpro/DNMT3Amut55/hic_format/DNMT3Amut55_hg38.hic "./" DNMT3A-WT.merged.hiccups_loop_list
    '''
}

// Run Moustache on all the provided HiC files
process moustache_score_genome {
    tag "${sample_name}"
    input:
        tuple val(sample_name), path(path_hic), path(peakachu_model), val(threshold)
	output:
        tuple val(sample_name), path("*_unfiltered.loops"), val(threshold)
	shell:
        '''
        tree
        hicfile=!{path_hic} # conversion of nextflow variable to bash variable to enable substitution
        # the balance option uses ICE/KR-balanced matrix. But for whatever reason, this does not work 
        # peakachu score_genome --balance --minimum-prob 0 -p ${hicfile} -r !{params.hic_resolution} -O !{sample_name}_unfiltered.loops -m !{peakachu_model}
        peakachu score_genome --minimum-prob 0 -p ${hicfile} -r !{params.hic_resolution} -l !{params.min_anchor_distance} -u !{params.max_anchor_distance} -O !{sample_name}_unfiltered.loops -m !{peakachu_model}
        '''
}
*/
