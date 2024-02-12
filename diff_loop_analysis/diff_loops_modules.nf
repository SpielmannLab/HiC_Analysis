/* 
This does two things:
1. Gather the max probababilites within each group for each loop identified by the process merge_loop_files
2. Find unique loops - either manually, or by using Gausian Mixture Model
Both scripts adapted from diffPeakachu
*/

// Get the union of all the loops passed on as a list of loop_files. Merge neighbours using slop_by in bp
process get_union_of_loops {
    tag "${meta.groupname}"
    input:
        tuple val(meta), path(loop_files)
        val slop_by
    output:
        tuple val(meta), path("*merged.loops"), emit: loop_format 
        tuple val(meta), path("*merged_forIGV.bedpe"), emit: bedpe_format
    shell:
    '''
    files=(!{loop_files}) # convert to bash array
    cat ${files[0]} > !{meta.groupname}_merged.loops # Transfer contents of first file to the output file for accumulation
    unset files[0] # remove first file from array

    # Merge the loops without duplication
    for file in "${files[@]}"; do
        bedtools pairtopair -is -slop !{slop_by} -type notboth -a $file -b !{meta.groupname}_merged.loops > tmp.loops
        cat tmp.loops >> !{meta.groupname}_merged.loops
        done

    # Remove score fields to ensure proper display in IGV
    cat !{meta.groupname}_merged.loops | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6}' > !{meta.groupname}_merged_forIGV.bedpe 
    '''
}

process do_diffpeakachu {
    tag "${meta_filtered_merged_1.groupname}"+"vs"+"${meta_filtered_merged_2.groupname}"
    input:
        tuple val(meta_filtered_merged_1), path(union_of_loops1)
        tuple val(meta_filtered_merged_2), path(union_of_loops2)
        tuple val(meta_unfiltered_1), path("group1_unfiltered_??.loops")
        tuple val(meta_unfiltered_2), path("group2_unfiltered_??.loops")
    output:
        path "${meta_filtered_merged_1.groupname}-${meta_filtered_merged_2.groupname}.merged.loops"
        tuple val(meta_filtered_merged_1), path("*${meta_filtered_merged_1.groupname}.unique.loops"), emit: unique_loops1
        tuple val(meta_filtered_merged_2), path("*${meta_filtered_merged_2.groupname}.unique.loops"), emit: unique_loops2
        path "*.png"
    shell:
    '''
    # Combine the probabilities from the repeats, staged as file_group{1,2}_*. Uses the maximum of the available repeats
    echo "Running pair-probs.py to gather max probabilities from the repeats from the two groups"
    pair-probs.py \
        !{union_of_loops1} \
        !{union_of_loops2} \
        "group1_unfiltered_" \
        "group2_unfiltered_" \
        !{meta_filtered_merged_1.groupname}-!{meta_filtered_merged_2.groupname}.merged.loops
    echo "pair-probs.py output:"
    head !{meta_filtered_merged_1.groupname}-!{meta_filtered_merged_2.groupname}.merged.loops

    # Calculate the diffPeakachu scores
    diffPeakachu.py !{union_of_loops1} !{union_of_loops2} !{meta_filtered_merged_1.groupname}-!{meta_filtered_merged_2.groupname}.merged.loops
    '''
}

process annotate_loops {
	  input:
        tuple val(meta), path(loops)
        path "genes.bed"
    output:
        path("overlapping_genes_in*")
    shell:
        '''
        cat !{loops} | grep -v ^# | awk 'BEGIN { OFS="\t" } {print $1, $2, $3, $7, $8}' > !{loops}_pe1.bed
        cat !{loops} | grep -v ^# | awk 'BEGIN { OFS="\t" } {print $4, $5, $6, $7, $8}' > !{loops}_pe2.bed
        bedtools intersect -wa -wb -a genes.bed -b !{loops}_pe1.bed !{loops}_pe2.bed > overlapping_genes_in_!{loops}
        '''
}
