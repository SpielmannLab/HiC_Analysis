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
        path "genome.fa.fai"
    output:
        path("genes_in_loops*")
        path("promoters_in_loops*")
    shell:
        '''
        #### GENE based annotation
        ## Genes overlapping with the loop anchors
        # To overlap the 5' anchor with genes.bed
        cat !{loops} | grep -v ^# | awk 'BEGIN { OFS="\t" } {print $0}' > pe1.bed
        bedtools intersect -wa -wb -a pe1.bed -b genes.bed | awk 'BEGIN { OFS="\t" } {print $0,"5p anchor" }' > pe1_with_genes.bed
        # To overlap the 3' anchor with genes.bed firs rearrange and then put back
        cat !{loops} | grep -v ^# | awk 'BEGIN { OFS="\t" } {chrx=$1; startx=$2; endx=$3; chry=$4; starty=$5; endy=$6; $1=chry; $2=starty; $3=endy; $4=chrx; $5=startx; $6=endx; print $0}' > pe2.bed
        bedtools intersect -wa -wb -a pe2.bed -b genes.bed | awk 'BEGIN { OFS="\t" } {chry=$1; starty=$2; endy=$3; chrx=$4; startx=$5; endx=$6; $1=chrx; $2=startx; $3=endx; $4=chry; $5=starty; $6=endy; print $0,"3p anchor"}' > pe2_with_genes.bed
        ## Genes enclosed by loops
        # But do this only if the chrx and chry are the same
        cat !{loops} | grep -v ^# | awk 'BEGIN { OFS="\t" } $1==$4{chrx=$1; startx=$2; endx=$3; chry=$4; starty=$5; endy=$6; $1=chrx; $2=endx; $3=starty; $4=chry; $5=startx; $6=endy; print $0}' > pe3.bed
        bedtools intersect -wa -wb -a pe3.bed -b genes.bed | awk 'BEGIN { OFS="\t" } {chrx=$1; endx=$2; starty=$3; chry=$4; startx=$5; endy=$6; $1=chrx; $2=startx; $3=endx; $4=chry; $5=starty; $6=endy; print $0,"inside loop"}' > pe3_with_genes.bed

        # combine the three and sort
        cat pe1_with_genes.bed pe2_with_genes.bed pe3_with_genes.bed > pe1_pe2_pe3_with_genes.bed
        bedtools sort -i pe1_pe2_pe3_with_genes.bed > genes_in_loops!{loops}
        
        #### PROMOTER based annotation
        ## Define promoters as the 500b upstream and 100b downstream of TSS, but only for genes in main chromosome scaffolds
        ## Only main scaffolds because the alternate contigs are not matching between genes.bed and genome.sizes
        bedtools flank -l 500 -r 0 -s -i genes.bed -g genome.fa.fai | bedtools slop -l 0 -r 100 -s -g genome.fa.fai > promoters.bed

        ## Promoters overlapping with loop anchors
        # To overlap the 5' anchor with promoters.bed
        bedtools intersect -wa -wb -a pe1.bed -b promoters.bed | awk 'BEGIN { OFS="\t" } {print $0,"5p anchor" }' > pe1_with_promoters.bed
        # To overlap the 3' anchor with promoters.bed
        bedtools intersect -wa -wb -a pe2.bed -b promoters.bed | awk 'BEGIN { OFS="\t" } {chry=$1; starty=$2; endy=$3; chrx=$4; startx=$5; endx=$6; $1=chrx; $2=startx; $3=endx; $4=chry; $5=starty; $6=endy; print $0,"3p anchor"}' > pe2_with_promoters.bed
        ## promoters enclosed by loops
        # But do this only if the chrx and chry are the same
        bedtools intersect -wa -wb -a pe3.bed -b promoters.bed | awk 'BEGIN { OFS="\t" } {chrx=$1; endx=$2; starty=$3; chry=$4; startx=$5; endy=$6; $1=chrx; $2=startx; $3=endx; $4=chry; $5=starty; $6=endy; print $0,"inside loop"}' > pe3_with_promoters.bed

        # combine the two and sort
        cat pe1_with_promoters.bed pe2_with_promoters.bed pe3_with_promoters.bed > pe1_pe2_pe3_with_promoters.bed
        bedtools sort -i pe1_pe2_pe3_with_promoters.bed > promoters_in_loops!{loops}
        '''
}

/*
process get_hiccups_scores_at_loops {
    input:
        tuple varname
    output:
        tuple varname, emit: label
    shell:
        '''
        Here is a tried and tested command
        srun -p shortterm -c 1 --mem 50GB --gres=gpu:1 --pty bash

        #  -m higher of 1024 causes error
        java -jar /data/humangen_external/HiC/installation/Juicer/juicer_tools_1.19.02.jar hiccups -m 500 -r 10000 -c 21 /data/humangen_external/HiC/steinhaeuser_hic/data_n_results/DNMT3A_clones/2023.10_fastq2HiC/juicer_1.19.02/DNMT3Amut55.allValidPairs.hic "./hiccups_mutant_unique/" /data/humangen_external/HiC/steinhaeuser_hic/data_n_results/DNMT3A_clones/2024.03_diff_loop_analysis_extraannotations/diffpeakachu/Mutant-WT.Mutant.unique.loops
        java -jar /data/humangen_external/HiC/installation/Juicer/juicer_tools_1.19.02.jar hiccups -m 500 -r 10000 -c 21 /data/humangen_external/HiC/steinhaeuser_hic/data_n_results/DNMT3A_clones/2023.10_fastq2HiC/juicer_1.19.02/DNMT3Amut55.allValidPairs.hic "./hiccups_merged/" /data/humangen_external/HiC/steinhaeuser_hic/data_n_results/DNMT3A_clones/2024.03_diff_loop_analysis_extraannotations/diffpeakachu/Mutant-WT.merged.loops

        # The output of these two commands are stored here:
        /data/humangen_external/test_area/varun/2024.08.19_hiccups_to_get_loop_scores
        '''
}
*/
