# Differential Peakachu analysis log

## Peakachu installation

Created peakachu env using installation instructions at their [github page](https://github.com/tariks/peakachu#installation):

## Analysis steps

1. After running the [loop analysis](loop_analysis/call_loops_core.sh) script for multiple conditions, we are now ready to called differential loops between the conditionsa.

2. Fill out the [parameters file](diff_loop_analysis/diff_loops_params.yaml). Make sure to supply both the filtered and unfiltered loops.

3. Run diff_loop_analysis to call peaks.

   The script is a slightly modified version of the [diffPeakachu](https://github.com/tariks/peakachu/tree/master/diffPeakachu). To run it, copy the script to anywhere in the omics cluster and submit:

   ```bash
   sbatch diff_loops_core.sh
   ```

4. The result will be of the following structure

   ```bash
   sbatch diff_loops_core.sh
   ├── diffpeakachu # Contains the differential loops unique to the two groups. Also has a PNG file showing the Gausian Mixture Model fit
   │   ├── Mutant-WT.Fold-GMM.png
   │   ├── Mutant-WT.Mutant.unique.loops
   │   └── Mutant-WT.WT.unique.loops
   ├── genes_in_loops # Gene annotations for the files in the other two folders
   │   ├── overlapping_genes_in_Mutant-WT.Mutant.unique.loops
   │   ├── overlapping_genes_in_Mutant-WT.WT.unique.loops
   │   ├── overlapping_genes_in_Mutant_merged.loops
   │   └── overlapping_genes_in_WT_merged.loops
   └── union_of_loops # Contains the union of all filtered loops from the samples belonging to the two groups
       ├── Mutant_merged.loops
       ├── Mutant_merged_forIGV.bedpe
       ├── WT_merged.loops
       └── WT_merged_forIGV.bedpe
   ```

   In the Group1-Group2.Group1(2).unique.loops file, there will be 8 column:

   - Column 1: chr1 - Chromosome of left side of the loop
   - Column 2: x1 - Left of the bin for the left side of the loop
   - Column 3: x2 - Right of the bin for the left side of the loop
   - Column 4: chr2 - Chromosome of right side of the loop
   - Column 2: y1 - Left of the bin for the right side of the loop
   - Column 3: y2 - Right of the bin for the right side of the loop
   - Column 7: score_this_group - The peakachu score of this loop for the Group1(2)
   - Column 8: score_other_group - The peakachu score of this loop for the Group2(1)

5. The next step is to run Aggregate Peak Analysis

   TODO
