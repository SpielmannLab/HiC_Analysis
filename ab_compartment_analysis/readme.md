# Perform AB compartment analysis

This pipeline contains two steps.

1. Use FANC to call the AB compartments and do some plotting.
   It is important to keep in mind that the eigen vectors might be of the wrong orientation and the assignment of the AB compartments may be switched. Moreover, in some cases, the 1st eigen vector may not have much to do with the compartments, but the chromosome arms instead. In this case, assignments might be completely off.

   Run this with: `sbatch ab_compartments_sbatch.sh`

2. Fix the orientation of the eigen vector using functional data.
   Use some data (accessibility, histone or gene density) relevant to the cell type to fix the orientation of the eigen vector in one reference sample. Then the eigen vectors and the domain assignments of all samples to this aligned reference eigen vectors. This also does plotting of enrichment of contacts between A and B compartments. And also provides you with compartment plots.

Run the pipeline. There seems to be not much happening with or without nromalization when calling AB compartments.

TODO
Important notes:
1. I am not removing the blacklisted regions of the genome. People do this to remove spurious signals. You may way to do this while doing PCA analysis using the AB eigenvectors.


## Utility scripts

It is often interesting to make a PCA plot that shows how the AB compartments are different between samples. The script [plot_pca_with_ev.py](utils/plot_pca_with_ev.py) uses the (aligned/corrected) eigen vectors across from all samples. 
