#!/usr/bin/env python
# submit as ./plot_enrichments.py
# do this in the varunkas/fanc:latest docker/singularity
import sys

import fanc
import fanc.plotting as fancplot
import matplotlib.pyplot as plt
import numpy as np

# plot whole genome or per chromosome?
ab = fanc.load("fanc.ab")
hic = fanc.load("fanc.hic")
ev = fanc.load("ev.bed")

# Since eigenvectors of some chromsomes are not corrected, the corrected eigen vector is not of the same size as the hic matrix.
# So subset the HiC matrix to only contain those chromosomes that are in the eigen vector
# This takes 5-10 min
if len(hic.regions) != len(ev.regions):
    chroms_in_ev = ev.chromosomes()
    print(
        "Subsetting the HiC to keep only chroms in the eigen vector:"
    )
    print(", ".join(chroms_in_ev))
    # The * operator expands the list elements into a function call.
    hic = hic.subset(*chroms_in_ev)

# Now create an enrichment plot for the entire genome
print("Calculating enrichment profile for the whole genome")
profile, cutoffs = ab.enrichment_profile(
    hic, percentiles=[10, 20, 30, 40, 50, 60, 70, 80, 90, 100],
    eigenvector=[region.score for region in ev.regions]
)
fig, axes = fancplot.saddle_plot(profile, cutoffs, vmin=-2, vmax=2)
fil = "genome_wide_enrichment.png"
plt.savefig(fil)
plt.close()
compartment_strength_file = "genome_wide_enrichment.txt"
a_a = 2 ** profile[0, 0]
b_b = 2 ** profile[profile.shape[1] - 1, profile.shape[1] - 1]
a_b = 2 ** profile[0, profile.shape[1] - 1]
s = np.log((a_a * b_b) / a_b**2)
with open(compartment_strength_file, "w") as o:
    o.write("AB-strength\t{}\n".format(s))
    o.write("AA\t{}\n".format(a_a))
    o.write("BB\t{}\n".format(b_b))
    o.write("AB\t{}\n".format(a_b))

# Create an enrichment plot for each chromosome, by excluding one a time
for chrom_to_remove in hic.chromosomes():
    chroms = hic.chromosomes()
    chroms.remove(chrom_to_remove)
    print(
        f"Calculating enrichment profile for chromosome: {
            chrom_to_remove}"
    )
    profile, cutoffs = ab.enrichment_profile(
        hic,
        percentiles=[10, 20, 30, 40, 50, 60, 70, 80, 90, 100],
        exclude_chromosomes=chroms,
        eigenvector=[region.score for region in ev.regions],
        # symmetric_at=0,
    )
    fig, axes = fancplot.saddle_plot(profile, cutoffs, vmin=-2, vmax=2)
    fil = f"chr{chrom_to_remove}_enrichment.png"
    plt.savefig(fil)
    plt.close()
    compartment_strength_file = f"{chrom_to_remove}_enrichment.txt"
    a_a = 2 ** profile[0, 0]
    b_b = 2 ** profile[profile.shape[1] - 1, profile.shape[1] - 1]
    a_b = 2 ** profile[0, profile.shape[1] - 1]
    s = np.log((a_a * b_b) / a_b**2)
    with open(compartment_strength_file, "w") as o:
        o.write("AB-strength\t{}\n".format(s))
        o.write("AA\t{}\n".format(a_a))
        o.write("BB\t{}\n".format(b_b))
        o.write("AB\t{}\n".format(a_b))

# Close the files
ab.close()
hic.close()
