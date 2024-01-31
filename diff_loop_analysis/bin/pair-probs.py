#!/usr/bin/env python

# The script stores the probabilites of the loop probabilities two comparative groups in a pairwise manner per co-ordinate.
# It automatically looks

# Usage: pair-probs.py group1_merged_loops group2_merged_loops outputfile
# ----- Description of input
# group{1,2}_merged_loops   These are the files containing all confident loops by merging all the repeats etc.
# outputfile                Name of the output file

import os
import sys


# This function gets the max probability at each location from individual unfiltered files
def parse_probs(file_list):
    D = {}
    for file in file_list:
        with open(file, "r") as source:
            for line in source:
                p = line.rstrip().split()
                chr1 = p[0]
                chr2 = p[3]
                loci1 = int(p[1])
                loci2 = int(p[4])
                prob = float(p[6])
                if (chr1, chr2, loci1, loci2) in D:
                    D[(chr1, chr2, loci1, loci2)] = max(
                        D[(chr1, chr2, loci1, loci2)], prob
                    )
                else:
                    D[(chr1, chr2, loci1, loci2)] = prob
    return D


# This function returns just a set of all distinct loops present in the data, ignoring the probabilities
def parse_loops(fil1, fil2):
    loop_set = set()
    with open(fil1, "r") as source:
        for line in source:
            p = line.rstrip().split()
            # This just gets all the columns except the loop probs(cols 6 & 7)
            loop_set.add(tuple(p[:6]))
    with open(fil2, "r") as source:
        for line in source:
            p = line.rstrip().split()
            # This just gets all the columns except the loop probs(cols 6 & 7)
            loop_set.add(tuple(p[:6]))
    return loop_set


# ---- parse arguments  -----
file_group1_merged = sys.argv[1]
file_group2_merged = sys.argv[2]
filepattern_group1_replicates = sys.argv[3]
filepattern_group2_replicates = sys.argv[4]
outfile = sys.argv[5]
# Gather the coordinates of all high confidence loops in both groups that passed the threshold. Called using "peakachu pool"
loop_set = parse_loops(file_group1_merged, file_group2_merged)

# Now gather the max probabilites for the above loops from amongst the unfiltered individual repeats
group1_replicate_files = [
    file for file in os.listdir() if file.startswith(filepattern_group1_replicates)
]
group2_replicate_files = [
    file for file in os.listdir() if file.startswith(filepattern_group2_replicates)
]
group1_probs = parse_probs(group1_replicate_files)
group2_probs = parse_probs(group2_replicate_files)

# Collect the max probabiliites from the two groups into a single list
pool = []
for loop in loop_set:
    key = (loop[0], loop[3], int(loop[1]), int(loop[4]))
    p1 = group1_probs.get(key, 0)
    p2 = group2_probs.get(key, 0)
    if p1 or p2:
        pool.append(loop + ("{0:.4g}".format(p1), "{0:.4g}".format(p2)))

# Write to the output file
with open(outfile, "w") as out:
    for line in pool:
        out.write("\t".join(list(line)) + "\n")
