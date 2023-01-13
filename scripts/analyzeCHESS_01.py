#!/bin/python

############################
### Code created by Nick Noel Machnik and published at
### https://github.com/vaquerizaslab/chess/blob/master/examples/dlbcl/example_analysis.ipynb
### Modified by K. Schultz, August 2021
############################

import numpy as np
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
import fanc
import fanc.plotting
from scipy import ndimage as ndi
import matplotlib.patches as patches
from scipy.ndimage import zoom
from mpl_toolkits.axes_grid1 import make_axes_locatable

##################
# Examine chess sim results
##################

winsize = "changeme.kb"
wdir = "./"

chess_results_file = "changeme.tsv"

region_pairs = "changeme.bed"

similarities = pd.read_csv(wdir + chess_results_file, sep='\t', index_col=0)
regions = pd.read_csv(wdir + region_pairs, sep='\t', header=None)

sim_field = "z_ssim"
### default filter thesholds
sn_thr = 0.5
zsim_thr = -1
sub_sim = similarities[(similarities["SN"]>= sn_thr) & (similarities["z_ssim"]<= zsim_thr)]

### obtaining regions of interest
regions2compare = regions.reindex(index = sub_sim.index)
regions2compare.to_csv(wdir + 'filtered_regions_{}.tsv'.format(winsize), '\t', index=False, header=False)
