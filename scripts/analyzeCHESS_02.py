#!/bin/python

############################
### Code created by Nick Noel Machnik and published at
### https://github.com/vaquerizaslab/chess/blob/master/examples/dlbcl/example_analysis.ipynb
### Modified by K. Schultz, August 2021
############################

import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import matplotlib
from matplotlib import pyplot as plt
import fanc
import fanc.plotting
from scipy import ndimage as ndi
import matplotlib.patches as patches
from scipy.ndimage import zoom
from mpl_toolkits.axes_grid1 import make_axes_locatable

prop_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
SMALL_SIZE = 13
MEDIUM_SIZE = 15
BIGGER_SIZE = 20

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


def clipped_zoom(img, zoom_factor, **kwargs):
    h, w = img.shape[:2]
    zoom_tuple = (zoom_factor,) * 2 + (1,) * (img.ndim - 2)
    if zoom_factor < 1:
        zh = int(np.round(h * zoom_factor))
        zw = int(np.round(w * zoom_factor))
        top = (h - zh) // 2
        left = (w - zw) // 2
        out = np.zeros_like(img)
        out[top:top+zh, left:left+zw] = zoom(img, zoom_tuple, **kwargs)
    elif zoom_factor > 1:
        zh = int(np.round(h / zoom_factor))
        zw = int(np.round(w / zoom_factor))
        top = (h - zh) // 2
        left = (w - zw) // 2
        out = zoom(img[top:top+zh, left:left+zw], zoom_tuple, **kwargs)
        trim_top = ((out.shape[0] - h) // 2)
        trim_left = ((out.shape[1] - w) // 2)
        out = out[trim_top:trim_top+h, trim_left:trim_left+w]
    else:
        out = img
    return out

def highlight_features(dataframe, region, color, a, axes):
    try:
        features = dataframe.loc[region].values.tolist()
        if type(features[0]) == int:
            _, x_min, x_max, y_min, y_max = features
            rect = patches.Rectangle((x_min,y_min),x_max-x_min,y_max-y_min,linewidth=1.2,
                                     edgecolor=color, facecolor='none')
            axes[a].add_patch(rect)
        else:
            for f in features:
                _, x_min, x_max, y_min, y_max = f
                rect = patches.Rectangle((x_min,y_min),x_max-x_min,y_max-y_min,linewidth=1.2,
                                         edgecolor=color, facecolor='none')
                axes[a].add_patch(rect)
                
    except KeyError:
        next

winsize = "changeme.kb"
wdir = "./"

chess_results_file = "changeme.tsv"

region_pairs = "changeme.bed"

similarities = pd.read_csv(wdir + chess_results_file, sep='\t', index_col=0)
regions = pd.read_csv(wdir + region_pairs, sep='\t', header=None)

sim_field = "z_ssim"
sn_thr = 0.5
zsim_thr = -1
sub_sim = similarities[(similarities["SN"]>= sn_thr) & (similarities[sim_field]<= zsim_thr)]

all_X = regions.loc[similarities.index, 1:2].mean(axis=1).values / 10 ** 6
X = regions.loc[sub_sim.index, 1:2].mean(axis=1).values / 10 ** 6
S = sub_sim[sim_field]
SN = sub_sim["SN"]
plt.figure(figsize=(10, 3))
plt.plot(all_X, similarities[sim_field], ":", alpha=0.4)
plt.hlines(zsim_thr, 0, max(all_X), linestyle=":", color="red")
plt.scatter(all_X, similarities[sim_field], facecolors='none', edgecolors='grey', alpha=0.1, s=20)
plt.scatter(X, S, c=SN, marker='.')
plt.ylabel(sim_field.replace("_", "-"))
plt.xlabel("window midpoint [Mb]")
c = plt.colorbar()
c.set_label("signal/noise")
plt.tight_layout()
plt.savefig(wdir + "/plots/{}_results_track.png".format(winsize), dpi=250)

### create links to workdir beforehand
patient_hic = fanc.load(wdir + "changeme.hic")
control_hic = fanc.load(wdir + "changeme.hic")

sub_sim['id'] = sub_sim.index 
sub_sim = sub_sim.dropna()
sub_sim = sub_sim.sort_values("z_ssim")

maxPlots = 20
candidates = (sub_sim.loc[:,'id']).head(maxPlots)
f = open(wdir + "ranking.txt", "a")
f.write("#rank\tregion_id\tposition\n")
for idx, candidate in enumerate(candidates): 

	region_id = candidate
	chromosome, window_start, window_end = regions.loc[region_id][0:3]

	region_string = "{}:{}-{}".format(chromosome, window_start, window_end)
	f.write(str(idx +1) + "\t" + str(region_id) + "\t" + region_string + "\n")

	patient_region_sub = patient_hic[region_string, region_string].data
	control_region_sub = control_hic[region_string, region_string].data

	min_v = min(
		[
			np.min(np.extract(patient_region_sub>0 , patient_region_sub)),
			np.min(np.extract(control_region_sub>0 , control_region_sub))
		]
	)

	patient_region_sub += min_v
	control_region_sub += min_v

	l2fcm = np.log2(patient_region_sub / control_region_sub)

	fig, axes = plt.subplots(1, 3, figsize=(9, 3))

	axes[0].set_title('patient')
	axes[1].set_title('control')
	axes[2].set_title('$log_2$(patient / control)')
	
	fig.suptitle(region_string, fontsize=10)
	
	m1 = axes[0].imshow(patient_region_sub, norm=matplotlib.colors.LogNorm(), cmap='germany')
	m2 = axes[1].imshow(control_region_sub, norm=matplotlib.colors.LogNorm(), cmap='germany')
	m3 = axes[2].imshow(l2fcm, cmap='seismic', vmax=5, vmin=-5)
	for m, ax in zip([m1, m2, m3], axes):
		ax.axis('off')
		divider = make_axes_locatable(ax)
		cax = divider.append_axes('bottom', size='5%', pad=0.05)
		fig.colorbar(m, cax=cax, orientation='horizontal')

	plt.tight_layout()
	plt.savefig(wdir + "/plots/{}_region.png".format(candidate), dpi=250)

	##################
	# Feature extraction
	##################

	### obtaining regions of interest
	regions2compare = regions.loc[sub_sim.index]
	regions2compare.to_csv('filtered_regions_{}.tsv'.format(winsize), '\t', index=False, header=False)

	### load gained and lost features
	try:
		gained = pd.read_csv(wdir + 'features/gained_features.tsv'.format(candidate), delimiter=',', usecols=[0, 1, 2, 3,4, 5], header=None, index_col=[0])
	except:
		print("No gained features.")
		gained = pd.DataFrame([])
	try:
		lost = pd.read_csv(wdir + 'features/lost_features.tsv'.format(candidate), delimiter=',', usecols=[0, 1, 2, 3, 4, 5], header=None, index_col=[0])
	except:
		print("No lost features.")
		lost = pd.DataFrame([])

	patient_region_sub = patient_hic[region_string, region_string].data
	control_region_sub = control_hic[region_string, region_string].data

	zml2 = clipped_zoom(l2fcm, 0.7)
	rot_l2 = ndi.rotate(zml2, 45, reshape=False)

	fig, axes = plt.subplots(1, 3, figsize=(16, 8))

	### clipped zoom and rotate patient and control and keep only half-matrix
	zm1 = clipped_zoom(control_region_sub, 0.7)
	rot_control = ndi.rotate(zm1, 45, reshape=False)

	zm2 = clipped_zoom(patient_region_sub, 0.7)
	rot_patient = ndi.rotate(zm2, 45, reshape=False)

	middle = int(np.shape(rot_control)[1]/ 2.)

	fig.suptitle(region_string, fontsize=20)

	m1 = axes[0].imshow(rot_patient[:middle, :], vmin=0, vmax=0.03, cmap='germany')
	m2 = axes[1].imshow(rot_control[:middle,:], vmin=0, vmax=0.03, cmap='germany')

	### per region check if identified features, to highlight
	if not gained.empty:
		highlight_features(gained, region_id, 'crimson', 0, axes)
	if not lost.empty:
		highlight_features(lost, region_id, 'royalblue', 1, axes)

	m3 = axes[2].imshow(rot_l2[:middle,:], cmap='seismic', vmax=5, vmin=-5)

	for m, ax in zip([m1, m2, m3], axes):
		ax.axis('off')
		divider = make_axes_locatable(ax)
		cax = divider.append_axes('bottom', size='5%', pad=0.05)
		fig.colorbar(m, cax=cax, orientation='horizontal')

	plt.tight_layout()
	plt.savefig(wdir + "/plots/{}_region_with_features.png".format(candidate), dpi=250)
	plt.close('all')

f.close()