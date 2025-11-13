# Comparison of Hi-C Experiments using Structural Similarity (CHESS) analysis

## Basics on CHESS analysis
Enables shotgun comparison between two HiC maps to detects differences in various chromatin structures (incl. TAD, Loops and Stripes). 
Check the documentation [at](https://chess-hic.readthedocs.io/en/latest/example_analysis.html#finding-differences-between-samples-of-different-conditions). Check out a more recent [paper](https://doi.org/10.1101/2021.10.18.464422) on choosing parameters.

## Basic scripts

    # Create bedpe file to compare regions between two HiC files
    chess pairs mm10 3000000 100000 ./mm10_chr2_3mb_win_100kb_step.bed --chromosome chr2

    # Run the search with CHESS sim command. Make sure to note the resolution and the normalization using the @
    chess sim /data/humangen_external/HiC/khandanpour_hic/data_n_results/2024.01_fastq2HiC/GFI136N_WT_KI/juicer_1.19.02/GFI136N_WT_KI.allValidPairs.hic@100kb@VC /data/humangen_external/HiC/khandanpour_hic/data_n_results/2024.01_fastq2HiC/GFI136S_KI_KI/juicer_1.19.02/GFI136S_KI_KI.allValidPairs.hic@100kb@VC mm10_chr2_3mb_win_100kb_step.bed output.tsv

    # Plotting in python
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
    
    # setup matplotlib
    %matplotlib inline
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

    winsize = "3mb"
    wdir = "./"
    chess_results_file = "output.tsv"

    region_pairs = "mm10_chr2_{}_win_100kb_step.bed".format(winsize)

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
    plt.savefig("plots/chr2_{}_results_track.png".format(winsize), dpi=250)

## Sort identified locations by the z_ssim
similarities[(similarities["SN"]>= sn_thr) & (similarities["z_ssim"] < -1)].sort_values("ID")

## Plotting of HiC data at these identified locations
patient_hic = fanc.load("/data/humangen_external/HiC/khandanpour_hic/data_n_results/2024.01_fastq2HiC/GFI136N_WT_KI/juicer_1.19.02/GFI136N_WT_KI.allValidPairs.hic@100kb@VC")
control_hic = fanc.load("/data/humangen_external/HiC/khandanpour_hic/data_n_results/2024.01_fastq2HiC/GFI136S_KI_KI/juicer_1.19.02/GFI136S_KI_KI.allValidPairs.hic@100kb@VC")

region_id = 860
window_start, window_end = regions.loc[region_id][1:3]
region_string = "2:{}-{}".format(window_start, window_end)
print(region_string)

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

m1 = axes[0].imshow(patient_region_sub, norm=matplotlib.colors.LogNorm(), cmap='germany')
m2 = axes[1].imshow(control_region_sub, norm=matplotlib.colors.LogNorm(), cmap='germany')
m3 = axes[2].imshow(l2fcm, cmap='seismic', vmax=5, vmin=-5)
for m, ax in zip([m1, m2, m3], axes):
    ax.axis('off')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('bottom', size='5%', pad=0.05)
    fig.colorbar(m, cax=cax, orientation='horizontal')

plt.tight_layout()
plt.savefig("plots/chr2_example_region.png", dpi=250)
