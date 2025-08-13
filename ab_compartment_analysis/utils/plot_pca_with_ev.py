#!/usr/bin/env python
import os
import sys
from itertools import repeat

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# Get all the plotting parameters
with open("input_files_for_pca.yaml", "r") as file:
    params = yaml.safe_load(file)
names = [sample["name"] for sample in params["samples"]]
labels = [
    sample["condition1"] + "_" + sample["condition2"] for sample in params["samples"]
]
files = [sample["ev_file"] for sample in params["samples"]]

# Load all the data
dfs = list()
for sample in params["samples"]:
    df = pd.read_csv(
        sample["ev_file"],
        sep="\t",
        index_col=False,
        names=[
            "chrom",
            "start",
            "end",
            "x",
            f"{sample["name"]}",
            "y",
        ],
    )
    df = df.iloc[:, [0, 1, 2, 4]].set_index(["chrom", "start", "end"])
    dfs.append(df)

df = dfs.pop(0)
df = df.join(dfs).transpose()
# df = StandardScaler().fit_transform(df)

# remove all colums that only has zeros
# df = df.loc[:, (df != 0).any(axis=0)]

# Do PCA
pca = PCA(n_components=12)
dimred = pca.fit_transform(df)

lb = list(set(labels))
lb.sort()
colors = ["red" if "DNMT3A" in label else "blue" for label in lb]
markers = ["o" if "Mut" in label else "s" for label in lb]
sizes = [10 if "Replicate" in label else 30 for label in lb]

fig, ax = plt.subplots()
for label in lb:
    plot_id = [x for x, value in enumerate(lb) if value == label]
    sample_id = [x for x, value in enumerate(labels) if value == label]
    ax.scatter(
        dimred[sample_id, 0],
        dimred[sample_id, 1],
        c=colors[plot_id[0]],
        label=label,
        marker=markers[plot_id[0]],
        s=sizes[plot_id[0]],
    )

ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
ax.legend()
# plt.show()
plt.savefig("PCA_DNMT3A_IDH1_AB_eigen_vector_python.png")
