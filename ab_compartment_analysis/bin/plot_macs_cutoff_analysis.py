#!/usr/bin/env python3
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("chippeaks_cutoff_analysis.txt", sep="\t")

fig, axes = plt.subplots(figsize=(3, 5), ncols=1,
                         nrows=3, layout="constrained")
axes[0].plot(df.score, df.npeaks)
axes[1].plot(df.score, df.lpeaks)
axes[2].plot(df.score, df.avelpeak)
titles = ["No. of peaks", "Total length of all peaks",
          "Average length of peaks"]

for idx in [0, 1, 2]:
    axes[idx].set_yscale("log")
    axes[idx].set_title(titles[idx])

axes[2].set_xlabel("Peak threshold")
plt.savefig("chippeaks_cutoff_analysis.png")
