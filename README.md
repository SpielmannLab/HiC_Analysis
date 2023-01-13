# Pipeline running HiC-Pro on SLURM

This pipeline uses Snakemake v6.  
Download example files

```bash
wget https://data.cyverse.org/dav-anon/iplant/home/kschultz/example.zip
```

# Alignment
For Alignment the pipeline starts HiC-Pro in parallel mode on a SLURM cluster.

```bash
### when run in separate steps use permanent work directory
snakemake -s hicpro.smk hicpro_parallel_step1 --cores $n --use-conda
snakemake -s hicpro.smk hicpro_parallel_step2 --cores $n --use-conda
### 
```

```bash
### start both steps
snakemake -s hicpro.smk all --cores $n --use-conda --config SCRATCH=$SCRATCH  
```

```bash
### pool all samples into one (adjust resources depending on the file size you're creating)
snakemake -s hicpro.smk pool --cores $n --use-conda --config SCRATCH=$SCRATCH  
```

# SV Analysis
For Structural Variant Calling we use HiC-Breakfinder and for Analysis NeoloopFinder.  
For HiC-Breakfinder please download associated files:
https://salkinstitute.box.com/s/m8oyv2ypf8o3kcdsybzcmrpg032xnrgx

```bash
### convert hic to cool format
snakemake -s neoloopFinder.smk all_cool --cores $n --use-conda --config SCRATCH=$SCRATCH  
```

```bash
### create neoloop tads and loops files
snakemake -s neoloopFinder.smk all --cores $n --use-conda --config SCRATCH=$SCRATCH  
```

```bash
### Plot breakpoints with
snakemake -s neoloopFinder.smk plot --cores $n --use-conda --config SCRATCH=$SCRATCH  
```

# FANC loops and insulation 
```bash
### detect loops and create plot of target gene regions
snakemake -s fanc.smk all_loops --cores $n --use-conda --config SCRATCH=$SCRATCH 
### insulation score for target gene regions
snakemake -s fanc.smk insulation_plots --cores $n --use-conda --config SCRATCH=$SCRATCH 
### aggregate analysis 
snakemake -s fanc.smk plot_APA --cores $n --use-conda --config SCRATCH=$SCRATCH 
snakemake -s fanc.smk matrix_APA --cores $n --use-conda --config SCRATCH=$SCRATCH 
```


# A/B Compartments
Compartment Analysis with FANC (works with fanc.yml)
```bash
### for every sample create chromosome overview figures
snakemake -s Compartment_Analysis.smk genome_overview --cores $n --use-conda --config SCRATCH=$SCRATCH  

### create a plot for the region around target genes
snakemake -s Compartment_Analysis.smk target_genes --cores $n --use-conda --config SCRATCH=$SCRATCH

### compare compartment shifts to a control group
snakemake -s Compartment_Analysis.smk compareToCtrls --cores $n --use-conda --config SCRATCH=$SCRATCH
```

# CHESS HiC Analysis
```bash
snakemake -s CHESS.smk all --cores $n --use-conda --config SCRATCH=$SCRATCH  
```


# HiChIP
FitHiChIP requires aligned validPairs files (from hicpro "hic_results" ).

```bash
### run fithichip
snakemake -s HiChIP.smk all --cores $n --use-conda --config SCRATCH=$SCRATCH  
```

