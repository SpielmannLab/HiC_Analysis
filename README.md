# Pipeline running HiC-Pro on SLURM

## General
This pipeline was created using Snakemake version 5, but will also work with version 6.  
Download example (incomplete) files

```bash
wget https://data.cyverse.org/dav-anon/iplant/home/kschultz/example.zip
```

There are multiple separate workflow steps. 
While most rules in a workflow use conda environments that Snakemake will automatically create, some tools, such as HiC-Pro, are not available with conda and require a manual installation.
Please follow the instructions provided by the authors. 
Some rules also need access to some supplementary files, that are provided by the original distributor.
\
Most workflows consist of a Snakefile \<something\>.smk and a corresponding configuration file \<something\>.yml
You can, however, use different config files or overwrite a variable set in in the .yml file.

```bash
### use a different config file
snakemake -s something.smk --configfile something_else.yml
### overwrite a set variable 
snakemake -s something.smk --config SAMPLES=test2
```

All workflows are using a SCRATCH directory for storing temp data. You can set any directory as scratch directory in the configuration (.yml) files.
If your scratch directory is not automatically removed after the job has finished, please remember to delete the temp files.

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

## Quality Control

```bash
### quality control with hicrep creates QC folder with stratum adjusted correlation coefficient scores
snakemake -s hicpro.smk hicrep --cores $n --use-conda --config SCRATCH=$SCRATCH  
```

## Denoise
This pipeline runs a denoise functionally provided by the tools HiCorr and DeepLoop. 
The output data can not get visualized like a .hic or .cool file, but one can generate plots of selected regions (here target gene regions).
In this pipeline the tools HiCorr and DeepLoop mix and mingle and some scripts might not get the correct paths. At the time of this writing HiCorr was published three months ago, so a lot of things might change soon.
For now setting up this pipeline might require some troubleshooting with where to put the downloaded files.
Please visit [HiCorr](https://github.com/JinLabBioinfo/HiCorr) and [DeepLoop](https://github.com/JinLabBioinfo/DeepLoop) for more information.

```bash
### quality control with hicrep creates QC folder with stratum adjusted correlation coefficient scores
snakemake -s denoise.smk  plot_all --cores $n --use-conda --config SCRATCH=$SCRATCH  
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

# Loop Calling 
## Peakachu
```bash
### peakachu version 1
snakemake -s peakachu.smk peakachu_v1 --cores $n --use-conda --config SCRATCH=$SCRATCH
### peakachu version 2
snakemake -s peakachu.smk peakachu_v2 --cores $n --use-conda --config SCRATCH=$SCRATCH    
```
Peakachu offers scripts to identify differential chromatin loops between samples.
In the configuration peakachu.yml one has to select a pool threshold. The best threshold may vary between samples. 
For best results try multiple and select the ones best working for each sample.
Create two directories peakachudiff\_input\_samples and peakachudiff\_input\_controls. Copy your selection of pooled loops calls (.txt files) into the respective directory.
All samples of configuration SAMPLES_exp1 must be present in peakachudiff\_input\_samples.

```bash
### merge control loops and compare each sample to merged set of loops
snakemake -s peakachu.smk peakachuDiff_mergedCtrl --cores $n --use-conda --config SCRATCH=$SCRATCH
### compare each sample to each control
snakemake -s peakachu.smk peakachuDiff_eachCtrl --cores $n --use-conda --config SCRATCH=$SCRATCH
```

## Other Loop Caller
Rules hiccups and HiC LDNet start new batch jobs with gcpu request.
The rule itself only creates the output directory, not the files. Please wait for the jobs to finish.
```bash
### run hiccups and HiC_LDNet and wait for the job to finish
snakemake -s loopcalling.smk hiccups_all --cores $n --use-conda --config SCRATCH=$SCRATCH
snakemake -s loopcalling.smk HiC_LDNet_all --cores $n --use-conda --config SCRATCH=$SCRATCH  
###
snakemake -s loopcalling.smk arrowhead_all --cores $n --use-conda --config SCRATCH=$SCRATCH    
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

