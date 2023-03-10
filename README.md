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
We're offering two pipeline to align HiC Fastq files.

HiC-Pro creates .validPairs files and statistical information that is required for some other tools. 
The validPairs can get converted to .hic, .cool or .matrix format.

The Juicer scripts create .hic files and also offer scaffolding. Scaffolding is useful when references are limited or when dealing with many large structural variants (e.g. chromothripsis).

## HiC-Pro
For Alignment the pipeline starts HiC-Pro in parallel mode on a SLURM cluster.
Please install HiC-Pro by following the [instructions](https://github.com/nservant/HiC-Pro). The path to the installation directory in required in the configuration file.
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

## Juicer
Clone or download Aiden Lab [Juicer repo](https://github.com/aidenlab/juicer).

Adjust juicer.sh depending on Cluster System. 
For example when working with the omics cluster in Luebeck you may add the line *isOMICS=1 .* to /juicer/SLURM/juicer.sh and extent the following if statement with

```bash
elif [ $isOMICS -eq 1 ]
then    
    load_bwa="module load bwa/v0.7.17"
    queue="shortterm"
    queue_time="0-12:00:00"
    long_queue="longterm"
    long_queue_time="3-00:00:00"
    load_gpu="module load nvidia-cuda/11.1-native"
```

For Scaffolding 
Clone or download Aiden Lab [3d-dna repo](https://github.com/aidenlab/3d-dna).
The workflow will use a singularity image of [w2rap-contigger](https://github.com/bioinfologics/w2rap-contigger), but not with the snakemake container setting.
If you prefer to use a local installation, remove the "singularity exec {params.img}" in Juicer.smk.

```bash
### Alignment to reference, i.a. creates .bam and .hic files i.a. 
### The script spawns multiple jobs. Once they've finished the output files are in the output folder
snakemake -s Juicer.smk juicerPipeline --cores $n --config SCRATCH=$SCRATCH 
```


```bash
### Scaffold Alignment, creates .bam 
### The script spawns multiple jobs. Wait for them to finish before starting rule asmPipeline. 
snakemake -s Juicer.smk juicerAssembly --cores $n --config SCRATCH=$SCRATCH   
### asmPipeline creates a .hic and .assembly file. 
snakemake -s Juicer.smk asmPipeline --cores $n --config SCRATCH=$SCRATCH 
```

## Quality Control
The tool hicrep works with samples in the .cool format. The hicpro.smk contains rules for converting .hic to .cool, which will be included in the workflow when calling the rule *hicrep*.

```bash
### quality control with hicrep creates QC folder with stratum adjusted correlation coefficient scores
snakemake -s hicpro.smk hicrep --cores $n --use-conda --config SCRATCH=$SCRATCH  
```

## Denoise
This pipeline runs a denoise functionality provided by the tools HiCorr and DeepLoop. 
The output data can **not** get visualized like a .hic or .cool file, but one can generate plots of selected regions (here target gene regions).
In this pipeline the tools HiCorr and DeepLoop mix and mingle and some scripts might not get the correct paths. At the time of this writing HiCorr was published three months ago, so a lot of things might change soon.
For now setting up this pipeline might require some troubleshooting with where to put the downloaded files.
Please visit [HiCorr](https://github.com/JinLabBioinfo/HiCorr) and [DeepLoop](https://github.com/JinLabBioinfo/DeepLoop) for more information.

```bash
### quality control with hicrep creates QC folder with stratum adjusted correlation coefficient scores
snakemake -s denoise.smk  plot_all --cores $n --use-conda --config SCRATCH=$SCRATCH  
```

# Analysis

![Recommended Process](https://data.cyverse.org/dav-anon/iplant/home/kschultz/flowchart_HiC.png)

## SV Analysis
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

## FANC loops and insulation 
```bash
### detect loops and create plot of target gene regions
snakemake -s fanc.smk all_loops --cores $n --use-conda --config SCRATCH=$SCRATCH 
### insulation score for target gene regions
snakemake -s fanc.smk insulation_plots --cores $n --use-conda --config SCRATCH=$SCRATCH 
### aggregate analysis 
snakemake -s fanc.smk plot_APA --cores $n --use-conda --config SCRATCH=$SCRATCH 
snakemake -s fanc.smk matrix_APA --cores $n --use-conda --config SCRATCH=$SCRATCH 
```

## A/B Compartments
Compartment Analysis with FANC (works with fanc.yml). Start an interactive job with xx cores and xx GB memory. It uses the conda environment *env/fanc.yml*. But for some reason it throws errors. Therefore, create the conda environment first and then run it. The script has been adapted accordingly.

```bash
### for every sample create chromosome overview figures
snakemake -s Compartment_Analysis.smk genome_overview --cores $n --use-conda --config SCRATCH=$SCRATCH  

### create a plot for the region around target genes
snakemake -s Compartment_Analysis.smk target_genes --cores $n --use-conda --config SCRATCH=$SCRATCH

### compare compartment shifts to a control group
snakemake -s Compartment_Analysis.smk compareToCtrls --cores $n --use-conda --config SCRATCH=$SCRATCH
```

## Loop Calling 
### Peakachu
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

### Other Loop Caller
Rules hiccups and HiC LDNet start new batch jobs with gcpu request.
The rule itself only creates the output directory, not the files. Please wait for the jobs to finish.
```bash
### run hiccups and HiC_LDNet and wait for the job to finish
snakemake -s loopcalling.smk hiccups_all --cores $n --use-conda --config SCRATCH=$SCRATCH
snakemake -s loopcalling.smk HiC_LDNet_all --cores $n --use-conda --config SCRATCH=$SCRATCH  
###
snakemake -s loopcalling.smk arrowhead_all --cores $n --use-conda --config SCRATCH=$SCRATCH    
```


## CHESS HiC Analysis
```bash
snakemake -s CHESS.smk all --cores $n --use-conda --config SCRATCH=$SCRATCH  
```


# HiChIP
FitHiChIP requires aligned validPairs files (from hicpro "hic_results" ).

```bash
### run fithichip
snakemake -s HiChIP.smk all --cores $n --use-conda --config SCRATCH=$SCRATCH  
```

