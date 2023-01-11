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
For Structural Variant Calling we use HiC-Breakfinder and for Analysis NeoloopFinder. HiC-Breakfinder is currently limited to using hg19 resources. 
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

# HiChIP
FitHiChIP requires aligned validPairs files (from hicpro "hic_results" ).

```bash
### run fithichip
snakemake -s HiChIP.smk all --cores $n --use-conda --config SCRATCH=$SCRATCH  
```

