# Pipeline running HiC-Pro on SLURM

This pipeline uses Snakemake v6. 


# For Alignment the pipeline starts HiC-Pro in parallel mode on a SLURM cluster.

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
### pool all samples into one (adjust ressources depending on the file size you're creating)
snakemake -s hicpro.smk pool --cores $n --use-conda --config SCRATCH=$SCRATCH  
```

