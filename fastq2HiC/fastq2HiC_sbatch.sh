#!/usr/bin/bash

#SBATCH -p shortterm
#SBATCH -t 3-0:0:0
#SBATCH -c 3
#SBATCH --mem=10GB

# Load your necessary modules:
module load nextflow/v22.04.1

# copy the scrits to $WORK directory - otherwise there is a jav aerror
mkdir -p $WORK/hic_nextflow_launchdir
find . -not -name "slurm*" -exec cp -r {} $WORK/hic_nextflow_launchdir/ \;
cd $WORK/hic_nextflow_launchdir

nextflow run fastq2HiC.nf -params-file fastq2HiC_params.yaml -with-trace
