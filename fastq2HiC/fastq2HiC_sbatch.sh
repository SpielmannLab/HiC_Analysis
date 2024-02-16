#!/usr/bin/bash

#SBATCH -p shortterm
#SBATCH -t 3-0:0:0
#SBATCH -c 3
#SBATCH --mem=10GB

# Load your necessary modules:
module load nextflow/v22.04.1

# copy the scrits to $WORK directory - otherwise there is a jav aerror
mkdir -p $WORK/hic_nextflow_launchdir
rsync -r --update * $WORK/hic_nextflow_launchdir/
cd $WORK/hic_nextflow_launchdir
rm slurm*.out

module load nextflow/v22.04.1

# if -resume option not present, then clean the nextflow
if [[ $1 = -resume ]]; then
    echo "Resuming the previous nextflow run"
else
    echo "Cleaning up all the previous nextflow runs"
    echo "If you wanted to resume the previous run, use the '-resume' option"
    while [[ $(nextflow log | wc -l) -gt 1 ]]; do
        nextflow clean -f
    done
fi

nextflow run fastq2HiC.nf -params-file fastq2HiC_params.yaml -with-trace -resume
