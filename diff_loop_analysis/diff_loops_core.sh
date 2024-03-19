#! /bin/bash

### Submit this Script with: sbatch <peakachu_core.sh> ###
#SBATCH --partition=shortterm
#SBATCH --nodes=1
#SBATCH -c 1
#SBATCH -J "diff_loops"
#SBATCH --mem=10GB
#set slurm file output nomenclature
#SBATCH --output "slurm-%x-%j.out"

# Load your necessary modules:
module load nextflow/v22.04.1

# copy the scrits to $WORK directory - otherwise there is a jav aerror
mkdir -p "$WORK"/diff_loops_nextflow_launchdir
rsync -r --update * "$WORK"/diff_loops_nextflow_launchdir/
cd "$WORK"/diff_loops_nextflow_launchdir
chmod +x bin/*
rm slurm*.out

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

# Submit the Nextflow Script:
nextflow run diff_loops.nf -params-file diff_loops_params.yaml -resume
