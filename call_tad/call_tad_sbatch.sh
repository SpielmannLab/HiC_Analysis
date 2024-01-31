#! /bin/bash

### Submit this Script with: sbatch <peakachu_core.sh> ###
#SBATCH --partition=shortterm
#SBATCH --nodes=1
#SBATCH -c 1
#SBATCH --mem=10GB
#SBATCH -J "call_tad"
#set slurm file output nomenclature
#SBATCH --output "slurm-%x-%j.out"

# Load your necessary modules:
module load nextflow/v22.04.1

# copy the scrits to $WORK directory - otherwise there is a jav aerror
mkdir -p "${WORK}/call_tad_nextflow_launchdir"
rsync -r --update * "${WORK}/call_tad_nextflow_launchdir"
cd "${WORK}/call_tad_nextflow_launchdir"
rm slurm*.out

# if -resume option not present, then clean the nextflow
if [[ $1 = -resume ]]; then
    echo "Resuming the previous nextflow run"
else
    echo "Cleaning up all the previous nextflow runs"
    echo "If you wanted to resume the previous run, use the '-resume' option when calling sbatch"
    while [[ $(nextflow log | wc -l) -gt 1 ]]; do
        nextflow clean -f
    done
fi

# Submit the Nextflow Script:
nextflow run call_tad_main.nf -params-file call_tad_params.yaml -resume
