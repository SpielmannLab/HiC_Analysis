#! /bin/bash

### Submit this Script with: sbatch <hic_to_cool_sbatch.sh> ###
#SBATCH --partition=shortterm
#SBATCH --nodes=1
#SBATCH -c 1
#SBATCH -J "hic_to_cool"
#SBATCH --mem=2GB
#set slurm file output nomenclature
#SBATCH --output "slurm-%x-%j.out"

# Load your necessary modules:
module load nextflow/v22.04.1
#
# copy the scrits to $WORK directory - otherwise there is a jav aerror
mkdir -p $WORK/hic_to_cool_nextflow_launchdir
rsync -rc * $WORK/hic_to_cool_nextflow_launchdir/
cd $WORK/hic_to_cool_nextflow_launchdir
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
nextflow run hic_to_cool.nf -params-file hic_to_cool_params.yaml -resume
