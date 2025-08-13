#! /bin/bash

### Submit this Script with: sbatch <peakachu_core.sh> ###
#SBATCH --partition=shortterm
#SBATCH --nodes=1
#SBATCH -c 1
#SBATCH -J "ab_correction"
#SBATCH --mem=10GB
#set slurm file output nomenclature
#SBATCH --output "slurm-%x-%j.out"

module load nextflow
module load singularity

# copy the scrits to $WORK directory - otherwise there is a jav aerror
mkdir -p $WORK/ab_correction_launchdir
rsync --exclude "slurm*" -rc --update * $WORK/ab_correction_launchdir/
cd $WORK/ab_correction_launchdir

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
nextflow -C ab_correction_nextflow.config run ab_correction_main.nf \
  -params-file ab_correction_params.yaml \
  -resume -ansi-log false \
  -profile omics \
