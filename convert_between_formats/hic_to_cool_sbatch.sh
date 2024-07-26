#! /bin/bash

### Submit this Script with: sbatch <hic_to_cool_sbatch.sh> ###
#SBATCH --partition=shortterm
#SBATCH --nodes=1
#SBATCH -c 1
#SBATCH -J "hic_to_cool"
#SBATCH --mem=2GB
#set slurm file output nomenclature
#SBATCH --output "slurm-%x-%j.out"

PATH=$WORK/.omics/anaconda3/bin:$PATH #add the anaconda installation path to the bash path
source $WORK/.omics/anaconda3/etc/profile.d/conda.sh # some reason conda commands are not added by default

conda activate nextflow
#
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
