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
mkdir -p $WORK/call_loops_nextflow_launchdir
cp -r * $WORK/call_loops_nextflow_launchdir/
cd $WORK/call_loops_nextflow_launchdir

# Submit the Nextflow Script:
nextflow run call_loops.nf -params-file call_loops_params.yaml
