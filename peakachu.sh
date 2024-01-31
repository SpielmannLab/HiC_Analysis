#! /bin/bash

#SBATCH --partition=shortterm
#SBATCH -c 3
#SBATCH --mem=200GB
#SBATCH --mail-type=END
#SBATCH --job-name=peakachu
##SBATCH --tmp=40GB


PATH=$WORK/.omics/anaconda3/bin:$PATH #add the anaconda installation path to the bash path
source $WORK/.omics/anaconda3/etc/profile.d/conda.sh # some reason conda commands are not added by default

cat peakachu.yml

# Run peakachu version 2
snakemake -s peakachu.smk peakachu_v2 --cores 3 --use-conda --config SCRATCH=$SCRATCH --conda-prefix "$WORK/hic_conda_envs/" --conda-frontend "mamba"