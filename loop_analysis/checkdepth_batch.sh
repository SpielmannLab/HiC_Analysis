#!/usr/bin/bash

#SBATCH -c 1
#SBATCH --mem=2GB
#SBATCH -p debug 
#SBATCH -t 0:10:0 

# initiate conda environment
PATH=$WORK/.omics/anaconda3/bin:$PATH #add the anaconda installation path to the bash path
source $WORK/.omics/anaconda3/etc/profile.d/conda.sh # some reason conda commands are not added by default
conda activate peakachu

hic_file=$1
echo "Checking depth of $hic_file using peakachu depth"

# tee writes to stdout and to a file
peakachu depth -p $hic_file | tee ${hic_file/.hic/_depth.txt}
