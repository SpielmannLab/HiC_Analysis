#! /bin/bash

#SBATCH --partition=shortterm
#SBATCH -c 2
#SBATCH --mem=4GB
#SBATCH --time=3-0:0:0
#SBATCH --mail-type=END
#SBATCH --mail-user=k.schultz@uni-luebeck.de
#SBATCH --job-name=hicpro

###later add #SBATCH --tmp=500GB

PATH=$WORK/.omics/anaconda3/bin:$PATH #add the anaconda installation path to the bash path
source $WORK/.omics/anaconda3/etc/profile.d/conda.sh # some reason conda commands are not added by default

cat hicpro.yml


# Note: version hicpro 3.1.0 
# change line 92 in bowtie_pairing.sh from 
# merge_pairs $sample_dir $R1 $R2 &> ${ldir}/mergeSAM.log
# to 
# merge_pairs $sample_dir $R1 $R2 ${ldir}

## for parallel mode 
snakemake -s hicpro.smk hicpro_parallel_step1 --cores $(nproc) --use-conda --conda-prefix $WORK/hic_conda_envs/
### wait for jobs to finish ?
#snakemake -s hicpro.smk hicpro_parallel_step2 --cores $(nproc) --use-conda 

snakemake -s hicpro.smk pool --cores $(nproc) --use-conda --conda-prefix $WORK/hic_conda_envs/

