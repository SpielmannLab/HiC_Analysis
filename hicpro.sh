#! /bin/bash

#SBATCH --partition=shortterm
#SBATCH -c 24
#SBATCH --mem=200GB
#SBATCH --mail-type=END
#SBATCH --job-name=hicpro
##SBATCH --tmp=40GB

conda activate HiC_pro

PATH=$WORK/.omics/anaconda3/bin:$PATH #add the anaconda installation path to the bash path
source $WORK/.omics/anaconda3/etc/profile.d/conda.sh # some reason conda commands are not added by default

cat hicpro.yml


# Note: version hicpro 3.1.0 
# change line 92 in bowtie_pairing.sh from 
# merge_pairs $sample_dir $R1 $R2 &> ${ldir}/mergeSAM.log
# to 
# merge_pairs $sample_dir $R1 $R2 ${ldir}


## for parallel mode 
## snakemake -s hicpro.smk hicpro_parallel_step1 --cores $(nproc) --use-conda --conda-prefix $WORK/hic_conda_envs/
### wait for jobs to finish ?
## snakemake -s hicpro.smk hicpro_parallel_step2 --cores $(nproc) --use-conda --conda-prefix $WORK/hic_conda_envs/

# or just run it all.
snakemake -s hicpro.smk all --cores $(nproc) --use-conda --conda-prefix $WORK/hic_conda_envs/ #--config SCRATCH=$SCRATCH #Check if $SCRATCH works.

# Do QC #Does not work due to some kind of buffer error.
# snakemake -s hicpro.smk hicrep --cores $(nproc) --use-conda --conda-prefix $WORK/hic_conda_envs/ --config SCRATCH=$SCRATCH

# For pooling multiple samples together.
#snakemake -s hicpro.smk pool --cores $(nproc) --use-conda --conda-prefix $WORK/hic_conda_envs/

