# Instructions to use fastq2HiC
*Varun K A Sreenivasan, 10th October, 2023*

This pipeline enables the primary analysis of HiC fastq files using [HiC-Pro](https://github.com/nservant/HiC-Pro) and outputs the contact matrix as *.hic file (using JuicerTools Pre) to be used with [Juicebox and JuicerTools](https://github.com/aidenlab/JuicerTools). Most of the internal scripts have been adapted from the previous work by [Kristin Schultz](https://github.com/SpielmannLab/HiC_Analysis.), but now implemented using NextFlow.

The user starts with filling out the [params.yaml file](fastq2HiC/fastq2HiC_params.yaml) and then runs using:

    # submit the following command. It will first copy all the scripts to $WORK/hic_nextflow_launchdir/ and then run the pipeline
    sbatch fastq2HiC_sbatch.sh

# Input data and installation of tools
The following data and tools are required. 

1. Genome reference

Several files are needed for HiC-Pro: 

    - The *.fa file containing full sequence
    - The *fa.fai file containing samtools index
    - The bowtie2 index. Note, bowtie*1* index is not compatible for the alignment

These files are best downloaded from illumina's [iGenomes](https://emea.support.illumina.com/sequencing/sequencing_software/igenome.html). Suggest the UCSC version. Make sure that the file does not containg *"_ALT"* scaffolds as this can draw out mappings and result in reduced number of mapped reads although bowtie2 aligner implemented by HiC-Pro is supposed to be ALT-aware. Better be safe.

Note - If the genome of interest is not available from iGenomes, then use the build_reference_indices process in the defined in the [modules file](fastq2HiC/fastq2HiC_modules.nf).

2. HiC-Pro

Check out the [github page](https://github.com/nservant/HiC-Pro/tree/master), the nice [paper](https://doi.org/10.1186/s13059-015-0831-x), as well as the indepth documentation [here](http://nservant.github.io/HiC-Pro/QUICKSTART.html) on HiC-pro. It comes with several handy scripts that are also used in this pipeline

Note that we not only require the HiC-Pro scripts, obtained by downloading the git repo, but also the dependency conda environment with the name: "HiC_pro" (Make sure that this is the exact name used. If not, change in the [nextflow.config](fastq2HiC/nextflow.config) file accordingly). The yaml file to install this conda environment is available in the HiC-Pro github repo.

3. Juicer tools

**Note**: Juicer tools and Juicebox are very version ustable. For example, the default normalizations have changed. KR is no more included by default and the hicpro2juicebox.sh script only uses the default option. This can be manually edited if needed. As of 10th October 2023, the output of juicer_tools.2.20.00.jar is readable with Juicebox version 2.15. However, these are development versions and may get updates. Peakachu (the peakcaller used downstream of this pipeline) is only compatible until juicer_1.22.01. So, download as many versions as needed and provide the details in the [params.yaml file](fastq2HiC/fastq2HiC_params.yaml) file. 

# Running the fastq2HiC
Fill out the [params.yaml file](fastq2HiC/fastq2HiC_params.yaml) file and run using the following script. You should not have to touch the config-hicpro.txt, since this is automatically edited within the nextflow script. Nor should you have to touch the the nextflow.config file. 

    sbatch fastq2HiC_sbatch.sh

You can follow the progress of the job by watching the slurm-xxx.out file generated in the current directory. The script spawns several sub-slurm jobs. You can follow these using these commands:

    watch -n 10 squeue -u $USER
    watch -n 10 tail -n 30 slurm*.out

# Expected output
When run successfully, several output files are generated in the output directory specified in the fastq2HiC_params.yaml file. The main file that can be opened using JuiceBox is located in "samplename/juicer_xxx/*.hic". You might also want to look at the "samplename/results/hic_results/stats" and "samplename/restuls/hic_results/pic" for quality check and sequencing depth check.

# Debugging when run unsuccesful
Always look at the completion of the jobs using the code below. Sometime a job may not have run succesfull due to memory issues. If so, there would be entries saying "OUT OF MEMORY". Sometimes, a particular FASTQ file (after splitting) does not align well and causes TIMEOUT of the sbatch script in the \"hic_pro_in_parallel_mode\" process. This needs to be looked at - as to why this happens sometimes. Could be a node-specific issue.
   
   sacct -u $USER
   or 
   sacct -u $USER --starttime now-13hours -l | column -ts $'\t' | vi
   or 
   sacct -u $USER --starttime 2023-12-13T15:00 -l | column -ts $'\t' | vi

In some cases, the workflow runs well, except for the process "allValidPairs_to_juicertools_hic". In this case, it is not necessary to rerun everything, but a subworkflow. Ask Varun how to do it.

# Next steps
The next step is the A/B compartment analysis and loop calling. For this either use Peakachu or HICCUPS. Differential loops can be called using DiffPeakachu and HICCUPSDiff.
