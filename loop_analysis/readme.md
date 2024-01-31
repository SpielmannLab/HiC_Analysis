# Peakachu analysis log

## Peakachu installation

Created peakachu env using installation instructions at their [github page](https://github.com/tariks/peakachu#installation):

## Analysis steps

1. Start with juicebox \*.hic files. For peakachu-based loop calling, use \*.hic file cureated using juicer_tools_1.19.02.jar (later verions are not compatible with Peakachu, unless you use a tool like HiCLift to convert to \*.cool or \*.mcool format).

2. Check the model to use:

   sbatch /data/humangen_external/HiC/steinhaeuser_hic/Steinhaeser_HiC_analysis/loop_analysis/checkdepth_batch.sh \*.hic

3. Based on the output, download the \"pickled models\" from the Peakachu github page. For example, 30M and 50M intra-chromosomal contacts at 10kb resolution can be downloaded as follows:

   cd /data/humangen_external/HiC/results/peakachu/picked_models
   wget http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.300million.10kb.w6.pkl
   wget http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.350million.10kb.w6.pkl

4. Run peakachu to call peaks.
   The script here is modification of the snakemake script in the [HiC_Analysis repo](https://github.com/SpielmannLab/HiC_Analysis/tree/Varuns_edits). The configuration file to call the peaks is [here](/Steinhaeser_Hi/HiC_Analysis_params/peakachu.yml). Threshold was chosen to be 0.97 to get confident peaks only. As can be seen the peaks are called without the "--balance" option. This is because peakachu does not seem to be able to read the normalized matrices within the \*.hic file properly. The only solution for this is to convert hic files to .mcool or .cool format using tools like HiCLift. However, upon preliminary testing, Peakachu-based peak calling seems to get worse when using the "--balance" option.

First fill out the [call_loops_params.yaml](Steinhaeser_HiC_analysis/loop_analysis/call_loops_params.yaml) file and then submit the following code from the headnode:

    ```bash
    sbatch call_loops_core.sh
    ```

5. Run diffpeakachu
   This is being done using new NextFlow script [here](/peakachu/nextflow_scripts/).

6. Run HICCUPS
