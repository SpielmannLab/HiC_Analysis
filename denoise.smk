################################################
#
# Denoise with HiCorr and DeepLoop. 
# Date: January 2023
# Author: Kristin Schultz, k.schultz@uni-luebeck.de
#
################################################

import os

configfile: "denoise.yml"

lib="HiCorr/bin/preprocess"

try:
	os.mkdir("%s/logs/"%(resultdir))
	print("Directory " , "%s/logs/"%(resultdir) ,  " created.") 
except FileExistsError:
	print("Directory " , "%s/logs/"%(resultdir) ,  " already exists.")

resultdir=config['OUT_DIR']

### list of all chromosomes to consider
CHR=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']	

### create denoise output data for all samples
rule all:
	input:
		expand("%s/HiCorr_DeepLoop_comb/chr{chrom}.raw_HiCorr_DeepLoop"%(resultdir), sample=config['SAMPLES'], chrom=CHR)

### replication of documentation at https://github.com/JinLabBioinfo/HiCorr/blob/master/documents/validPairs_2_fragloop.md
rule validpairs2fragmentpairs:
	input:
		"%s/{sample}/hic_results/data/{sample}.allValidPairs"%(config['HICPRO_INDIR'])
	output:
		"%s/frag_loop.{sample}.trans"%("results/{sample}/fragmentpairs"),
		"%s/frag_loop.{sample}.cis"%("results/{sample}/fragmentpairs")
	params:
		lib=lib,
		fragbed=config['fragbed'],
		loopbed=config['DeepLoopBed'],
		dir="results/{sample}/fragmentpairs",
		scratch=config['SCRATCH'],
		readlength=50
	threads: 3
	shell:
		"""
		mkdir -p {params.dir}
		cat {input} | cut -f2-7 | perl {params.lib}/reads_2_cis_frag_loop.pl {params.fragbed} {params.readlength} {params.dir}/{wildcards.sample}.loop.inward {params.dir}/{wildcards.sample}.loop.outward {params.dir}/{wildcards.sample}.loop.samestrand {params.dir}/{wildcards.sample}.summary.frag_loop.read_count - &
		cat {input} | cut -f2-7 | perl {params.lib}/reads_2_trans_frag_loop.pl {params.fragbed} {params.readlength} {params.dir}/{wildcards.sample}.loop.trans - 

		origin=$(pwd)
		cd {params.dir}/
		for file in {wildcards.sample}.loop.outward {wildcards.sample}.loop.samestrand; do
				cat $file | $origin/{params.lib}/summary_sorted_frag_loop.pl $origin/{params.loopbed}  > temp.$file &
		done
		wait
		for file in {wildcards.sample}.loop.inward {wildcards.sample}.loop.outward {wildcards.sample}.loop.samestrand; do
				perl $origin/{params.lib}/resort_by_frag_id.pl $origin/{params.loopbed} temp.$file &
		done
		wait
		cd $origin
		
		cat {params.dir}/{wildcards.sample}.loop.trans | {params.lib}/summary_sorted_trans_frag_loop.pl - > {params.dir}/temp.{wildcards.sample}.loop.trans
		
		perl {params.lib}/merge_sorted_frag_loop.pl {params.dir}/temp.{wildcards.sample}.loop.samestrand <(cat {params.dir}/temp.{wildcards.sample}.loop.inward | awk '{{if($4>1000)print $0}}') <(cat {params.dir}/temp.{wildcards.sample}.loop.outward | awk '{{if($4>5000)print $0}}') > {params.dir}/frag_loop.{wildcards.sample}.cis &
		perl {params.lib}/merge_sorted_frag_loop.pl {params.dir}/temp.{wildcards.sample}.loop.trans > {params.dir}/frag_loop.{wildcards.sample}.trans 
		wait
		cat {params.dir}/frag_loop.{wildcards.sample}.cis <(cat {params.dir}/frag_loop.{wildcards.sample}.cis | awk '{{print $2 "\t" $1 "\t" $3 "\t" $4}}') | sed s/"frag_"//g | sort -k1,2n -k2,2n | awk '{{print "frag_"$1 "\t" "frag_"$2 "\t" $3 "\t" $4}}' > {params.dir}/frag_loop.{wildcards.sample}.cis.tmp &
		cat {params.dir}/frag_loop.{wildcards.sample}.trans <(cat {params.dir}/frag_loop.{wildcards.sample}.trans | awk '{{print $2 "\t" $1 "\t" $3 "\t" $4}}') | sed s/"frag_"//g | sort -k1,2n -k2,2n | awk '{{print "frag_"$1 "\t" "frag_"$2 "\t" $3 "\t" $4}}' > {params.dir}/temp.frag_loop.{wildcards.sample}.trans 
		wait
		"""


rule HiCorr_all:
	input:
		HiCorr_output=expand("%s/HiCorr_output/anchor_2_anchor.loop.chr{chrom}"%(resultdir), chrom=CHR, sample=config['SAMPLES'])

### all HiCorr output goes to folder HiCorr_output so we can only run samples in parallel, when initiating HiCorr from a different directory.
### Since HiCorr works with relative paths we're copying all scripts to scratch/sample temp directory and start from there.
rule HiCorr:
	input:
		trans="%s/frag_loop.{sample}.trans"%("results/{sample}/fragmentpairs"),
		cis="%s/frag_loop.{sample}.cis"%("results/{sample}/fragmentpairs")
	output:
		HiCorr_output=expand("%s/HiCorr_output/anchor_2_anchor.loop.chr{chrom}"%("results/{{sample}}"), chrom=CHR)
	params:
		enzyme=config['enzyme'],
		genome=config['genome'],
		scratch=config['SCRATCH']
	shell:
		"""
		origin=$(pwd)
		mkdir -p {params.scratch}/{wildcards.sample}/data
		cp -r HiCorr {params.scratch}/{wildcards.sample}/
		cp -r data/HiCorr {params.scratch}/{wildcards.sample}/data/
		cd {params.scratch}/{wildcards.sample}/HiCorr
		bash HiCorr {params.enzyme} $origin/{input.cis} $origin/{input.trans} {wildcards.sample} {params.genome}
		mkdir -p $origin/results/{wildcards.sample}/HiCorr_output/
		mv HiCorr_output/* $origin/results/{wildcards.sample}/HiCorr_output/
		"""

### parallelization for each chromosome
rule DeepLoopGenome:
	input:
		expand("%s/DeepLoopOutput/chr{chrom}.denoised.anchor.to.anchor"%(resultdir), chrom=CHR, sample=config['SAMPLES'])

### DeepLoop requires gpu access so this rule creates a new script and sbatches
rule DeepLoop:
	input:
		"%s/HiCorr_output/anchor_2_anchor.loop.chr{chrom}"%(resultdir)
	output:
		DL_output=("%s/DeepLoopOutput/chr{chrom}.denoised.anchor.to.anchor"%(resultdir))
	params:
		lib=lib,
		DeepLoopPath=config['DeepLoopPath'],
		enzyme=config['enzyme'],
		genome=config['genome'],
		DL_output=("%s/DeepLoopOutput/"%(resultdir)),
		scratch=config['SCRATCH'],
		dir=resultdir
	conda: "envs/DeepLoop.yml"
	log: "%s/logs/HiCorr_output_{sample}_chr{chrom}.log"%(resultdir)
	threads: 1
	shell:
		"""
		origin=$(pwd)
		mkdir -p {params.scratch}
		cd {params.scratch}
		script={wildcards.sample}_chr{wildcards.chrom}.slurm
		
		echo '#!/bin/bash' >$script
		echo "PATH=$""PATH:/home/schultz/miniconda3/bin" >>$script
		echo "source activate DeepLoop" >>$script
		echo "source /etc/profile.d/modules.sh; module load nvidia-cuda/11.1-native;" >> $script
		echo "module load R/3.5.1" >>$script
		echo "module load samtools/v1.3.1 " >>$script
		echo "module load bowtie2/v2.3.0" >>$script
		echo "python3 $origin/DeepLoop/prediction/predict_chromosome.py --full_matrix_dir $origin/{params.dir}/HiCorr_output/ \
			--input_name $origin/{input} \
			--h5_file $origin/{params.DeepLoopPath}/CPGZ_trained/50M.h5 \
			--out_dir $origin/{params.DL_output} \
			--anchor_dir $origin/{params.DeepLoopPath}/ref/{params.genome}_{params.enzyme}_anchor_bed/  \
			--chromosome chr{wildcards.chrom} --small_matrix_size 128  --step_size 128 --dummy 5
			" >>$script
		
		sbatch  --gres=gpu:1 --time=0-00:20:00 --partition=shortterm -J "DeepLoop_"{wildcards.chrom}{wildcards.sample} $script
		sleep 22m
		"""		

rule CombineGenome:
	input:
		expand("%s/HiCorr_DeepLoop_comb/chr{chrom}.raw_HiCorr_DeepLoop"%(resultdir), chrom=CHR, sample=config['SAMPLES'])

### combine HiCorr and DeepLoop output 
rule combine:
	input:
		HiCorr_output="%s/HiCorr_output/anchor_2_anchor.loop.chr{chrom}"%(resultdir),
		DL_output=("%s/DeepLoopOutput/chr{chrom}.denoised.anchor.to.anchor"%(resultdir))
	output:
		"%s/HiCorr_DeepLoop_comb/chr{chrom}.raw_HiCorr_DeepLoop"%(resultdir)
	params:
		lib=lib,
		dir="%s/HiCorr_DeepLoop_comb"%(resultdir)
	conda: "envs/DeepLoop.yml"
	shell:
		"""
		origin=$(pwd)
		mkdir -p {params.dir}
		perl $origin/DeepLoop/OtherAnalysis_in_DeepLoop_paper/lib/pair.HiCorr_DeepLoop.pl  {input.HiCorr_output} {input.DL_output}  5 > {output}
		"""

### create target regions 125 bins (pixel) in each direction starting from target genes
rule target_regions:
	output:
		"%s/Plots/target_regions.txt"%(resultdir)
	params: 
		genes=config['TARGET_GENES'],
		list=config['GENELIST'],
		res=config['RESOLUTION']
	shell:
		"""
		for gene in {params.genes} ; do
			chr=$(cat {params.list} | grep -w $gene | head -1 | awk 'BEGIN {{FS="\\t"}} {{print $1}}') || true
			if [[ $chr == "" ]] ; then
				echo $gene" not found in "{params.list}
				break
			fi
			genestart=$(cat {params.list} | grep -w $gene | head -1 | awk 'BEGIN {{ FS="\\t"}} {{print $2}}')
			win_start=$(($genestart - 125 * {params.res}))
			win_start_kb=$((($win_start / {params.res}) *{params.res} ))
			win_end=$(($genestart + 125 * {params.res}))
			win_end_kb=$((($win_end / {params.res}) *{params.res} ))
			chr=$(echo $chr | sed 's/chr//g')
			echo -e $chr'\\t'$win_start_kb'\\t'$win_end_kb >> {output}
		done			
		"""



rule plot_all:
	input:
		expand("%s/Plots/{gene}.png"%(resultdir), gene=config['TARGET_GENES'], sample=config['SAMPLES'])

rule plot:
	input: 
		expand("%s/HiCorr_DeepLoop_comb/chr{chrom}.raw_HiCorr_DeepLoop"%("results/{{sample}}"), chrom=CHR),
		list="%s/Plots/target_regions.txt"%(resultdir)
	output:
		expand("%s/Plots/{gene}.png"%("results/{{sample}}"), gene=config['TARGET_GENES'])
	params:
		dir="%s/Plots"%(resultdir),
		lib=lib,
		genome=config['genome'],
		enzyme=config['enzyme'],
		genes=config['TARGET_GENES']
	shell:
		"""
		IFS=' ' read -r -a genes <<< "{params.genes}"
		i=0
		while read line; do
			while IFS=$'\\t' read -r -a tmp ; do
  				chr="${{tmp[0]}}"
  				start="${{tmp[1]}}"
    			end="${{tmp[2]}}"
    			
    			perl HiCorr/bin/generate.raw.expt.ratio.matrix.pl data/DeepLoop/DeepLoop_models/ref/{params.genome}_{params.enzyme}_anchor_bed/chr$chr.bed results/{wildcards.sample}/HiCorr_DeepLoop_comb/chr$chr.raw_HiCorr_DeepLoop chr$chr $start $end {params.dir}/${{genes[$i]}}.${{chr}}_${{start}}_${{end}}	
				
				for f in expt.matrix  ratio.matrix  raw.matrix ; do 
					Rscript --vanilla HiCorr/bin/plot.heatmap.r {params.dir}/${{genes[$i]}}.${{chr}}_${{start}}_${{end}}.$f
				done 
				
				convert $(ls {params.dir}/${{genes[$i]}}.${{chr}}_${{start}}_${{end}}*.png) +append {params.dir}/${{genes[$i]}}.png
				
				i=$(($i+1))
    			
    		done <<< $line	
    	done <{input.list}
    	  	
		"""
		
		
