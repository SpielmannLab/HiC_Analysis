################################################
#
# This workflow uses the tool CHESS for intra chromosomal analysis of hic data.
# Date: January, 2023
# Author: Kristin Schultz, k.schultz@uni-luebeck.de
#
################################################

import os
    
configfile: "config.yml"

rule all:
	input:
		"%s/summary.txt"%config['OUTDIR']

rule chess_pairs:
	output:
		"%s/regions.bed"%config['OUTDIR']
	params:
		ref=config['REFERENCE_NAME'],
		win=config['WINDOW_SIZE'],
		step=config['STEP_SIZE'],
		prefix1=config['CONTROL_PREFIX'],
		prefix2=config['SAMPLE_PREFIX'],
		scratch=config['SCRATCH']
	conda:
		"envs/chess-hic.yml"
	shell:
		"""
		chess pairs {params.ref} {params.win} {params.step} {params.scratch}/tmp.bed
		sed -i -r "s|chr||g" {params.scratch}/tmp.bed
		awk -v OFS='\\t' '$1="{params.prefix1}"$1, $4="{params.prefix2}"$4' {params.scratch}/tmp.bed > {output}
		rm {params.scratch}/tmp.bed
		"""

rule chess_sim:
	input:
		"%s/regions.bed"%config['OUTDIR']
	output:
		"%s/%s_%s.tsv"%(config['OUTDIR'], config['SAMPLE'], config['CONTROL'])
	params:
		win=config['WINDOW_SIZE'],
		step=config['STEP_SIZE'],
		res=config['RESOLUTION'],
		obsFile=config['SAMPLE_FILE'],
		conFile=config['CONTROL_FILE'],
		norm=config['NORM']
	conda:
		"envs/chess-hic.yml"
	shell:
		"""
		chess sim -p $(nproc) {params.conFile}@{params.res} {params.obsFile}@{params.res} {input} {output}
		"""
		
rule filter_regions:
	input:
		regions="%s/regions.bed"%config['OUTDIR'],
		filtered="%s/%s_%s.tsv"%(config['OUTDIR'], config['SAMPLE'], config['CONTROL'])
	output:
		"%s/filtered_regions_%s.tsv"%(config['OUTDIR'], config['WINDOW_SIZE'])
	params:
		win=config['WINDOW_SIZE'],
		outdir=config['OUTDIR'],
		sn_thr=config['SIM_THRESHOLD'],
		zsim_thr=config['Z_SIM_THRESHOLD']
	conda:
		"envs/chess-hic.yml"
	shell:
		"""
		regions=$(basename {input.regions})
		filtered=$(basename {input.filtered})
		cp scripts/analyzeCHESS_01.py {params.outdir}/
		ANALYZE={params.outdir}/analyzeCHESS_01.py
		sed -i -r "s|wdir = .*|wdir = '{params.outdir}/'|g" $ANALYZE
		sed -i -r "s|winsize =.*|winsize = {params.win}|g" $ANALYZE
		sed -i -r "s|chess_results_file =.*|chess_results_file = '$filtered'|g" $ANALYZE
		sed -i -r "s|region_pairs =.*|region_pairs = '$regions'|g" $ANALYZE
		sed -i -r "s|sn_thr = .*|sn_thr = {params.sn_thr}|g" $ANALYZE
		sed -i -r "s|zsim_thr = .*|zsim_thr = {params.zsim_thr}|g" $ANALYZE
		python $ANALYZE
		"""
		
rule chess_extract:
	input:
		"%s/filtered_regions_%s.tsv"%(config['OUTDIR'], config['WINDOW_SIZE'])
	output:
		gain="%s/features/gained_features.tsv"%config['OUTDIR'],
		loss="%s/features/lost_features.tsv"%config['OUTDIR'],
	params:
		res=config['RESOLUTION'],
		obsFile=config['SAMPLE_FILE'],
		conFile=config['CONTROL_FILE'],
		out="%s/features/"%config['OUTDIR'],
		norm=config['NORM']
	conda:
		"envs/chess-hic.yml"
	shell:
		"""
		chess extract {input} {params.conFile}@{params.res} {params.obsFile}@{params.res}  {params.out}
		"""

rule chess_crosscorrelate:
	input:
		gain="%s/features/gained_features.tsv"%config['OUTDIR'],
		loss="%s/features/lost_features.tsv"%config['OUTDIR'],
		regions="%s/filtered_regions_%s.tsv"%(config['OUTDIR'], config['WINDOW_SIZE'])
	output:
		directory("%s/correlations/"%config['OUTDIR'])
	conda:
		"envs/chess-hic.yml"
	shell:
		"""
		mkdir -p {output}/gain/
		mkdir -p {output}/loss/
		chess crosscorrelate {input.gain} {input.regions} {output}/gain/
		chess crosscorrelate {input.loss} {input.regions} {output}/loss/
		"""

rule plot_results:
	input:
		regions="%s/regions.bed"%config['OUTDIR'],
		filtered="%s/%s_%s.tsv"%(config['OUTDIR'], config['SAMPLE'], config['CONTROL']),
		pairs="%s/%s_%s.tsv"%(config['OUTDIR'], config['SAMPLE'], config['CONTROL']),
		gain="%s/features/gained_features.tsv"%config['OUTDIR']
	output:
		dynamic("%s/plots/{N}_region_with_features.png"%config['OUTDIR'])
	params:
		win=config['WINDOW_SIZE'],
		res=config['RESOLUTION'],
		obsFile=config['SAMPLE_FILE'],
		conFile=config['CONTROL_FILE'],
		numPlots=config['N_PLOTS'],
		outdir=config['OUTDIR'],
		sn_thr=config['SIM_THRESHOLD'],
		zsim_thr=config['Z_SIM_THRESHOLD']
	conda:
		"envs/chess-hic.yml"
	shell:
		"""
		regions=$(basename {input.regions})
		filtered=$(basename {input.filtered})
		cp scripts/analyzeCHESS_02.py {params.outdir}
		ANALYZE={params.outdir}/analyzeCHESS_02.py
		sed -i -r "s|wdir = .*|wdir = '{params.outdir}/'|g" $ANALYZE
		sed -i -r "s|winsize =.*|winsize = {params.win}|g" $ANALYZE
		sed -i -r "s|chess_results_file =.*|chess_results_file = '$filtered'|g" $ANALYZE
		sed -i -r "s|region_pairs =.*|region_pairs = '$regions'|g" $ANALYZE
		sed -i -r "s|sn_thr = .*|sn_thr = {params.sn_thr}|g" $ANALYZE
		sed -i -r "s|zsim_thr = .*|zsim_thr = {params.zsim_thr}|g" $ANALYZE
		sed -i -r "s|patient_hic = fanc.load.*|patient_hic = fanc.load('{params.obsFile}@{params.res}')|g" $ANALYZE
		sed -i -r "s|control_hic = fanc.load.*|control_hic = fanc.load('{params.conFile}@{params.res}')|g" $ANALYZE
		sed -i -r "s|maxPlots = .*|maxPlots = {params.numPlots}|g" $ANALYZE
		python $ANALYZE
		"""

rule report:
	input:
		dynamic("%s/plots/{N}_region_with_features.png"%config['OUTDIR'])
	output:
		"%s/summary.txt"%config['OUTDIR']
	params:
		outdir=config['OUTDIR'],
		win=config['WINDOW_SIZE'],
		step=config['STEP_SIZE'],
		res=config['RESOLUTION'],
		obsFile=config['SAMPLE_FILE'],
		conFile=config['CONTROL_FILE'],
		numPlots=config['N_PLOTS'],
		sn_thr=config['SIM_THRESHOLD'],
		zsim_thr=config['Z_SIM_THRESHOLD'],
		norm=config['NORM']
	conda:
		"envs/chess-hic.yml"
	shell:
		"""
		echo "Chess version" $(chess --version) > {output}
		echo "Processed sample {params.obsFile} (patient) and sample {params.conFile} (control), both {params.norm} normalized." >> {output}
		echo "Run with window size {params.win} and stepsize {params.step} at resolution {params.res}." >> {output}
		echo "Similarity Threshold was set to >=  {params.sn_thr} and z-score was set to <= {params.zsim_thr}."  >> {output}
		echo "Top {params.numPlots} structural differences :"
		cat {params.outdir}/ranking.txt >> {output}
		echo "Run by $USER and finished at $(date)" >> {output}
		"""
	