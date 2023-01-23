################################################
#
# This workflow uses HiCPro to align HiC reads
# and offers to transform them into the .hic format
# Date: October, 2022
# Author: Kristin Schultz, k.schultz@uni-luebeck.de
#
################################################

import os
import glob
configfile: "hicpro.yml"

from random import random
run_nr=int(random()*1000000000)

n_Lanes=config['LANES']
lanes=list()
if n_Lanes:
	lst=list(range(1,n_Lanes+1))
	for i in lst:
		lanes.append('L00'+str(i))

### if rename not true ruse sample names as names
if config['RENAME'] != "TRUE" :
	config['NAMES'] = config['SAMPLES']
	
rule all:
	input:
		directory(expand("%s/{name}/{out}"%config['OUT_DIR'], name=config['NAMES'], out=config['OUTS'])),
		directory(expand("%s/{name}/logs"%config['OUT_DIR'], name=config['NAMES']))
	params:
		scratch=config['SCRATCH']
	shell:
		"""
		rm -r {params.scratch}/*
		"""

### Hicpro configfile will be generated for each run according to input in HiC-Pro.yml
###
rule configuration:
	input:
		sizes="data/%s.sizes"%(config['REFERENCE_NAME']),
		fraq="data/%s_resfraq_%s.bed"%(config['ENZYME'],config['REFERENCE_NAME'])
	output:
		"%s/%s-config.yml"%(config['SCRATCH'],run_nr)
	params:
		ref=config['REFERENCE_NAME'],
		bw_idx=config['REF_IDX_DIR'],
		ligation=config['LIGATION_SITE'],
		dir=config['HICPRO_INSTALL_DIR']
	shell:
		"""
		cp {params.dir}/config-hicpro.txt {output}

		sed -i "s@REFERENCE_GENOME =.*@REFERENCE_GENOME = {params.ref}@g" {output}
		sed -i "s@GENOME_SIZE =.*@GENOME_SIZE = {input.sizes}@g" {output}
		sed -i "s@GENOME_FRAGMENT =.*@GENOME_FRAGMENT = {input.fraq}@g" {output}
		sed -i "s@BOWTIE2_IDX_PATH =.*@BOWTIE2_IDX_PATH = {params.bw_idx}@g" {output}
		sed -i "s@LIGATION_SITE =.*@LIGATION_SITE ={params.ligation}@g" {output}

		sed -i "s@PAIR1_EXT =.*@PAIR1_EXT = _R1@g" {output}
		sed -i "s@PAIR2_EXT =.*@PAIR2_EXT = _R2@g" {output}
		"""

### You may change and add parameters like this
#
# 		sed -i "s/N_CPU =.*/N_CPU = 4/g" {output}
# 		sed -i "s/JOB_NAME =.*/JOB_NAME = $USER/g" {output}
# 		sed -i "s/JOB_QUEUE =.*/JOB_QUEUE = shortterm/g" {output}
# 		sed -i "s/JOB_MEM =.*/JOB_MEM = 8000/g" {output}
# 		sed -i "s/JOB_WALLTIME =.*/JOB_WALLTIME = 1-00:00:00/g" {output}
# 		sed -i "s/SORT_RAM =.*/SORT_RAM = 4000/g" {output}
###


### link reference to local 
### and create index
rule link_ref:
	params:
		path=config['REFERENCE_FILE'],
		ref=config['REFERENCE_NAME']
	output:
		fa="data/ref/%s.fa"%(config['REFERENCE_NAME']),
		idx="data/ref/%s.fa.fai"%(config['REFERENCE_NAME'])
	shell:
		"""
		if [ ! -e {output.fa} ] ; then
			cp {params.path} {output.fa}
		fi

		samtools faidx -o {output.idx} {output.fa}
		"""
### create fragment file (reference fragmented by config enzyme)
### hicpro uses this file to check for valid pairs
### index .fai must be in same directory as reference
rule ref_fragments:
	input:
		config['REFERENCE_FILE'],
	output:
		"data/%s_resfraq_%s.bed"%(config['ENZYME'],config['REFERENCE_NAME'])
	params:
		motif=config['RES_MOTIF'],
		dir=config['HICPRO_INSTALL_DIR']
	conda: "envs/HiC-Pro.yml"
	shell:
		"""
		python {params.dir}/bin/utils/digest_genome.py -r {params.motif} -o {output} {input}
		"""

### create file that lists the size of each chromosome
rule ref_sizes:
	input:
		"data/ref/%s.fa.fai"%(config['REFERENCE_NAME'])
	output:
		"data/%s.sizes"%(config['REFERENCE_NAME'])
	shell:
		"""
		cut -f1-2 {input} > {output}
		"""
### create soft link to fastq file
### the link follows the pipeline nomenclature, while the fastq may contain additional strings
rule link_files:
	input:
		R1=lambda wildcards: glob.glob("%s/%s*%s_R1*.f*q.gz"%(config['SAMPLES_PATH'], wildcards.sample, wildcards.lane)),
		R2=lambda wildcards: glob.glob("%s/%s*%s_R2*.f*q.gz"%(config['SAMPLES_PATH'], wildcards.sample, wildcards.lane)),
	output:
		R1="%s/rawData/{sample}-{lane}_R1.fastq.gz"%config['SCRATCH'],
		R2="%s/rawData/{sample}-{lane}_R2.fastq.gz"%config['SCRATCH']
	params:
		scratch=config['SCRATCH']
	shell:
		"""
		for f in {input.R1} ; do
			ln -s $f {params.scratch}/rawData/{wildcards.sample}-{wildcards.lane}_R1.fastq.gz
		done
		for f in {input.R2} ; do
			ln -s $f {params.scratch}/rawData/{wildcards.sample}-{wildcards.lane}_R2.fastq.gz
		done
		"""

### rename the softlinks 
rule all_rename_R1:
	input:
		expand("%s/samples/{name}/{name}-{lane}_R1.fq.gz"%config['SCRATCH'], name=config['NAMES'], lane=lanes)

rule rename_R1:
	input:
		expand("%s/rawData/{sample}-{{lane}}_R1.fastq.gz"%config['SCRATCH'], sample=config['SAMPLES'])
	output:
		expand("%s/samples/{name}/{name}-{{lane}}_R1.fq.gz"%config['SCRATCH'], name=config['NAMES'])
	params:
		name=config['NAMES'],
		samples=config['SAMPLES'],
		scratch=config['SCRATCH']
	shell:
		"""
		file_array=({params.samples})
		filepath_array=({input})
		name_array=({params.name})

		count=${{#file_array[@]}}
		for i in `seq 1 $count` ; do
			echo "RENAME "${{filepath_array[$i-1]}}" to "{params.scratch}/samples/${{name_array[$i-1]}}/${{name_array[$i-1]}}-{wildcards.lane}_R1.fq.gz
			mv ${{filepath_array[$i-1]}} {params.scratch}/samples/${{name_array[$i-1]}}/${{name_array[$i-1]}}-{wildcards.lane}_R1.fq.gz
		done

		"""

rule all_rename_R2:
	input:
		expand("%s/samples/{name}/{name}-{lane}_R2.fq.gz"%config['SCRATCH'], name=config['NAMES'], lane=lanes)

rule rename_R2:
	input:
		expand("%s/rawData/{sample}-{{lane}}_R2.fastq.gz"%config['SCRATCH'], sample=config['SAMPLES'])
	output:
		expand("%s/samples/{name}/{name}-{{lane}}_R2.fq.gz"%config['SCRATCH'], name=config['NAMES'])
	params:
		name=config['NAMES'],
		samples=config['SAMPLES'],
		scratch=config['SCRATCH']
	shell:
		"""
		file_array=({params.samples})
		filepath_array=({input})
		name_array=({params.name})

		count=${{#file_array[@]}}
		for i in `seq 1 $count` ; do
			echo "RENAME "${{filepath_array[$i-1]}}" to "{params.scratch}/samples/${{name_array[$i-1]}}/${{name_array[$i-1]}}-{wildcards.lane}_R2.fq.gz
			mv ${{filepath_array[$i-1]}} {params.scratch}/samples/${{name_array[$i-1]}}/${{name_array[$i-1]}}-{wildcards.lane}_R2.fq.gz
		done
		"""

### use hic-pro script to split fastq files into smaller fastq files (used for parallelization) 
rule all_split_sample:
	input:
		directory(expand("%s/splits/{name}/"%config['SCRATCH'], name=config['NAMES']))

rule split_sample:
	input:
		expand("%s/samples/{{name}}/{{name}}-{lane}_R1.fq.gz"%config['SCRATCH'], lane=lanes),
		expand("%s/samples/{{name}}/{{name}}-{lane}_R2.fq.gz"%config['SCRATCH'], lane=lanes)
	output:
		directory("%s/splits/{name}/"%config['SCRATCH'])
	params: dir=config['HICPRO_INSTALL_DIR']
	conda: "envs/HiC-Pro.yml"
	shell:
		"""
		for f in {input}
		do
			python {params.dir}/bin/utils/split_reads.py --results_folder {output} $f
		done
		"""

### run HiC-Pro in parallel mode (-p flag)
### this creates inputfiles_$USER.txt, which lists all input files
### using sed command to change inputfiles_$USER.txt because the paths are always messed up
### this may be fixed in later HiC-Pro versions 
### HiC-Pro also generates two scripts: HiCPro_step1_$USER.sh and HiCPro_step2_$USER.sh
### starting HiCPro_step1_$USER.sh in this rule
rule hicpro_parallel_step1:
	input:
		directory(expand("%s/splits/{name}/"%config['SCRATCH'], name=config['NAMES'])),
		conf="%s/%s-config.yml"%(config['SCRATCH'],run_nr)
	output:
		directory(expand("%s/results/bowtie_results/bwt2/{name}/"%config['SCRATCH'], name=config["NAMES"])),
		data=directory(expand("%s/results/hic_results/data/{name}/"%config['SCRATCH'], name=config["NAMES"])),
		logs=directory(expand("%s/results/logs/{name}/"%config['SCRATCH'], name=config["NAMES"]))
	conda: "envs/HiC-Pro.yml"
	params:
		scratch=config['SCRATCH'],
		dir=config['HICPRO_INSTALL_DIR']
	threads: 1
	shell:
		"""
		set +o pipefail
		PATH={params.dir}/bin:$PATH
		yes | HiC-Pro -i {params.scratch}/splits -o {params.scratch}/results -c {input.conf} -p
		sed -i "s@^@../splits/@" {params.scratch}/results/inputfiles_$USER.txt
		cd {params.scratch}/results
		srun HiCPro_step1_$USER.sh
		"""

### run the 2nd step HiCPro_step2_$USER.sh
rule hicpro_parallel_step2:
	input:
		directory(expand("%s/results/bowtie_results/bwt2/{name}/"%config['SCRATCH'], name=config["NAMES"])),
	output:
		plots=directory(expand("%s/results/hic_results/pic/{name}/"%config['SCRATCH'], name=config["NAMES"])),
		stats=directory(expand("%s/results/hic_results/stats/{name}/"%config['SCRATCH'], name=config["NAMES"])),
		raw=directory(expand("%s/results/hic_results/matrix/{name}/iced"%config['SCRATCH'], name=config["NAMES"])),
		iced=directory(expand("%s/results/hic_results/matrix/{name}/raw"%config['SCRATCH'], name=config["NAMES"]))
	params:
		scratch=config['SCRATCH']
	shell:
		"""
		set +o pipefail
		cd {params.scratch}/results
		srun HiCPro_step2_$USER.sh
		"""

### HiC-pro results are stored in "batches" (each HiC-Pro run creates hic_results with all samples)
### moving files to sample based folder structure
# OUTS:
#     - "matrix"
#     - "bowtie_results"
#     - "hic_results" <<<<<
#     - "hic_format"
rule move_hicresults:
	input:
		plots=directory("%s/results/hic_results/pic/{name}/"%config['SCRATCH']),
		stats=directory("%s/results/hic_results/stats/{name}/"%config['SCRATCH']),
		data=directory("%s/results/hic_results/data/{name}/"%config['SCRATCH']),
	output:
		directory("%s/{name}/hic_results"%config['OUT_DIR'],)
	params:
		scratch=config['SCRATCH']
	shell:
		"""
		mkdir -p {output}/pic/
		mkdir -p {output}/stats/
		mkdir -p {output}/data/
		mv -v {input.plots}* {output}/pic/
		mv -v {input.stats}* {output}/stats/
		mv -v {input.data}* {output}/data/
		"""

rule move_logs:
	input:
		logs=directory("%s/results/logs/{name}/"%config['SCRATCH'])
	output:
		directory("%s/{name}/logs"%config['OUT_DIR'],)
	shell:
		"""
		mkdir -p {output}
		mv -v {input}* {output}
		"""

# OUTS:
#     - "matrix"
#     - "bowtie_results" <<<<<
#     - "hic_results"
#     - "hic_format"
rule move_bams:
	input:
		directory("%s/results/bowtie_results/bwt2/{name}/"%config['SCRATCH'])
	output:
		directory("%s/{name}/bowtie_results"%config['OUT_DIR'],)
	shell:
		"""
		mkdir -p {output}
		mv -v {input}* {output}/
		"""
		
# OUTS:
#     - "matrix" <<<<<
#     - "bowtie_results"
#     - "hic_results"
#     - "hic_format"
rule move_icematrix:
	input:
		raw=directory("%s/results/hic_results/matrix/{name}/iced"%config['SCRATCH']),
		iced=directory("%s/results/hic_results/matrix/{name}/raw"%config['SCRATCH'])
	output:
		directory("%s/{name}/matrix"%config['OUT_DIR'])
	shell:
		"""
		mkdir -p {output}
		mv -v {input.raw} {output}/
		mv -v {input.iced} {output}/
		"""

### pools ALL SAMPLES into one
### remove samples that should not be merged from configfile hicpro.yml
rule pool:
	input:
		expand(directory("%s/{name}/hic_results"%config['OUT_DIR'],), name=config['NAMES'])
	output:
		"%s/%s/cool_format/%s.5000.cool"%(config['OUT_DIR'],config['POOL'], config['POOL'])
	params:
		dir=config['HICPRO_INSTALL_DIR'],
		sizes="data/%s.sizes"%(config['REFERENCE_NAME']),
		scratch="%s"%config['SCRATCH']
	conda: "envs/Cooler.yml"
	shell:
		"""
		cool_files=""
		for path in {input} ; do
			vp=$path/data/*.allValidPairs
			name=$(basename $vp .allValidPairs)
			mkdir -p {params.scratch}/$name/
			bash {params.dir}/hicpro2higlass.sh -i $vp -r 5000 -c {params.sizes} -n -o {params.scratch}/$name/
			cool_files=$cool_files" "$(ls {params.scratch}/$name/*.cool)
		done

		echo "Merge "$cool_files" to "{output}"."
		cooler merge {output} $cool_files
		"""

rule cool2mcool:
	input:
		"%s/%s/cool_format/%s.5000.cool"%(config['OUT_DIR'],config['POOL'], config['POOL'])
	output:
		"%s/%s/cool_format/%s.mcool"%(config['OUT_DIR'], config['POOL'], config['POOL'])
	params:
		resolutions="5000,10000,25000,50000,100000,500000,10000000,2500000"
	conda: "envs/Cooler.yml"
	threads: 4
	shell:
		"""
		cooler zoomify -o {output} -n {threads} -r {params.resolutions} ${out}-both_wt.$res.cool {input}
		"""

rule cool2hic:
	input:
		cool="%s/%s/cool_format/%s.5000.cool"%(config['OUT_DIR'], config['POOL'], config['POOL'])
	output:
		"%s/%s/hic_format/%s.hic"%(config['OUT_DIR'],config['POOL'], config['POOL'])
	params:
		sizes="data/%s.sizes"%(config['REFERENCE_NAME']),
		fraq="data/%s_resfraq_%s.bed"%(config['ENZYME'],config['REFERENCE_NAME']),
		scratch="%s"%config['SCRATCH'],
		resolutions="5000,10000,25000,50000,100000,500000,10000000,2500000",
		name=config['POOL']
	conda: "envs/HiCExplorer.yml"
	threads: 4
	shell:
		"""
		hicConvertFormat -m {input} -o {params.scratch}/{params.name}.ginteractions --inputFormat cool --outputFormat ginteractions
		awk -F "\\t" '{{print 0, $1, $2, 0, 0, $4, $5, 1, $7}}'  {params.scratch}/{params.name}.ginteractions >  {params.scratch}/{params.name}.ginteractions.tsv.short
		sort -k2,2d -k6,6d  {params.scratch}/{params.name}.ginteractions.tsv.short >  {params.scratch}/{params.name}.ginteractions.tsv.short.sorted
		sed -r "s|chr||g" {params.scratch}/{params.name}.ginteractions.tsv.short.sorted
		TMP={params.scratch}/tmp
		mkdir -p $TMP
		export _JAVA_OPTIONS="-Xmx79g -Xms40g"
		java -Xmx79g -Xms40g -jar $juicer pre -r {params.resolution} -j {threads} -t $TMP  -f {params.fraq} {params.scratch}/{params.name}.ginteractions.tsv.short.sorted {output} {params.sizes}
		"""


rule all_create_hic_from_validpairs:
	input:
		directory(expand("%s/{name}/hic_format"%config['OUT_DIR'], name=config['NAMES']))

rule create_hic_from_validpairs:
	input:
		dir=directory("%s/{name}/hic_results"%config['OUT_DIR']),
		sizes="data/%s.sizes"%(config['REFERENCE_NAME'])
	output:
		tmp=directory("%s/validPairs2hic/{name}"%config['SCRATCH']),
		save=directory("%s/{name}/hic_format"%config['OUT_DIR'],)
	params:
		ref=config['REFERENCE_NAME'],
		tools=config['JUICER_PATH'],
		scratch=config['SCRATCH'],
		res=config['RESOLUTION'],
		dir="%s/scripts/"%config['HICPRO_INSTALL_DIR']
	threads: 6
	shell:
		"""
		mkdir -p {output.save}
		mkdir -p {output.tmp}
		mkdir -p {params.scratch}/tmp/
		bash {params.dir}/bin/utils/hicpro2juicebox.sh -o {output.tmp}/ -t {params.scratch}/tmp/ -i {input.dir}/data/{wildcards.name}.allValidPairs -g {input.sizes} -j {params.tools}
		wait
		mv {output.tmp}/{wildcards.name}.allValidPairs.hic {output.tmp}/{wildcards.name}_{params.ref}.hic
		java -jar {params.tools} addNorm -j {threads} -w {params.res} {output.tmp}/{wildcards.name}_{params.ref}.hic
		cp -r {output.tmp}/{wildcards.name}_{params.ref}.hic {output.save}
		"""

rule create_hic_noChr:
	input:
		dir=directory("%s/{name}/hic_results"%config['OUT_DIR'],),
		sizes="data/%s.sizes"%(config['REFERENCE_NAME'])
	output:
		tmp=directory("%s/validPairs2hic/{name}"%config['SCRATCH']),
		save="%s/{name}/hic_format/{name}_%s.noChr.hic"%(config['OUT_DIR'],config['REFERENCE_NAME'])
	params:
		ref=config['REFERENCE_NAME'],
		tools=config['JUICER_PATH'],
		scratch=config['SCRATCH'],
		res=config['RESOLUTION'],
		dir="%s/scripts/"%config['HICPRO_INSTALL_DIR']
	shell:
		"""
		mkdir -p {output.tmp}
		sed -r "s|chr||g" {input.dir}/data/{wildcards.name}.allValidPairs > {output.tmp}/{wildcards.name}.fixed.allValidPairs
		sed -r "s|chr||g" {input.sizes} > {output.tmp}/noChr.sizes

		bash {params.dir}/bin/utils/hicpro2juicebox.sh -o {output.tmp}/ -t {params.scratch}/tmp/ -i {output.tmp}/{wildcards.name}.fixed.allValidPairs -g {output.tmp}/noChr.sizes -j {params.tools}
		wait
		mv {output.tmp}/{wildcards.name}.fixed.allValidPairs.hic {output.tmp}/{wildcards.name}_{params.ref}.noChr.hic
		java -jar {params.tools} addNorm -j $(nproc) -w {params.res} {output.tmp}/{wildcards.name}_{params.ref}.noChr.hic
		cp -r {output.tmp}/{wildcards.name}_{params.ref}.noChr.hic {output.save}
		"""
