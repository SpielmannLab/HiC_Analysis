###############################################
#
# 
# Hi-C SV analysis with hic_breakfinder and neoloop
# Date: December, 2022
# 
################################################


import os 
from glob import glob
configfile: "neoloopFinder.yml"

hicpro_dir="results"

rule all:
	input:
		expand("results/{sample}/neoloopfinder/neoTads", sample=config['SAMPLES']),
		expand("results/{sample}/neoloopfinder/neoLoops", sample=config['SAMPLES'])

rule plot:
	input:
		directory(expand("results/{sample}/neoloopfinder/plots",sample=config['SAMPLES']))

rule all_cool:
	input:
		expand("results/{sample}/cool_format/{sample}_%s_%s.cool"%(config['REFERENCE'], config['RESOLUTION']), sample=config['SAMPLES'])


rule hic2cool:
	input: "%s/{sample}/hic_format/{sample}_%s.hic"%(hicpro_dir,config['REFERENCE'])
	output:
		"results/{sample}/cool_format/{sample}_%s_%s.cool"%(config['REFERENCE'], config['RESOLUTION'])
	params:
		scratch=config['SCRATCH'],
		res=config['RESOLUTION']
	conda:
		"envs/hic2cool.yml"
	threads: 4
	shell:
		"""
		mkdir -p {params.scratch}/{wildcards.sample}
		cp $input {params.scratch}/{wildcards.sample}/
		hic2cool convert {params.scratch}/{wildcards.sample}/$(basename {input}) {params.scratch}/{wildcards.sample}/tmp -p {threads} -r {params.res}
		mv {params.scratch}/tmp_{params.res}.cool {output}
		"""

rule calculate_cnv:
	input:
		"results/{sample}/cool_format/{sample}_%s_%s.cool"%( config['REFERENCE'], config['RESOLUTION'])
	output:
		cnv="results/{sample}/neoloopfinder/cnv",
		log="results/{sample}/neoloopfinder/logs/cnv-calculation.log"
	params: 
		ref=config['REFERENCE'],
		enz=config['ENZYME']
	conda:
		"envs/neoloopfinder.yml"
	shell:
		"""
		mkdir -p cache
		calculate-cnv -H {input} -g {params.ref} -e {params.enz} --output {output.cnv} \
		--cachefolder cache --logFile {output.log}
		"""

rule segment_cnv:
	input:
		cnv="results/{sample}/neoloopfinder/cnv"
	output:
		sCNVs="results/{sample}/neoloopfinder/seg-cnv",
		logSeg="results/{sample}/neoloopfinder/logs/cnv-seg.log"
	params: 
		res=config['RESOLUTION']
	conda:
		"envs/neoloopfinder.yml"
	threads: 2
	shell:
		"""
		segment-cnv --cnv-file {input.cnv} --binsize {params.res} --output {output.sCNVs} \
		--nproc {threads} --logFile {output.logSeg}
		"""

rule correct_cnv:
	input:
		cool="results/{sample}/cool_format/{sample}_%s_%s.cool"%(config['REFERENCE'], config['RESOLUTION']),
		sCNVs="results/{sample}/neoloopfinder/seg-cnv"
	output:
		cnvcool="results/{sample}/cool_format/{sample}_CNVnorm_%s_%s.cool"%(config['REFERENCE'], config['RESOLUTION'])
	log:
		"results/{sample}/neoloopfinder/logs/correct-cnv.log"
	conda:
		"envs/neoloopfinder.yml"
	params:
		scratch=config['SCRATCH']
	threads: 4
	shell:
		"""
		cp {input.cool} {params.scratch}/{wildcards.sample}.tmp.cool
		cooler balance -p {threads} --force {params.scratch}/{wildcards.sample}.tmp.cool
		correct-cnv -H {params.scratch}/{wildcards.sample}.tmp.cool --cnv-file {input.sCNVs} --nproc {threads} -f --logFile {log}
		cp {params.scratch}/{wildcards.sample}.tmp.cool {output.cnvcool}
		"""

rule preprocess_breaks:
	input:
		"results/{sample}/hic_breakfinder/{sample}.breaks.txt"
	output:
		"results/{sample}/neoloopfinder/{sample}.breaks.txt"
	conda:
		"envs/neoloopfinder.yml"
	shell:
		"""
		python ../utils/prepare-SV-breakpoints.py {input} {output}
		"""

rule assemble_complexSV:
	input:
		sv_bp="results/{sample}/neoloopfinder/{sample}.breaks.txt",
		cnvcool="results/{sample}/cool_format/{sample}_CNVnorm_%s_%s.cool"%(config['REFERENCE'], config['RESOLUTION']),
		log="results/{sample}/neoloopfinder/logs/correct-cnv.log"
	output:
		sv="results/{sample}/neoloopfinder/complexSVs/{sample}.assemblies.txt"
	log:
		"results/{sample}/neoloopfinder/logs/assembleSVs.log"
	params:
		dir="results/{sample}/neoloopfinder/complexSVs/{sample}"
	conda:
		"envs/neoloopfinder.yml"
	threads: 2
	shell:
		"""
		assemble-complexSVs -H {input.cnvcool} --output {params.dir} -B {input.sv_bp} \
		--nproc {threads} --logFile {log}
		"""

rule neoloop_caller:
	input:
		cnvcool="results/{sample}/cool_format/{sample}_CNVnorm_%s_%s.cool"%(config['REFERENCE'], config['RESOLUTION']),
		sv="results/{sample}/neoloopfinder/complexSVs/{sample}.assemblies.txt"
	output:
		"results/{sample}/neoloopfinder/neoLoops"
	log:
		"results/{sample}/neoloopfinder/logs/neoloop.log"
	conda:
		"envs/neoloopfinder.yml"
	threads: 2
	shell:
		"""
		neoloop-caller  -H {input.cnvcool} --output {output} --assembly {input.sv} \
		--nproc {threads} --logFile {log}
		"""

rule neotad_caller:
	input:
		cnvcool="results/{sample}/cool_format/{sample}_CNVnorm_%s_%s.cool"%(config['REFERENCE'], config['RESOLUTION']),
		sv="results/{sample}/neoloopfinder/complexSVs/{sample}.assemblies.txt"
	output:
		"results/{sample}/neoloopfinder/neoTads"
	log:
		"results/{sample}/neoloopfinder/logs/neotad.log"
	conda:
		"envs/neoloopfinder.yml"
	threads: 4
	shell:
		"""
		neotad-caller  -H {input.cnvcool} --output {output} --assembly {input.sv} \
		--nproc {threads} --logFile {log}
		"""

rule neo_plots:
	input:
		cnvcool="results/{sample}/cool_format/{sample}_CNVnorm_%s_%s.cool"%(config['REFERENCE'], config['RESOLUTION']),
		sv="results/{sample}/neoloopfinder/complexSVs/{sample}.assemblies.txt",
		tads="results/{sample}/neoloopfinder/neoTads",
		loops="results/{sample}/neoloopfinder/neoLoops"
	output:
		directory("results/{sample}/neoloopfinder/plots")
	conda:
		"envs/neoloopfinder.yml"
	params:
		scratch=config['SCRATCH']
	shell:
		"""
		mkdir -p {params.scratch}/{wildcards.sample}/plots
		dir=$PWD
		cp scripts/plotNeoLoops.py {params.scratch}/{wildcards.sample}/plots
		cd {params.scratch}/{wildcards.sample}/plots
		python plotNeoLoops.py $dir/{input.cnvcool} $dir/{input.loops} $dir/{input.sv} {wildcards.sample}
		mkdir -p $dir/{output}/neoLOOP/
		for f in *.pdf ; do mv $f $dir/{output}/neoLOOP/ ; done
		python plotNeoLoops.py $dir/{input.cnvcool} $dir/{input.tads} $dir/{input.sv} {wildcards.sample}
		mkdir -p $dir/{output}/neoTAD/
		for f in *.pdf ; do mv $f $dir/{output}/neoTAD/ ; done
		cd $dir
		"""

rule assemblies2bedpe:
	input:
		"results/{sample}/neoloopfinder/complexSVs/{sample}.assemblies.txt"
	output:
		TRA="results/{sample}/SVs/{sample}_translocations.bedpe",
		INV="results/{sample}/SVs/{sample}_inversions.bedpe",
		DEL="results/{sample}/SVs/{sample}_deletions.bedpe",
		DUP="results/{sample}/SVs/{sample}_duplications.bedpe"
	shell:
		"""
		prefix="chr"
		cat {input}| awk 'BEGIN {OFS="\t"}  {print $2,$(NF-1),$(NF)}' | sed 's/,/\t/g' | grep "translocation" | sort -k1,1 -k2,2n | awk -v prefix=$prefix 'BEGIN {OFS="\t"} {id++} {print prefix$2,$3,$9,prefix$5,$6,$11,"TRA",".",$4,$7,".","."}'  > {output.TRA}
		cat {input}| awk 'BEGIN {OFS="\t"}  {print $2,$(NF-1),$(NF)}' | sed 's/,/\t/g' | grep "inversion" | sort -k1,1 -k2,2n | awk -v prefix=$prefix 'BEGIN {OFS="\t"} {id++} {print prefix$2,$3,$9,prefix$5,$6,$11,"INV",".",$4,$7,".","."}'  > {output.INV}
		cat {input}| awk 'BEGIN {OFS="\t"}  {print $2,$(NF-1),$(NF)}' | sed 's/,/\t/g' | grep "deletion" | sort -k1,1 -k2,2n | awk -v prefix=$prefix 'BEGIN {OFS="\t"} {id++} {print prefix$2,$3,$9,prefix$5,$6,$11,"DEL",".",$4,$7,".","."}'  > {output.DEL}
		cat {input}| awk 'BEGIN {OFS="\t"}  {print $2,$(NF-1),$(NF)}' | sed 's/,/\t/g' | grep "duplication" | sort -k1,1 -k2,2n | awk -v prefix=$prefix 'BEGIN {OFS="\t"} {id++} {print prefix$2,$3,$9,prefix$5,$6,$11,"DUP",".",$4,$7,".","."}'  > {output.IDUP}
		cat {output} | sort -k1,1 -k2,2n > results/{wildcards.sample}/SVs/all.bedpe
		"""

rule merge_PE:
	input:
		R1=lambda wildcards: sorted(glob("%s/%s/bowtie_results/*R1*.bam"%(hicpro_dir, wildcards.sample))),
		R2=lambda wildcards: sorted(glob("%s/%s/bowtie_results/*R2*.bam"%(hicpro_dir, wildcards.sample)))
	output:
		dynamic("%s/{sample}/{n}.bwt2merged.withSingles.mapq30.bam"%config['SCRATCH'])
	params:
		scratch=config['SCRATCH']
	conda: "envs/hicpro.yml"
	shell:
		"""
		mkdir -p {params.scratch}/{wildcards.sample}/
		
		IFS=' '
		read -ra bw_strands <<< "{input.R2}"
		read -ra fw_strands <<< "{input.R1}"
		
		for i in "${{!fw_strands[@]}}"; do
			echo "R1: ${{fw_strands[$i]}}"
			echo "R2: ${{bw_strands[$i]}}"

			python2 HiCnv-master/scripts/mergeSAM-singletons.py -f "${{fw_strands[$i]}}" -r "${{bw_strands[$i]}}"  -o {params.scratch}/{wildcards.sample}/$i.bwt2merged.withSingles.mapq30.bam -v --single -q 30
		done
		
		"""


rule merge_lanes:
	input: 
		dynamic("%s/{sample}/{n}.bwt2merged.withSingles.mapq30.bam"%config['SCRATCH'])
	output:
		merged="%s/{sample}/{sample}.bwt2merged.merged.withSingles.mapq30.bam"%config['SCRATCH'],
		sorted="results/{sample}/input/{sample}.bwt2merged.sorted.withSingles.mapq30.bam"
	params:
		scratch=config['SCRATCH'],
		hicprodir=hicpro_dir
	threads: 4
	shell:
		"""
		samtools merge -@ {threads} {output.merged} {input}
		samtools sort -@ {threads} -m 2000M -o {output.sorted} {output.merged}
		samtools index {output.sorted}
		"""

rule hic_breakfinder:
	input:
		"results/{sample}/input/{sample}.bwt2merged.sorted.withSingles.mapq30.bam"
	output:
		"results/{sample}/hic_breakfinder/{sample}.breaks.txt"
	params:
		breakfinder=config['HIC_BREAKFINDER_IMG']
	shell:
		"""
		mkdir -p results/{wildcards.sample}/hic_breakfinder/
		singularity exec {params.breakfinder} hic_breakfinder --bam-file {input} \
			--exp-file-inter resources/inter_expect_1Mb.hg19.txt \
			--exp-file-intra resources/intra_expect_100kb.hg19.txt \
			--name results/{wildcards.sample}/hic_breakfinder/{wildcards.sample}
		"""
