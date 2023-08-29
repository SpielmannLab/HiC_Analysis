

import os 
from glob import glob

configfile: "peakachu.yml"

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

result_dir=config['OUTDIR']

SAMPLES = config['SAMPLES']
SAMPLES_diff = config['SAMPLES_ctrl'] + config['SAMPLES_exp1']

rule peakachu_v2:
	input:
		expand("%s/{sample}/peakachu_v2/loops/%s/{sample}_%s.loops"%(result_dir, config['PEAKACHU_MODEL'], config['MIN_CONFIDENCE']), sample=SAMPLES)

rule peakachu_v1:
	input: 
		expand("%s/{sample}/peakachu_v1/loops/merged/{sample}_%s.loops"%(result_dir, config['MIN_CONFIDENCE']), sample=SAMPLES)

### note peakachu version 1 throws errors because scikit-learn=0.20.2  and not scikit-learn>0.21 ?
### but version 0.20.2 is requested on github page
rule peakachu1:
	input:
		hic="%s/{sample}/hic_format/{sample}_%s.hic"%(config['HICPRO_INDIR'], config['REFERENCE'])
	output:
		ctcf="%s/{sample}/peakachu_v1/loops/ctcf/{sample}.loops.txt"%(result_dir),
		histone="%s/{sample}/peakachu_v1/loops/h3k27ac/{sample}.loops.txt"%(result_dir)
	conda:
		"envs/fanc-peakachu.yml"
	params:
		scratch=config['SCRATCH'],
		model=config['PEAKACHU_MODEL'],
		dir=result_dir,
		res=config['RESOLUTION']
	shell:
		"""
		mkdir -p {params.scratch}/{wildcards.sample}/peakachu_v1/{params.res}/scores/ctcf/
		peakachu score_genome -p {input}  -r {params.res} -O {params.scratch}/{wildcards.sample}/peakachu_v1/{params.res}/scores/ctcf/  -m {params.model}.ctcf.pkl  || true
		
		touch {output.ctcf} 
		cat {params.scratch}/{wildcards.sample}/peakachu_v1/*/scores/ctcf/* > {output.ctcf} 
		
		mkdir -p {params.scratch}/{wildcards.sample}/peakachu_v1/{params.res}/scores/h3k27ac/
		peakachu score_genome -p {input} -r {params.res} -O {params.scratch}/{wildcards.sample}/peakachu_v1/{params.res}/scores/h3k27ac/  -m {params.model}.h3k27ac.pkl  || true
		
		touch {output.histone}
		cat {params.scratch}/{wildcards.sample}/peakachu_v1/*/scores/h3k27ac/* > {output.histone} 
		
		"""
		
rule peakachu2:
	input:
		hic="%s/{sample}/hic_format/{sample}_%s.hic"%(config['HICPRO_INDIR'], config['REFERENCE'])
	output:
		scores="%s/{sample}/peakachu_v2/loops/%s/{sample}.loops.txt"%(result_dir, config['PEAKACHU_MODEL'])
	conda:
		"envs/peakachu2.yml"
	params:
		scratch=config['SCRATCH'],
		model=config['PEAKACHU_MODEL'],
		dir=result_dir,
		res=config['RESOLUTION']
	shell:
		"""
		mkdir -p {params.scratch}/{wildcards.sample}/peakachu/{params.res}/scores/
		peakachu score_genome --minimum-prob 0 -p {input}  -r {params.res} -O {output}  -m {params.model}
		"""

### bedpe format for 2D annotation in JuiceBox
rule pool_peakachu2:
	input:
		scores="%s/{sample}/peakachu_v2/loops/%s/{sample}.loops.txt"%(result_dir, config['PEAKACHU_MODEL'])
	output:
		loops="%s/{sample}/peakachu_v2/loops/%s/{sample}_%s.loops"%(result_dir, config['PEAKACHU_MODEL'], config['MIN_CONFIDENCE']),
		bedpe="%s/{sample}/peakachu_v2/loops/%s/{sample}_%s.bedpe"%(result_dir, config['PEAKACHU_MODEL'], config['MIN_CONFIDENCE'])
	conda: "envs/peakachu2.yml"
	params:
		th=config['MIN_CONFIDENCE'],
		scratch=config['SCRATCH'],
		res=config['RESOLUTION']
	shell:
		"""
		peakachu pool -i {input} -t {params.th}  -o {output.loops} -r {params.res}
		cat {output.loops} | awk ' BEGIN {{OFS="\\t"}} {{print $1,$2,$5,$1,$2,$5,".",".",".",".","0,255,255"}}' > {output.bedpe}
		"""

rule pool_peakachu1:
	input:
		ctcf="%s/{sample}/peakachu_v1/loops/ctcf/{sample}.loops.txt"%result_dir,
		histon="%s/{sample}/peakachu_v1/loops/h3k27ac/{sample}.loops.txt"%result_dir
	output:
		ctcf="%s/{sample}/peakachu_v1/loops/ctcf/{sample}_%s.txt"%(result_dir, config['MIN_CONFIDENCE']),
		histon="%s/{sample}/peakachu_v1/loops/h3k27ac/{sample}_%s.txt"%(result_dir, config['MIN_CONFIDENCE'])
	conda: "envs/fanc-peakachu.yml"
	params:
		th=config['MIN_CONFIDENCE']
	threads: 2
	shell:
		"""
		peakachu pool -i {input.ctcf} -t {params.th} > {output.ctcf} 
		peakachu pool -i {input.histon} -t {params.th} > {output.histon} 
		"""

### it is recommended to check for the best thresholds each and merge them 
### for simplicity here the same threshold is (MIN_CONFIDENCE) is used for both models
rule merge_peakachu1:
	input:
		ctcf="%s/{sample}/peakachu_v1/loops/ctcf/{sample}_%s.txt"%(result_dir, config['MIN_CONFIDENCE']),
		histon="%s/{sample}/peakachu_v1/loops/h3k27ac/{sample}_%s.txt"%(result_dir, config['MIN_CONFIDENCE'])
	output:
		merged="%s/{sample}/peakachu_v1/loops/merged/{sample}_%s.loops"%(result_dir, config['MIN_CONFIDENCE']),
		bedpe="%s/{sample}/peakachu_v1/loops/merged/{sample}_%s.bedpe"%(result_dir, config['MIN_CONFIDENCE'])
	conda: "envs/fanc-peakachu.yml"
	params:
		scratch=config['SCRATCH']
	shell:
		"""
		id=$RANDOM
		bedtools pairtopair -is -slop 25000 -type notboth -a {input.histon} -b {input.ctcf} > {params.scratch}/$id
		cat {input.ctcf} {params.scratch}/$id > {output.merged}
		
		awk ' BEGIN {{OFS="\\t"}} {{print $1,$2,$5,$1,$2,$5,".",".",".",".","0,255,255"}}' {output.merged} > {output.bedpe}
		
		rm  {params.scratch}/$id
		"""





#### select thresholds for peakachu pool (see rule filter_peakachu_scores) and 
#### collect the best results for each control sample in directory peakachudiff_input_controls
#### collect the best results for each sample in directory peakachudiff_input

rule merge_control_loops:
	input:
		glob("peakachudiff_input_controls/*.txt") 
	output:
		"%s/comparison/control.loops"%(result_dir)
	conda: "../envs/peakachu2.yml"
	params:
		scratch=config['SCRATCH']
	shell: 
		"""
		IFS=' ' read -r -a files <<< "{input}"
		file1="${{files[0]}}"
		for j in "${{files[@]}}" ; do
			file2=$j
			if [[ $file2  != $file1 ]] ; then
				bedtools pairtopair -is -slop 25000 -type notboth -a $file1 -b $file2 > {params.scratch}/tmp
				cat $file2 {params.scratch}/tmp >> {params.scratch}/tmp.collection
			fi
			done
		cat {params.scratch}/tmp.collection | sort | uniq > {output}
		"""

rule compare2allCtrlloops:
	input:
		ctrl="%s/comparison/control.loops"%(result_dir),
		samples= glob("peakachudiff_input_samples/*.txt") 
	output:
		"%s/PeakachuDiff/merged_control/control.unique.loops"%(result_dir),
		expand("%s/PeakachuDiff/merged_control/{sample}.control.unique.loops"%(result_dir), sample=config['SAMPLES_exp1'])
	conda: "../envs/peakachu2.yml"
	params:
		scratch=config['SCRATCH'],
		dir=result_dir
	shell: 
		"""
		cp diffPeakachu/pair-probs.py {params.scratch}/pair-probs.py
			
		sed  '/D = /i# peakachu2 workaround\\n    fil1=cell\\n    fil2=cell'  work/pair-probs.py {params.scratch}/pair-probs.py
			
		for file in {input.samples} ; do 
			IFS='.' read -ra ADDR <<< $(basename $file )
			name="${{ADDR[0]}}"
			python {params.scratch}/pair-probs.py $file {input.ctrl} {params.scratch}/${{name}}-control.merged.loops
			python diffPeakachu/diffPeakachu.py $file {input.ctrl}  {params.scratch}/${{name}}-control.merged.loops
			cp {params.scratch}/$(basename $file .txt)-control.unique.loops {params.dir}/PeakachuDiff/${{name}}-control.unique.loops
		done

		cat $(ls {params.scratch}/control-*.unique.loops) > {params.dir}/PeakachuDiff/merged_control/control.unique.loops
		"""
rule compare2eachCtrlloop:
	input:
		ctrl=glob("peakachudiff_input_controls/*.txt"), 
		samples=glob("peakachudiff_input_samples/*.txt") 
	output:
		expand("%s/PeakachuDiff/{sample}.{ctrl}.unique.loops"%(result_dir), sample=config['SAMPLES_exp1'],ctrl=config['SAMPLES_ctrl']),
		expand("%s/PeakachuDiff/{ctrl}.{sample}.unique.loops"%(result_dir), sample=config['SAMPLES_exp1'],ctrl=config['SAMPLES_ctrl'])
	conda: "../envs/peakachu2.yml"
	params:
		scratch=config['SCRATCH'],
		dir=result_dir
	shell: 
		"""
		for ctrl in {input.ctrl} ; do 
			for file in {input.samples} ; do 
				IFS='.' read -ra ADDR <<< $(basename $file )
				name="${{ADDR[0]}}"
				IFS='.' read -ra ADDR <<< $(basename $ctrl )
				c_name="${{ADDR[0]}}"
				echo "sample is "$name
				echo "control is "$c_name
				python diffPeakachu/pair-probs.py $file $ctrl {params.scratch}/${{name}}-${{c_name}}.merged.loops
				python diffPeakachu/diffPeakachu.py $file $ctrl {params.scratch}/${{name}}-${{c_name}}.merged.loops
				
				if [[ -f {params.scratch}/$(basename $file .txt)-$(basename $ctrl .txt).unique.loops ]] ; then
					cp {params.scratch}/$(basename $file .txt)-$(basename $ctrl .txt).unique.loops {params.dir}/PeakachuDiff/${{name}}-${{c_name}}.unique.loops 
				fi
				if [[ -f {params.scratch}/$(basename $ctrl .txt)-$(basename $file .txt).unique.loops ]] ; then
					cp {params.scratch}/$(basename $ctrl .txt)-$(basename $file .txt).unique.loops {params.dir}/PeakachuDiff/${{c_name}}-${{name}}.unique.loops 
				fi				
			done
		done

		"""


rule annotateLoops_Ctrl:
	input: 
		"%s/PeakachuDiff/merged_control/control.unique.loops"%(result_dir)
	output:
		"%s/PeakachuDiff/merged_control/annotated/control.unique.loops.annotated"%(result_dir)
	params:
		genes=config['GENELIST']
	shell:
		"""
		T=$(printf '\\t')
		echo "gene_chr${{T}}gene_start${{T}}gene_end${{T}}genesymbol${{T}}loop_start_chr${{T}}loop_start_pos1${{T}}loops_start_pos2${{T}}loop_end_chr${{T}}loop_end_pos1${{T}}loop_end_pos2" > {output}
		bedtools intersect -loj -a {params.genes} -b {input}| awk '!/-1/' >> {output}
		"""

rule annotateLoops_mergedCtrl:
	input: 
		"%s/PeakachuDiff/merged_control/{sample}.control.unique.loops"%(result_dir),
		
	output:
		"%s/PeakachuDiff/merged_control/annotated/{sample}.control.unique.loops.annotated"%(result_dir)
	params:
		genes=config['GENELIST']
	shell:
		"""
		T=$(printf '\\t')
		echo "gene_chr${{T}}gene_start${{T}}gene_end${{T}}genesymbol${{T}}loop_start_chr${{T}}loop_start_pos1${{T}}loops_start_pos2${{T}}loop_end_chr${{T}}loop_end_pos1${{T}}loop_end_pos2" > {output}
		bedtools intersect -loj -a {params.genes} -b {input}| awk '!/-1/' >> {output}
		"""

rule annotateLoops_eachCtrl:
	input: 
		"%s/PeakachuDiff/{sample}.{ctrl}.unique.loops"%(result_dir)
	output:
		"%s/PeakachuDiff/annotated/{sample}.{ctrl}.unique.loops.annotated"%(result_dir)
	params:
		genes=config['GENELIST']
	shell:
		"""
		T=$(printf '\\t')
		echo "gene_chr${{T}}gene_start${{T}}gene_end${{T}}genesymbol${{T}}loop_start_chr${{T}}loop_start_pos1${{T}}loops_start_pos2${{T}}loop_end_chr${{T}}loop_end_pos1${{T}}loop_end_pos2" > {output}
		bedtools intersect -loj -a {params.genes} -b {input}| awk '!/-1/' >> {output}
		"""


rule peakachuDiff_mergedCtrl:
	input:
		ctrl="%s/PeakachuDiff/merged_control/annotated/control.unique.loops.annotated"%(result_dir),
		samples=expand("%s/PeakachuDiff/merged_control/annotated/{sample}.control.unique.loops.annotated"%(result_dir), sample=config['SAMPLES_exp1'])
	output:
		ctrl="%s/PeakachuDiff/merged_control/uniqueGenes/control.unique.genes.txt"%(result_dir)
	params:
		scratch=config['SCRATCH'],
		dir=result_dir
	shell:
		"""
		cat {input.ctrl} | awk '{{ print $4 }}' | sort | uniq > {output.ctrl}
		for sample in {input.samples} ; do
			name=$(basename $sample .control.unique.loops.annotated)
			output={params.dir}/PeakachuDiff/merged_control/uniqueGenes/$name.unique.genes.txt
			echo "#$name" > $output
			cat $sample | awk '{{ print $4 }}' | sort | uniq >> $output
		done
		
		"""
		
rule  peakachuDiff_eachCtrl:
	input:
		ctrl=expand("%s/PeakachuDiff/annotated/{ctrl}.{sample}.unique.loops.annotated"%(result_dir), sample=config['SAMPLES_exp1'], ctrl=config['SAMPLES_ctrl']),
		samples=expand("%s/PeakachuDiff/annotated/{sample}.{ctrl}.unique.loops.annotated"%(result_dir), sample=config['SAMPLES_exp1'], ctrl=config['SAMPLES_ctrl'])
	output:
		directory("%s/PeakachuDiff/uniqueGenes"%(result_dir))
	params:
		scratch=config['SCRATCH'],
		dir=result_dir
	shell:
		"""
		for file in {input} ; do
			name=$(basename $file .unique.loops.annotated)
			output={params.dir}/PeakachuDiff/uniqueGenes/$name.unique.genes.txt
			echo "#$name" > $output
			cat $file | awk '{{ print $4 }}' | sort | uniq >> $output
		done
		
		"""



