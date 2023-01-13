

import os 
from glob import glob

configfile: "peakachu.yml"

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

result_dir=config['OUTDIR']

SAMPLES=config['SAMPLES_ctrl'] + config['SAMPLES_exp1']

rule peakachu_v1:
	input:
		expand("%s/{sample}/peakachu_v2/loops/%s/{sample}_%s.loops"%(result_dir, config['PEAKACHU_MODEL'], config['MIN_CONFIDENCE']), sample=SAMPLES)

rule peakachu_v2:
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
