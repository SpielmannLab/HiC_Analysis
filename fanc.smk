import os 
from glob import glob

configfile: "fanc.yml"

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

result_dir=config['OUTDIR']

SAMPLES=config['SAMPLES_ctrl'] + config['SAMPLES_exp1']

rule target_regions:
	output:
		region="%s/cohort/target_regions_r%s.txt"%(config['SCRATCH'],config['RANGE'])
	params: 
		genes=config['TARGET_GENES'],
		list=config['GENELIST'],
		res=config['RESOLUTION'],
		range=config['RANGE']
	shell:
		"""
		echo -n > {output.fancplot}
		echo -n > {output.LDNet}
		echo -n > {output.bed}
		
		for gene in {params.genes} ; do
			chr=$(cat {params.list} | grep -w $gene | head -1 | awk 'BEGIN {{FS="\\t"}} {{print $1}}') || true
			if [[ $chr == "" ]] ; then
				echo $gene" not found in "{params.list}
				break
			fi
			genestart=$(cat {params.list} | grep -w $gene | head -1 | awk 'BEGIN {{ FS="\\t"}} {{print $2}}')
			win_start=$(($genestart - {params.range} * {params.res}))
			win_start_kb=$((($win_start / {params.res}) *{params.res} ))
			win_end=$(($genestart + {params.range} * {params.res}))
			win_end_kb=$((($win_end / {params.res}) *{params.res} ))
			
			echo -e  $chr":"$win_start_kb"-"$win_end_kb >> {output.region}
		done			
		"""

rule all_loops:
	input:
		expand(directory("results/{sample}/fanc/loops/plots"), sample=SAMPLES)

rule loops:
	input:
		hic="%s/{sample}/hic_format/{sample}_%s.hic"%(config['HICPRO_INDIR'], config['REFERENCE'])
	output:
		smallLoops="results/{sample}/fanc/loops/{sample}.10kb.loops",
		bigLoops="results/{sample}/fanc/loops/{sample}.25kb.loops",
		allLoops="results/{sample}/fanc/loops/{sample}.all.loops"
	conda: "envs/fanc-hic.yml"
	shell:
		"""
		fanc loops {input.hic}@10000 {output.smallLoops}
		fanc loops {input.hic}@25000 {output.bigLoops}
		cat {output.smallLoops}  {output.bigLoops} >  {output.allLoops}
		"""	

rule loops_plots:
	input:
		smallLoops="results/{sample}/fanc/loops/{sample}.10kb.loops",
		bigLoops="results/{sample}/fanc/loops/{sample}.25kb.loops",
		hic="%s/{sample}/hic_format/{sample}_%s.hic"%(config['HICPRO_INDIR'], config['REFERENCE']),
		region="%s/cohort/target_regions_r%s.txt"%(result_dir, config['RANGE'])
	output:
		directory("results/{sample}/fanc/loops/plots")
	conda: "envs/fanc-hic.yml"
	params: genes=config['TARGET_GENES']
	shell:
		"""
		i=0
		IFS=' ' read -r -a genes <<< "{params.genes}"
		
		while read region ; do
			fancplot  --pdf-text-as-font  -o {output}/"${{genes[$i]}}".pdf "$region" --plot square  {input.smallLoops}
			fancplot  --pdf-text-as-font  -o {output}/"${{genes[$i]}}".pdf "$region" --plot square  {input.bigLoops}
			i=$(($i+1))
		done < {input.region}
		"""	

### fanc insulation does their own normalization so we're forwarding no norm ( @NONE )
rule insulation:
	input:
		"%s/{sample}/hic_format/{sample}_%s.hic"%(config['HICPRO_INDIR'], config['REFERENCE'])
	output:
		all="%s/{sample}/fanc/insulation/{sample}_10kb.insulation"%(result_dir)
	conda: "envs/fanc-hic.yml"
	shell:
		"""
		fanc insulation {input}@10000@NONE {output.all} -w 50000 100000 250000 --impute
		fanc insulation {output.all} -o bed -w 50kb 100kb 250kb
		"""

rule insulation2bed:
	input:
		"%s/{sample}/fanc/insulation/{sample}_10kb.insulation"%(result_dir)
	output:
 		s="%s/{sample}/fanc/insulation/{sample}_10kb.insulation_50kb.bed"%(result_dir),
 		m="%s/{sample}/fanc/insulation/{sample}_10kb.insulation_100kb.bed"%(result_dir),
 		l="%s/{sample}/fanc/insulation/{sample}_10kb.insulation_250kb.bed"%(result_dir)
	conda: "envs/fanc-hic.yml"
	shell:
		"""
		fanc insulation {input} -o bed -w 50kb 250kb 100kb 
		"""

rule boundaries:
	input:
		"%s/{sample}/fanc/insulation/{sample}_10kb.insulation"%(result_dir)
	output:
		"%s/{sample}/fanc/boundaries/{sample}_10kb_25kb.bed"%(result_dir),
		"%s/{sample}/fanc/boundaries/{sample}_10kb_50kb.bed"%(result_dir),
		"%s/{sample}/fanc/boundaries/{sample}_10kb_100kb.bed"%(result_dir),
		"%s/{sample}/fanc/boundaries/{sample}_10kb_250kb.bed"%(result_dir)
	conda: "../envs/fanc.yml"
	params:
		scratch=config['SCRATCH'],
		dir=result_dir
	shell:
		"""
		mkdir -p {params.dir}/{wildcards.sample}/fanc/boundaries/
		fanc boundaries {input} {params.dir}/{wildcards.sample}/fanc/boundaries/{wildcards.sample}_10kb
		"""


### change the max_saturation if necessary 
### some .hic files behave a bit off ... 
rule insulation_plots:
	input:
		hic=expand("%s/{sample}/hic_format/{sample}_%s.hic"%(config['HICPRO_INDIR'], config['REFERENCE']), sample=SAMPLES),
		scores=expand("%s/{sample}/fanc/insulation/{sample}_10kb.insulation_100kb.bed"%(result_dir), sample=SAMPLES),
		region="%s/cohort/target_regions_r%s.txt"%(result_dir,config['RANGE'])
	output:
		expand("%s/cohort/plots/fanc/insulation/{gene}.svg"%(result_dir), gene=config['TARGET_GENES'])
	conda: "envs/fanc-hic.yml"
	params:
		dir=result_dir,
		samples= SAMPLES,
		ref=config['REFERENCE'],
		scratch=config['SCRATCH'],
		prefix=config['PREFIX'],
		genes=config['TARGET_GENES']
	shell:
		"""	
		mapString=""
		for hic in {input.hic} ; do
			name=$(basename $hic {params.ref}.hic)
			max_saturation=20
			mapString=$mapString" -p triangular -vmax $max_saturation "$hic@10000
		done

		mkdir -p {params.scratch}/insulation_scores/
		for f in {input.scores} ; do
			name=$(basename $f)
 			cat $f | sed 's/chr//g' | awk 'BEGIN {{OFS="\\t"}}  $1="{params.prefix}"$1' > {params.scratch}/insulation_scores/$name 
		done
		echo $mapString
		i=0
		IFS=' ' read -r -a genes <<< "{params.genes}"
		while read region ; do
			fancplot  --pdf-text-as-font  -o {params.dir}/cohort/plots/fanc/insulation/"${{genes[$i]}}".svg "$region"  \
 				$mapString \
				-p line $(ls {params.scratch}/insulation_scores/*.bed) -l {params.samples}
				
			fancplot  --pdf-text-as-font  -o {params.dir}/cohort/plots/fanc/insulation/"${{genes[$i]}}".png "$region"  \
 				$mapString \
				-p line $(ls {params.scratch}/insulation_scores/*.bed) -l {params.samples}
			i=$(($i+1))
		done < {input.region}
		"""	

### haven't tried that but it should work like this

rule plot_APA:
	input:
		hic="%s/{sample}/hic_format/{sample}_%s.hic"%(config['HICPRO_INDIR'], config['REFERENCE']),
		region="%s/cohort/target_regions_r%s.txt"%(result_dir, config['RANGE'])
	output:
		expand("%s/cohort/plots/fanc/APA/{gene}.png"%(result_dir), gene=config['TARGET_GENES'])	
	conda: "envs/fanc-hic.yml"
	params:
		dir=result_dir,
	shell:
		"""
		i=0
		IFS=' ' read -r -a genes <<< "{params.genes}"
		
		while read region ; do
			 fanc aggregate -p {params.dir}/cohort/plots/fanc/APA/"${{genes[$i]}}".png {input.hic} $region
		done < {input.region}
		"""

rule matrix_APA:
	input:
		hic="%s/{sample}/hic_format/{sample}_%s.hic"%(config['HICPRO_INDIR'], config['REFERENCE']),
		region="%s/cohort/target_regions_r%s.txt"%(result_dir, config['RANGE'])
	output:
		expand("%s/cohort/plots/fanc/APA/{gene}.matrix"%(result_dir), gene=config['TARGET_GENES'])	
	conda: "envs/fanc-hic.yml"
	params:
		dir=result_dir,
	shell:
		"""
		i=0
		IFS=' ' read -r -a genes <<< "{params.genes}"
		
		while read region ; do
			 fanc aggregate -m {params.dir}/cohort/plots/fanc/APA/"${{genes[$i]}}".matrix {input.hic} $region
		done < {input.region}
		"""

