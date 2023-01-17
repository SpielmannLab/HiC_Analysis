import os 
from glob import glob
from random import random
from itertools import permutations as p
configfile: "loopcalling.yml"

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

result_dir=config['OUTDIR']

SAMPLES=config['SAMPLES_ctrl'] + config['SAMPLES_exp1']

###ldl net....
rule target_regions:
	output:
		LDNet="%s/cohort/LDNet_target_regions_r%s.txt"%(result_dir,config['RANGE'])
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
			
			chr=$(echo $chr | sed 's/chr//g')
			echo -e $chr'\\t'$win_start_kb'\\t'$win_end_kb >> {output.LDNet}
		done			
		"""

rule merge:
	input:
		expand("%s/{sample}/{sample}_TAD_calls.bedpe"%(result_dir), sample=SAMPLES)

rule arrowhead_all:
	input:
		expand("%s/{sample}/arrowhead/{sample}_arrowhead.bedpe"%result_dir, sample=SAMPLES)

rule arrowhead:
	input:
		"%s/{sample}/hic_format/{sample}_%s.hic"%(config['HICPRO_INDIR'], config['REFERENCE'])
	output:
		"%s/{sample}/arrowhead/{sample}_arrowhead.bedpe"%result_dir
	params:
		juicer=config['JUICER'],
		dir=result_dir
	threads: 4
	resources:
		mem_mb=13*1000
	shell:
		"""
		export _JAVA_OPTIONS="-Xmx12g -Xms5g"
		java -jar {params.juicer} arrowhead -r 10000 -k KR --threads {threads} --ignore-sparsity {input} {params.dir}/{wildcards.sample}/arrowhead
		java -jar {params.juicer} arrowhead -r 25000 -k KR --threads {threads} --ignore-sparsity {input} {params.dir}/{wildcards.sample}/arrowhead
		java -jar {params.juicer} arrowhead -r 50000 -k KR --threads {threads} --ignore-sparsity {input} {params.dir}/{wildcards.sample}/arrowhead
		java -jar {params.juicer} arrowhead -r 100000 -k KR --threads {threads} --ignore-sparsity {input} {params.dir}/{wildcards.sample}/arrowhead
		cat $( ls {params.dir}/{wildcards.sample}/arrowhead/*0_blocks.bedpe) > {params.dir}/{wildcards.sample}/arrowhead/all_blocks.bedpe
		awk 'BEGIN {{OFS="\\t"}} {{print $1,$2,$3,$4,$5,$6,".",".",".",".","255,255,0"}}' {params.dir}/{wildcards.sample}/arrowhead/all_blocks.bedpe > {output}
		"""
		
rule hiccups_all:
	input:
		expand(directory("%s/{sample}/hiccups/"%result_dir), sample=SAMPLES)
		
rule hiccups:
	input:
		"%s/{sample}/hic_format/{sample}_%s.hic"%(config['HICPRO_INDIR'], config['REFERENCE'])
	output:
		directory("%s/{sample}/hiccups/"%result_dir)
	params:
		dir=result_dir,
		juicer=config['JUICER'],
		scratch=config['SCRATCH']
	threads: 1
	shell:
		"""
		mkdir -p {output}
		job=$RANDOM
		echo "#!/bin/bash" > {params.scratch}/$job.sh
		echo "source /etc/profile.d/modules.sh; module load nvidia-cuda/11.1-native; java -jar {params.juicer} hiccups --threads 4 -r 10000,25000 -k KR --ignore-sparsity $PWD/{input} $PWD/{params.dir}/{wildcards.sample}/hiccups/" >> {params.scratch}/$job.sh
		sbatch --gres=gpu:1 --job-name=hiccups-$job {params.scratch}/$job.sh
		"""


rule  HiC_LDNet_all:
	input:
		expand(directory("%s/{sample}/HiC_LDNet/regions"%(result_dir)), sample=SAMPLES)

rule  HiC_LDNet:
	input:
		list="%s/cohort/LDNet_target_regions_r%s.txt"%(result_dir, config['RANGE']),
		hic="%s/{sample}/hic_format/{sample}_%s.hic"%(config['HICPRO_INDIR'], config['REFERENCE'])
	output:
		directory("%s/{sample}/HiC_LDNet/regions"%(result_dir))
	conda: "../envs/HiC-LDNet.yml"
	params:
		scratch=config['SCRATCH'],
		model=config['LDNet_MODEL'],
		cohort=result_dir,
		prefix=config['PREFIX']
	shell:
		"""
		echo '#!/bin/bash' > {params.scratch}/LDNet-script.sh
		echo "PATH=$PATH:/home/schultz/miniconda3/bin" >> {params.scratch}/LDNet-script.sh
		echo "source activate /data/humangen_external/schultz/projects/HiC/pipeline/Loop_Calling/.snakemake/conda/966a9654" >> {params.scratch}/LDNet-script.sh
		echo "source /etc/profile.d/modules.sh" >> {params.scratch}/LDNet-script.sh
		echo "module load nvidia-cuda/11.1-native" >> {params.scratch}/LDNet-script.sh
		mkdir -p {output}
		while read line; do
  			id=$RANDOM
  			cp {params.scratch}/LDNet-script.sh {params.scratch}/LDNet-script_$id.sh
  			while IFS=$'\\t' read -r -a tmp ; do
  				chr="${{tmp[0]}}"
  				start="${{tmp[1]}}"
    			end="${{tmp[2]}}"
				echo "mkdir -p {output}/${{chr}}_${{start}}_${{end}}" >> {params.scratch}/LDNet-script_$id.sh
    			echo "python HiC-LDNet/detect_loops.py --hic_file={input.hic} --output_path={output}/${{chr}}_${{start}}_${{end}}  --model_path="{params.model}" --start=$start --end=$end  --chr=$chr --prefix={params.prefix} " >> {params.scratch}/LDNet-script_$id.sh
			done <<< $line
			sbatch --gres=gpu:1 --job-name=LDNet-$id --time=0-0:20:0 {params.scratch}/LDNet-script_$id.sh
		done <{input.list}
		"""
			
rule hiccups2bedpe:
	input:
		"%s/{sample}/hiccups/merged_loops.bedpe"%result_dir
	output:
		"%s/{sample}/hiccups/{sample}_hiccups.bedpe"%result_dir
	params:
		scratch=config['SCRATCH']
	shell:
		"""
		awk ' BEGIN {{OFS="\\t"}} {{print $1,$2,$3,$4,$5,$6,".",".",".",".","255,0,255"}}' {input} > {output}
		"""

rule bedpe2callset:
	input:
		hiccups="%s/{sample}/hiccups/{sample}_hiccups.bedpe"%result_dir,
		arrowhead="%s/{sample}/arrowhead/{sample}_arrowhead.bedpe"%result_dir
#		peakachu="%s/{sample}/peakachu/loops/%s/{sample}.loops.filtered%s.bedpe"%(result_dir, config['PEAKACHU_MODEL'], config['MIN_CONFIDENCE'])
	output:
		"%s/{sample}/{sample}_TAD_calls.bedpe"%result_dir
	params:
		scratch=config['SCRATCH']
	shell:
		"""
		cat {input} | grep -iv '#' | sort -k1,1 -k2,2n > {params.scratch}/{wildcards.sample}_TAD_calls.bedpe
		
		sed -i -r 's/chr//g' {params.scratch}/{wildcards.sample}_TAD_calls.bedpe
		
		awk 'BEGIN {{OFS="\t"}}  $1="chr"$1, $4="chr"$4 ' {params.scratch}/{wildcards.sample}_TAD_calls.bedpe > {output}
		"""