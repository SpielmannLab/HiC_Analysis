import os 
from glob import glob

configfile: "Compartment_Analysis.yml"

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

result_dir=config['OUTDIR']

SAMPLES=config['SAMPLES_ctrl'] + config['SAMPLES_exp1'] 

res=["50000", "500000", "100000"] # bin size for compartment analysis


rule genome_overview:
	input:
		expand(directory("%s/{sample}/fanc/compartments/genome/"%result_dir), sample=SAMPLES)

rule target_genes:
	input:
		expand(directory("%s/{sample}/fanc/compartments/targetRegions/r%s/"%(result_dir, config['RANGE'])), sample=SAMPLES)

rule compareToCtrls:
	input:
		expand("%s/cohort/ABcompartments/{sample}/{sample}.log"%result_dir, sample=SAMPLES)

rule fanc_compartments:
	input:
		expand("%s/{sample}/fanc/compartments/matrix/{sample}.{res}.matrix"%result_dir, sample=SAMPLES, res=res),
		expand("%s/{sample}/fanc/compartments/matrix/{sample}.{res}.AB.bed"%result_dir, sample=SAMPLES, res=res),
		expand("%s/{sample}/fanc/compartments/matrix/{sample}.{res}.ev.bed"%result_dir, sample=SAMPLES, res=res)

rule fanc_compartments_n_kb:
	input:
		"%s/{sample}/hic_format/{sample}_%s.hic"%(config['HICPRO_INDIR'], config['REFERENCE'])
	output:
		matrix="%s/{sample}/fanc/compartments/matrix/{sample}.{res}.matrix"%result_dir,
		domain="%s/{sample}/fanc/compartments/matrix/{sample}.{res}.AB.bed"%result_dir,
		ev="%s/{sample}/fanc/compartments/matrix/{sample}.{res}.ev.bed"%result_dir
	params:
		ref=config['REFERENCE_FILE']
	conda: "envs/fanc.yml"
	shell:
		"""
		fanc compartments -d {output.domain} -v {output.ev} -g {params.ref} {input}@{wildcards.res} {output.matrix}
		"""

### get AB bins from control samples
rule get_AB_control:
	input:
		expand("%s/{sample}/fanc/compartments/matrix/{sample}.50000.AB.bed"%result_dir, sample=config['SAMPLES_ctrl'])
	output:
		A=expand("%s/{sample}.50000.A.bed"%config['SCRATCH'], sample=config['SAMPLES_ctrl']),
		B=expand("%s/{sample}.50000.B.bed"%config['SCRATCH'], sample=config['SAMPLES_ctrl'])
	params:
		dir=result_dir,
		scratch=config['SCRATCH']
	shell:
		"""
		for f in {input}; do
			name=$(basename $f .50000.AB.bed)
			cat $f | awk '{{ OFS="\\t"}} {{ if ($4 == "A") {{ print $1,$2,$3 }} }}' > {params.scratch}/${{name}}.50000.A.bed
			cat $f | awk '{{ OFS="\\t"}} {{ if ($4 == "B") {{ print $1,$2,$3 }} }}' > {params.scratch}/${{name}}.50000.B.bed
		done
		"""

### get consensus bins for all control samples. these are our conserved A B bins.
rule consens_A_control:
	input:
		expand("%s/{sample}.50000.A.bed"%config['SCRATCH'], sample=config['SAMPLES_ctrl'])
	output:
		"%s/ctrl_consens_A.bed"%config['SCRATCH']
	params:
		dir=result_dir,
		scratch=config['SCRATCH']
	shell:	
		"""
		IFS=' ' read -r -a files <<< "{input}"
		len=${{#files[@]}}
		
		if [[ $len -gt 1 ]] ; then
			bedtools multiinter -i {input} > {params.scratch}/ctrl_a_consens.bed
			awk -v len=$len 'BEGIN {{OFS="\\t"}} {{ if ($4 == len ) {{print $1,$2,$3}} }}'  {params.scratch}/ctrl_a_consens.bed > {output}
		else
			cp {input} {output} 
		fi		
		"""

rule consens_B_control:
	input:
		expand("%s/{sample}.50000.B.bed"%config['SCRATCH'], sample=config['SAMPLES_ctrl'])
	output:
		"%s/ctrl_consens_B.bed"%config['SCRATCH']
	params:
		dir=result_dir,
		scratch=config['SCRATCH']
	shell:	
		"""
		IFS=' ' read -r -a files <<< "{input}"
		len=${{#files[@]}}
		
		if [[ $len -gt 1 ]] ; then
			bedtools multiinter -i {input} > {params.scratch}/ctrl_b_consens.bed
			awk -v len=$len 'BEGIN {{OFS="\\t"}} {{ if ($4 == len ) {{print $1,$2,$3}} }}'  {params.scratch}/ctrl_b_consens.bed > {output}
		else
			cp {input} {output} 
		fi	
		"""	

rule compare_AB_control:
	input:
		"%s/ctrl_consens_B.bed"%config['SCRATCH'],
		"%s/ctrl_consens_A.bed"%config['SCRATCH']
	output:
		final="%s/cohort/ABcompartments/conversed_compartments.bed"%result_dir
	params:
		dir=result_dir,
		scratch=config['SCRATCH']
	shell:	
		"""
		cat {input} | sort -k1,1 -k2,2n > {params.scratch}/ctrl_consens.bed
		bedtools merge -i {params.scratch}/ctrl_consens.bed > {output.final}
		"""		

rule compareAB_exp1_all:
	input:
		expand("%s/cohort/ABcompartments/{sample}/{sample}.log"%result_dir, sample=config['SAMPLES_exp1'])
		
### compare conserved bins to each experiment sample 	
rule compareAB_exp1:
	input:
		B="%s/ctrl_consens_B.bed"%config['SCRATCH'],
		A="%s/ctrl_consens_A.bed"%config['SCRATCH'],
		sample="%s/{sample}/fanc/compartments/matrix/{sample}.50000.AB.bed"%result_dir
	output:
		"%s/cohort/ABcompartments/{sample}/{sample}.log"%result_dir
	params:
		dir=result_dir,
		scratch=config['SCRATCH'],
		genes=config['GENELIST']
	shell:	
		"""
		name=$(basename {input.sample} .50000.AB.bed)
		mkdir -p {params.dir}/cohort/ABcompartments/$name
		
		cat {input.sample} | awk '{{ OFS="\\t"}} {{ if ($4 == "B") {{ print $1,$2,$3 }} }}' > {params.scratch}/$name.B.bed
		cat {input.sample} | awk '{{ OFS="\\t"}} {{ if ($4 == "A") {{ print $1,$2,$3 }} }}' > {params.scratch}/$name.A.bed
		
		bedtools intersect -a {params.scratch}/$name.A.bed -b {input.A} > {params.scratch}/$name.AA.bed
		bedtools intersect -a {params.scratch}/$name.B.bed -b {input.B} > {params.scratch}/$name.BB.bed
		bedtools intersect -a {params.scratch}/$name.A.bed -b {input.B} > {params.scratch}/$name.BA.bed
		bedtools intersect -a {params.scratch}/$name.B.bed -b {input.A} > {params.scratch}/$name.AB.bed
		
		AA=$(cat {params.scratch}/$name.AA.bed | awk '{{ {{ s+= int(($3 - ($2) )/50000 + 0.5) }} }}END{{print s}}' 2>&1)
		BB=$(cat {params.scratch}/$name.BB.bed | awk '{{ {{ s+= int(($3 - ($2) )/50000 + 0.5) }} }}END{{print s}}' 2>&1)
		AB=$(cat {params.scratch}/$name.AB.bed | awk '{{ {{ s+= int(($3 - ($2) )/50000 + 0.5) }} }}END{{print s}}' 2>&1)
		BA=$(cat {params.scratch}/$name.BA.bed | awk '{{ {{ s+= int(($3 - ($2) )/50000 + 0.5) }} }}END{{print s}}' 2>&1)
		c_bins=$(cat {input.A} {input.B} | awk '{{ {{ s+= int(($3 - ($2) )/50000 + 0.5) }} }}END{{print s}}' 2>&1)
		all_bins=$(cat {input.sample}  | awk '{{ {{ s+= int(($3 - ($2) )/50000 + 0.5) }} }}END{{print s}}' 2>&1)	
		
		echo "Number of bins at 50000 bin size: "$all_bins >> {params.dir}/cohort/ABcompartments/$name/$name.log
		echo "Number of bins conserved in control samples: "$c_bins >> {params.dir}/cohort/ABcompartments/$name/$name.log
		echo "Number of reciprocal A bins in control and sample "$name" : "$AA >> {params.dir}/cohort/ABcompartments/$name/$name.log
		echo "Number of reciprocal B bins in control and sample "$name" : "$BB >> {params.dir}/cohort/ABcompartments/$name/$name.log
		echo "Number of bin that switched from A to B in sample "$name" : "$AB >> {params.dir}/cohort/ABcompartments/$name/$name.log
		echo "Number of bin that switched from B to A in sample "$name" : "$BA >> {params.dir}/cohort/ABcompartments/$name/$name.log
		
		bedtools intersect -a {params.scratch}/$name.AA.bed -b {params.genes} -loj | awk '{{print $7 }}' | sort | uniq > {params.dir}/cohort/ABcompartments/$name/AA.genes.txt
		bedtools intersect -a {params.scratch}/$name.BB.bed -b {params.genes} -loj | awk '{{print $7 }}' | sort | uniq > {params.dir}/cohort/ABcompartments/$name/BB.genes.txt
		bedtools intersect -a {params.scratch}/$name.AB.bed -b {params.genes} -loj | awk '{{print $7 }}' | sort | uniq > {params.dir}/cohort/ABcompartments/$name/AB.genes.txt
		bedtools intersect -a {params.scratch}/$name.BA.bed -b {params.genes} -loj | awk '{{print $7 }}' | sort | uniq > {params.dir}/cohort/ABcompartments/$name/BA.genes.txt
	
		"""		

### for each chromosome create an A B compartment overview		
rule fanc_compartments_genome:
	input:
		matrix="%s/{sample}/fanc/compartments/matrix/{sample}.500000.matrix"%result_dir
	output:
		plots=directory("%s/{sample}/fanc/compartments/genome/"%result_dir)
	params:
		ref=config['REFERENCE_FILE'],
		prefix=config['PREFIX']
	conda: "envs/fanc.yml"
	shell:
		"""
		mkdir -p {output.plots}
		for chr in $(seq 1 1 22) X Y ; do 
			fancplot  --pdf-text-as-font  -o {output.plots}/{wildcards.sample}_chr$chr.pdf {params.prefix}$chr --plot square {input.matrix} 
		done	
		"""

### target gene region
### 125 bins both ways from gene start
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


### create figures for target gene regions
rule fanc_compartments_region:
	input:
		matrix="%s/{sample}/fanc/compartments/matrix/{sample}.100000.matrix"%result_dir,
		region="%s/cohort/target_regions_r%s.txt"%(config['SCRATCH'],config['RANGE'])
	output:
		directory("%s/{sample}/fanc/compartments/targetRegions/r%s/"%(result_dir, config['RANGE']))
	params:
		genes=config['TARGET_GENES']
	conda: "envs/fanc.yml"
	shell:
		"""
		mkdir -p {output}
		i=0
		IFS=' ' read -r -a genes <<< "{params.genes}"
		while read region ; do
			fancplot  --pdf-text-as-font  -o {output}/{wildcards.sample}_"${{genes[$i]}}".pdf "$region" --plot square {input.matrix}
			i=$(($i+1))
		done < {input.region}
		
		"""
