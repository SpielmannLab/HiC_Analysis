import os

configfile: "Juicer.yml"


### creating juicer top dir
### Juicer top dir should contain certain folders, among them a scripts folder
### the scripts folder depends on the cluster system
rule set_up:
	output: 
		directory("%s/Juicer"%config['SCRATCH'])
	params:
		dir=config['JUICER_INSTALL_DIR'],
		scripts_dir=directory("%s/%s"%(config['JUICER_INSTALL_DIR'],config['CLUSTER']))
	shell:
		"""
		mkdir -p {output}
		ln -s {params.dir}/references {output}/
		ln -s {params.dir}/restriction_sites {output}/
		ln -s {params.dir}/misc {output}/
		ln -s {params.scripts_dir}/scripts {output}/
		"""


############ ALIGNMENT WITHOUT SCAFFOLDING
### this is the "normal" HiC aligning to a reference genome

### get chromsome sizes from reference genome
rule ref_sizes:
	output:
		"%s/ref.sizes"%config['SCRATCH']
	params:
		ref=config['REF_PATH']
	shell:
		"""
		samtools faidx {params.ref} 
		cut -f1-2 {params.ref}.fai > {output}
		"""

### depending on the restriction enzyme determine cut positions in the reference
rule generate_site_positions:
	output:
		"data/juicer/%s_%s.txt"%(config['REF_NAME'],config['ENZYME'] )
	params:
		ref=config['REF_NAME'],
		ref_fasta=config['REF_PATH'],
		enzyme=config['ENZYME'],
		dir=config['JUICER_INSTALL_DIR']
	shell:
		"""
		python {params.dir}/misc/generate_site_positions.py {params.enzyme} {params.ref} {params.ref_fasta}
		mv {params.ref}_{params.enzyme}.txt {output}
		"""

### the alignment and .hic generation
### juicer.sh will start many dispatched jobs, both "longterm"  and "shortterm"
### since longterm (flag -l) is < 3 days we set it to omics luebeck partition shorterm
### since shortterm (flag -q) is < 1 day we set it to omics luebeck partition debug
### as soon as the scripts have finished the files will be in the output folders
rule juicerPipeline:
	input:
		genome=config['REF_PATH'],
		sizes="%s/ref.sizes"%config['SCRATCH'],
		fraq="data/juicer/%s_%s.txt"%(config['REF_NAME'],config['ENZYME'] ),
		R1=config['R1'],
		R2=config['R2'],
		dir=directory("%s/Juicer"%config['SCRATCH'])
	output:
		directory("%s/splits"%config['OUTDIR']),
		directory("%s/aligned"%config['OUTDIR'])
	params:
		enzyme=config['ENZYME'],
		sample=config['NAME'],
		out=config['OUTDIR'],
		scratch=config['SCRATCH']
	shell:
		"""
		mkdir -p {params.scratch}/fastq/
		rm -r {params.scratch}/fastq/*
		ln -s {input.R1} {params.scratch}/fastq/
		ln -s {input.R2} {params.scratch}/fastq/
		
		bash {input.dir}/scripts/juicer.sh -D {input.dir} -d {params.out} -y {input.fraq} -a {params.sample} -z {input.genome} -f -q debug -l shortterm -p {input.sizes}
		"""

############ END OF ALIGNMENT WITHOUT SCAFFOLDING


############ ALIGNMENT WITH SCAFFOLDING
		 
### create a draft fasta file 
### and index with bwa
rule w2rapContigger:
	input:
		R1=config['R1'],
		R2=config['R2']
	output:
		"%s/draft/draft.fa"%config['OUTDIR']
	params:
		img=config['w2rap'],
		out=config['OUTDIR']
	threads: 16
	shell:
		"""	
 		singularity exec {params.img} w2rap-contigger -m 80 -t {threads} -r {input.R1},{input.R2} --dump_all 1 -o {params.out} -p a
 		mv {params.out}/a.lines.fasta {output}
 		bwa index {output}
 		"""

### create enzyme site positions for draft reference 
rule generate_site_positions_draft:
	input:
		draft="%s/draft/draft.fa"%config['OUTDIR']
	output:
		"%s/draft_%s.txt"%(config['OUTDIR'], config['ENZYME'])
	params:
		enzyme=config['ENZYME'],
		dir=config['JUICER_INSTALL_DIR']
	shell:
		"""
		job=$RANDOM
		python {params.dir}/misc/generate_site_positions.py {params.enzyme} $job {input.draft}
		mv ${{job}}_{params.enzyme}.txt {output}
		"""

### get "chromosome chucks" size from draft genome
rule draft_sizes:
	input:
		"%s/draft/draft.fa"%config['OUTDIR']
	output:
		"%s/draft.sizes"%config['OUTDIR']
	shell:
		"""
		samtools faidx {input} 
		cut -f1-2 {input}.fai > {output}
		"""

### align to the draft genome
### -D Juicer scripts parent dir, which should have scripts/ references/ and restriction_sites/ underneath it
### -d topdir topDir/fastq must contain the fastq file
rule juicerAssembly:
	input:
		R1=config['R1'],
		R2=config['R2'],
		draft="%s/draft/draft.fa"%config['OUTDIR'],
		sizes="%s/draft.sizes"%config['OUTDIR'],
		positions="%s/draft_%s.txt"%(config['OUTDIR'], config['ENZYME']),
		dir=directory("%s/Juicer"%config['SCRATCH'])
	output:
		directory("%s/juicerAssembly/fastq"%config['OUTDIR'])
	params:
		enzyme=config['ENZYME'],
		out=config['OUTDIR'],
		scratch=config['SCRATCH']
	shell:
		"""
		mkdir -p {params.out}/juicerAssembly/fastq
		ln -s {input.R1} {params.out}/juicerAssembly/fastq/
		ln -s {input.R2} {params.out}/juicerAssembly/fastq/
		bash {input.dir}/scripts/juicer.sh -s {params.enzyme} -z {input.draft} -y {input.positions} -D {input.dir} -d {params.out}/juicerAssembly -p {input.sizes} -e -j --assembly -q shortterm -l longterm
		"""



### create 3D de novo .assembly and .hic file with 3d-dna
### both .hic and .assembly are required for interactive session in JuiceBox
rule asmPipeline:
	input:
		reads="%s/juicerAssembly/aligned/merged_nodups.txt"%config['OUTDIR'],
		draft="%s/draft/draft.fa"%config['OUTDIR']
	output:
		assembly="%s/3DNA_Assembly/draft.assembly"%config['OUTDIR'],
		final="%s/hic_format/draft.final.hic"%config['OUTDIR'],
		resolved="%s/hic_format/draft.rawchrom.hic"%config['OUTDIR']
	params:
		scratch=config['SCRATCH'],
		dir=config['3D-DNA_INSTALL_DIR'],
		out=config['OUTDIR']
	shell:
		"""
		mkdir -p {params.scratch}/asm
		origin=$PWD
		cd {params.scratch}/asm
		bash {params.dir}/run-asm-pipeline.sh -q 0 -s stage $origin/{input.draft} $origin/{input.reads}  
		awk -f {params.dir}/utils/generate-assembly-file-from-fasta.awk $origin/{input.draft} > $origin/{output.assembly}
		bash {params.dir}/visualize/run-assembly-visualizer.sh -q 0 -m temp.draft.asm_mnd.txt  $origin/{output.assembly} draft.mnd.txt
		mv {draft_0,draft.polished,draft.split}.hic {params.out}/hic_format/
		mv draft.rawchrom.hic {output.final}
		mv draft.1.hic {output.resolved}
		"""








