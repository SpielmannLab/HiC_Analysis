################################################
#
# This workflow uses fithichip to generate chip peaks and tracks from HiCHiP data.
# Date: September, 2021
#
################################################


import os
    
configfile: "config.yml"

def modeInt2name(argument):
    switcher = {
        1: "Peak2Peak",
        2: "Peak2NonPeak",
        3: "Peak2ALL",
        4: "ALL2ALL",
        5: "['Peak2Peak', 'Peak2NonPeak', 'Peak2ALL', 'ALL2ALL']",
    }
    return switcher.get(argument, "Invalid Interaction type - 1: peak to peak 2: peak to non peak 3: peak to all (default) 4: all to all 5: everything from 1 to 4.")
    
modeName=modeInt2name(config['MODE'])

#"%s/pooled/allbedgraph.bdg"%config['OUTDIR']
# TODO change path /P2PBckgr_0/Coverage_Bias/ to /Coverage_Bias/ only when mode 4 is used ...

rule all:
	input:
		expand("%s/{sample}/fithichip/FitHiChIP_%s_b%s_L%s_U%s/Coverage_Bias/FitHiC_BiasCorr/FitHiChIP.interactions_FitHiC_Q%s.bed"%(config['OUTDIR'], modeName, config['BINSIZE'], config['LOW_TH'], config['UP_TH'], config['QVALUE']), sample=config['SAMPLES']),
		"%s/pooled/fithichip/FitHiChIP_%s_b%s_L%s_U%s/Coverage_Bias/FitHiC_BiasCorr/FitHiChIP.interactions_FitHiC_Q%s.bed"%(config['OUTDIR'], modeName, config['BINSIZE'], config['LOW_TH'], config['UP_TH'], config['QVALUE'])
		


rule get_tracks:
	output:
		"%s/{sample}/hichip-peaks/{sample}bedgraph.bdg"%config['OUTDIR']
	params:
		resfraq=config['RESFRAQ'],
		sizes=config['REF_SIZES'],
		outdir="%s/{sample}/hichip-peaks/"%config['OUTDIR'],
		indir="%s/{sample}/hic_results/data/"%config['HICPRO_INDIR']
	conda:
		"envs/hichip-peaks.yml"
	shell:
		"""
		peak_call -i {params.indir} -o {params.outdir} -r {params.resfraq} -a {params.sizes} -w $(nproc) -p {wildcards.sample} -d -x -k
		bgzip < {output} > {output}.gz
		tabix -s 1 -b 2 -e 3 {output}.gz
		"""

rule merge_tracks:
	input:
		expand("%s/{sample}/hichip-peaks/{sample}bedgraph.bdg"%config['OUTDIR'], sample=config['SAMPLES'])
	output:
		"%s/pooled/allbedgraph.bdg"%config['OUTDIR']
	shell:
		"""
		paste {input} | awk -v OFS='\\t'  '{{ print $1, $2, $3, $4 + $8 + $12 + $16 + $20; }}' >  {output}	
		bgzip < {output} > {output}.gz
		tabix -s 1 -b 2 -e 3 {output}.gz
		"""

rule call_peaks:
	output:
		"%s/{sample}/peak_calls/MACS2_input.bed"%config['OUTDIR'],
		"%s/{sample}/peak_calls/MACS2_ExtSize/out_macs2_peaks.narrowPeak"%config['OUTDIR']
	params:
		outdir="%s/{sample}/peak_calls/"%config['OUTDIR'],
		indir="%s/{sample}/hic_results/"%config['HICPRO_INDIR'],
		macs2=config['MACS2_CONFIG'],
		RL=config['READ_LENGTH'],
		image=config['IMAGE']
	shell:
		"""
		cp config.yml {params.outdir}
		singularity exec {params.image} cp /FitHiChIP/Imp_Scripts/PeakInferHiChIP.sh {params.outdir}
		sed -i -r "s|MACS2ParamStr='.*|MACS2ParamStr='{params.macs2}'|g" {params.outdir}PeakInferHiChIP.sh
		sed -i -r "s|ReadLength=75|ReadLength={params.RL}|g" {params.outdir}PeakInferHiChIP.sh
		singularity exec {params.image} bash {params.outdir}PeakInferHiChIP.sh -H {params.indir} -D {params.outdir}

		"""

rule ms_call_peaks:
	input:
		expand("%s/{sample}/peak_calls/MACS2_ExtSize/out_macs2_peaks.narrowPeak"%config['OUTDIR'], sample=config['SAMPLES'])
	output:
		"%s/pooled/mspc/ConsensusPeaks_mspc_peaks.txt"%config['OUTDIR'],
		"%s/pooled/mspc/ConsensusPeaks.bed"%config['OUTDIR']
	params:
		outdir="%s/"%config['OUTDIR']
	shell:
		"""
		./../utils/mspc/mspc -i $(ls {params.outdir}/*/peak_calls/MACS2_ExtSize/out_macs2_peaks.narrowPeak) -r bio -w 1e-4 -s 1e-8 -o {params.outdir}/pooled/mspc
		"""

rule fithichip:
	input:
		"%s/{sample}/peak_calls/MACS2_ExtSize/out_macs2_peaks.narrowPeak"%config['OUTDIR']
	output:
		"%s/{sample}/fithichip/FitHiChIP_%s_b%s_L%s_U%s/Coverage_Bias/FitHiC_BiasCorr/FitHiChIP.interactions_FitHiC_Q%s.bed"%(config['OUTDIR'], modeName, config['BINSIZE'], config['LOW_TH'], config['UP_TH'], config['QVALUE']),
	params:
		sizes=config['REF_SIZES'],
		outdir="%s/{sample}/fithichip/"%config['OUTDIR'],
		indir=config['HICPRO_INDIR'],
		mode=config['MODE'],
		bs=config['BINSIZE'],
		low=config['LOW_TH'],
		up=config['UP_TH'],
		qval=config['QVALUE'],
		samples=config['SAMPLES'],
		image=config['IMAGE']
	shell:
		"""
		singularity exec {params.image} cp /FitHiChIP/FitHiChIP_HiCPro.sh {params.outdir}

		sed -i -r "s|./Analysis/|/FitHiChIP/Analysis/|g" {params.outdir}FitHiChIP_HiCPro.sh
		sed -i -r "s|./src/|/FitHiChIP/src/|g" {params.outdir}FitHiChIP_HiCPro.sh
		sed -i -r "s|/FitHiChIP/src/FitHiC_SigInt.r|$PWD/../scripts/FitHiC_SigInt.r|g" {params.outdir}FitHiChIP_HiCPro.sh
		
		singularity exec {params.image} cp /FitHiChIP/configfile_BiasCorrection_CoverageBias {params.outdir}
		
		sed -i -r "s|ValidPairs=.*|ValidPairs=$PWD/{params.indir}/{wildcards.sample}/hic_results/data/{wildcards.sample}.allValidPairs|g" {params.outdir}configfile_BiasCorrection_CoverageBias
		sed -i -r "s|PeakFile=.*|PeakFile=$PWD/{input}|g" {params.outdir}configfile_BiasCorrection_CoverageBias
		sed -i -r "s|OutDir=.*|OutDir=$PWD/{params.outdir}|g" {params.outdir}configfile_BiasCorrection_CoverageBias
		sed -i -r "s|ChrSizeFile=.*|ChrSizeFile=$PWD/{params.sizes}|g" {params.outdir}configfile_BiasCorrection_CoverageBias
		sed -i -r "s|IntType=.*|IntType={params.mode}|g" {params.outdir}configfile_BiasCorrection_CoverageBias		
		sed -i -r "s|QVALUE=.*|QVALUE={params.qval}|g" {params.outdir}configfile_BiasCorrection_CoverageBias		
		sed -i -r "s|BINSIZEr=.*|BINSIZE={params.bs}|g" {params.outdir}configfile_BiasCorrection_CoverageBias
		sed -i -r "s|LowDistThr=.*|LowDistThr={params.low}|g" {params.outdir}configfile_BiasCorrection_CoverageBias
		sed -i -r "s|UppDistThr=.*|UppDistThr={params.up}|g" {params.outdir}configfile_BiasCorrection_CoverageBias
		
		singularity exec {params.image} bash {params.outdir}FitHiChIP_HiCPro.sh  -C {params.outdir}configfile_BiasCorrection_CoverageBias
		
		"""

# rule merge_validPairs:
# 	output:
# 		"%s/pooled/fithichip/pooled.allValidPairs"%config['OUTDIR']
# 	params:
# 		outdir="%s/pooled/fithichip/"%config['OUTDIR'],
# 		samples=config['SAMPLES'],
# 		indir=config['HICPRO_INDIR']
# 	shell:
# 		"""
# 		arr=($"{params.samples}")
# 		for s in $arr ; do 
# 			for f in {params.indir}/$s/hic_results/data/*.validPairs ; do
# 				cat $f >> {params.outdir}/pooled.allValidPairs
# 			done
# 		done
# 		"""

rule ms_fithichip:
	input:
		peaks="%s/pooled/mspc/ConsensusPeaks_mspc_peaks.txt"%config['OUTDIR'],
		validPairs="%s/{sample}/hic_results/data/{sample}.allValidPairs"%config['OUTDIR']
	output:
		"%s/pooled/fithichip/FitHiChIP_%s_b%s_L%s_U%s/Coverage_Bias/FitHiC_BiasCorr/FitHiChIP.interactions_FitHiC_Q%s.bed"%(config['OUTDIR'], modeName, config['BINSIZE'], config['LOW_TH'], config['UP_TH'], config['QVALUE']),
	params:
		sizes=config['REF_SIZES'],
		outdir="%s/pooled/fithichip/"%config['OUTDIR'],
		indir=config['HICPRO_INDIR'],
		mode=config['MODE'],
		bs=config['BINSIZE'],
		low=config['LOW_TH'],
		up=config['UP_TH'],
		qval=config['QVALUE'],
		samples=config['SAMPLES'],
		image=config['IMAGE']
	shell:
		"""
		singularity exec {params.image} cp /FitHiChIP/FitHiChIP_HiCPro.sh {params.outdir}

		sed -i -r "s|./Analysis/|/FitHiChIP/Analysis/|g" {params.outdir}FitHiChIP_HiCPro.sh
		sed -i -r "s|./src/|/FitHiChIP/src/|g" {params.outdir}FitHiChIP_HiCPro.sh
		sed -i -r "s|/FitHiChIP/src/FitHiC_SigInt.r|$PWD/../scripts/FitHiC_SigInt.r|g" {params.outdir}FitHiChIP_HiCPro.sh
		
		singularity exec {params.image} cp /FitHiChIP/configfile_BiasCorrection_CoverageBias {params.outdir}
		
		sed -i -r "s|ValidPairs=.*|ValidPairs=$PWD/{input.validPairs}|g" {params.outdir}configfile_BiasCorrection_CoverageBias
		sed -i -r "s|PeakFile=.*|PeakFile=$PWD/{input.peaks}|g" {params.outdir}configfile_BiasCorrection_CoverageBias
		sed -i -r "s|OutDir=.*|OutDir=$PWD/{params.outdir}|g" {params.outdir}configfile_BiasCorrection_CoverageBias
		sed -i -r "s|ChrSizeFile=.*|ChrSizeFile=$PWD/{params.sizes}|g" {params.outdir}configfile_BiasCorrection_CoverageBias
		sed -i -r "s|IntType=.*|IntType={params.mode}|g" {params.outdir}configfile_BiasCorrection_CoverageBias		
		sed -i -r "s|QVALUE=.*|QVALUE={params.qval}|g" {params.outdir}configfile_BiasCorrection_CoverageBias		
		sed -i -r "s|BINSIZEr=.*|BINSIZE={params.bs}|g" {params.outdir}configfile_BiasCorrection_CoverageBias
		sed -i -r "s|LowDistThr=.*|LowDistThr={params.low}|g" {params.outdir}configfile_BiasCorrection_CoverageBias
		sed -i -r "s|UppDistThr=.*|UppDistThr={params.up}|g" {params.outdir}configfile_BiasCorrection_CoverageBias
		
		singularity exec {params.image} bash {params.outdir}FitHiChIP_HiCPro.sh  -C {params.outdir}configfile_BiasCorrection_CoverageBias
		
		"""
			
rule bed2ucsc:
	input:
		"%s/pooled/fithichip/FitHiChIP_%s_b%s_L%s_U%s/Coverage_Bias/FitHiC_BiasCorr/FitHiChIP.interactions_FitHiC_Q%s.bed"%(config['OUTDIR'], modeName, config['BINSIZE'], config['LOW_TH'], config['UP_TH'], config['QVALUE'])
	output:
		"%s/ucsc/loop.bb"%config['OUTDIR']
	params:
		outdir="%s/ucsc/"%config['OUTDIR'],
		sizes=config['REF_SIZES']
	shell:
		"""
		awk -F'[\\t]' '{{if (NR>1) {{if ($NF>0) {{a=int(-log($NF)/log(10)) + 1}} else {{a=500}}; b=100; print $1"\\t"(($2+$3)/2-1)"\\t"(($5+$6)/2+1)"\\t.\\t"a"\\t"b"\\t.\\t#0000FF\\t"$1"\\t"(($2+$3)/2-1)"\\t"(($2+$3)/2+1)"\\t.\\t.\\t"$4"\\t"(($5+$6)/2-1)"\\t"(($5+$6)/2+1)"\\t.\\t."}}}}' {input} | sort -k1,1 -k2,2n -k3,3n > {params.outdir}/out_FitHiChIP_converted.bed
		../utils/bedToBigBed -as=../utils/interact.as -type=bed5+13 {params.outdir}/out_FitHiChIP_converted.bed {params.sizes} {output}
		"""