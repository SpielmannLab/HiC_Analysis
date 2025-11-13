#installation
mamba create --name hicexplorer hicexplorer=3.7.2 -c bioconda -c conda-forge

hicConvertFormat -m /data/humangen_mouse/STIGMA2/hiC/GSM3262966_D80_HiC_Rep1.hic --inputFormat hic --outputFormat cool -o GSM3262966_50k.cool --resolution 50000 --load_raw_values

hicNormalize --matrices GSM3262966_50k.cool GSM3262967_50k.cool --normalize smallest -o GSM3262966_50k_norm.cool GSM3262967_50k_norm.cool

hicCorrectMatrix diagnostic_plot --matrix GSM3262967_50k_norm.cool -o GSM3262967_50k_norm.png

hicCorrectMatrix diagnostic_plot --matrix GSM3262966_50k_norm.cool -o GSM3262966_50k_norm.png


hicCorrectMatrix correct --matrix GSM3262966_50k_norm.cool --correctionMethod KR --outFileName GSM3262966_50k_norm_corrected.cool --filterThreshold -3 3

hicCorrectMatrix correct --matrix GSM3262967_50k_norm.cool --correctionMethod KR --outFileName GSM3262967_50k_norm_corrected.cool --filterThreshold -3 3

hicPlotDistVsCounts --matrices GSM3262966_50k_norm_corrected.cool GSM3262967_50k_norm_corrected.cool -o plot_vs_counts.png

hicCorrelate --log1p --matrices GSM3262966_50k_norm_corrected.cool GSM3262967_50k_norm_corrected.cool --range 20000:500000 -oh between_matrix_cor_h.png -os between_matrix_cor_s.png

hicSumMatrices -m GSM3262966_50k_norm_corrected.cool GSM3262967_50k_norm_corrected.cool -o GSE116862_50k_sum.cool

hicFindTADs -m GSE116862_50k_sum.cool --outPrefix GSE116862_50k_sum --numberOfProcessors 16 --correctForMultipleTesting None
