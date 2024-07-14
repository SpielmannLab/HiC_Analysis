# TAD calling and TAD boundary detection

## Literature and conclusions

Here are some important resources and papers that deal with TAD and boundary detection.

- An up-to-date list of resources [Github page](https://github.com/mdozmorov/HiC_tools?tab=readme-ov-file#tad-detection-benchmarking)
- A comparison of 22 tools [Zufferey et al, Genome Biology, 2018](https://doi.org/10.1186/s13059-018-1596-9)
- An online TAD caller with multiple algorithms [Higgins et al, BMC Bioinformatics, 2022](https://doi.org/10.1186/s12859-022-05020-2)

Based on the comparison bapser by Zufferey et al, we should implement the following tools:

PCA grouping in Figure 4f of Zufferey et al. Ordered based on further assessments is below. The goal is to implement all of these. [x] mark shows the implemented ones.

Group 1:

1. [Arrowhead](https://github.com/aidenlab/juicer/wiki/Arrowhead#defaults) (need high resolution, best at 10kb resolution) [x]
2. Insulation Score (implementation using [FANC](https://fan-c.readthedocs.io/en/latest/fanc-executable/fanc-analyse-hic/domains.html))? [ ]

Group 2:

1. TopDom [ ]
2. CHDF (need high resolution) [ ]
3. CaTCH [ ]
4. HiCseg [ ]

Not assessed by the authors?:

- HiTAD [ ]

### TODO

The online TAD caller by Higgins et al requires data in a particular format. And is quite restrictive when uploading HiC files (500MB limit). Need to generate a pipeline to export hic files to mcool format etc. Eg using HiCLift

## Usage

You can also use juicer_tools validate command to check which normalizations are present
Fill out the information in the [parameters file](./call_tad_params.yaml) and submit the following script from OMICS headnode:

    sbatch call_tad_sbatch.sh

or, if you want to resume a previous job, then:

    sbatch call_tad_sbatch.sh -resume

The output for each of the implemented algorithms will be found at the `outdir/<algorithm>/`. Please also check the log file.

## Notes
