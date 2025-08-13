# Scripts to compare two HiC matrices using fanc compare:

1. Calculate difference matrix. Test output stored at: `/data/humangen_external/test_area/varun/2025.07_fanc_compare/`

```bash
srun -p shortterm -c 1 --mem 100GB --pty bash
cd $SCRATCH && \
cp /data/humangen_external/HiC/steinhaeuser_hic/data_n_results/DNMT3A_clones/2023.10_fastq2HiC/juicer_1.19.02/DNMT3Amut55.allValidPairs.hic ./ && \
cp /data/humangen_external/HiC/steinhaeuser_hic/data_n_results/DNMT3A_clones/2023.10_fastq2HiC/juicer_1.19.02/DNMT3Awtc13.allValidPairs.hic ./

conda activate fanc
fanc compare -c diff -u  DNMT3Amut55.allValidPairs.hic@100000 DNMT3Awtc13.allValidPairs.hic@100000 /data/humangen_external/test_area/varun/2025.07_fanc_compare/mut55_minus_ct13_100k.out
```

2. Plot some regions, that were identified to be different based on insulation scores.

```bash
# chr1_part_small
fancplot 1:500kb-10mb -o /data/humangen_external/test_area/varun/2025.07_fanc_compare/mut55_minus_ct13_500k_chr1_part_small.png -p triangular -vmin -5000 -vmax 5000 -c RdBu /data/humangen_external/test_area/varun/2025.07_fanc_compare/mut55_minus_ct13_500k.out 
fancplot 1:500kb-10mb -o /data/humangen_external/test_area/varun/2025.07_fanc_compare/mut55_chr1_part_small_vc.png -p triangular -l -vmin 1 -vmax 10000 DNMT3Amut55.allValidPairs.hic@500000@VC
fancplot 1:500kb-10mb -o /data/humangen_external/test_area/varun/2025.07_fanc_compare/wtc13_chr1_part_small_vc.png -p triangular -l -vmin 1 -vmax 10000 DNMT3Awtc13.allValidPairs.hic@500000@VC



# chr10_part_small
fancplot 10:500kb-100mb -o /data/humangen_external/test_area/varun/2025.07_fanc_compare/mut55_minus_ct13_500k_chr10_part_small.png -p triangular -vmin -5000 -vmax 5000 -c RdBu /data/humangen_external/test_area/varun/2025.07_fanc_compare/mut55_minus_ct13_500k.out 
fancplot 10:500kb-100mb -o /data/humangen_external/test_area/varun/2025.07_fanc_compare/mut55_chr10_part_small.png -p triangular -vmin 0 -vmax 10000 DNMT3Amut55.allValidPairs.hic@500000@VC
fancplot 10:500kb-100mb -o /data/humangen_external/test_area/varun/2025.07_fanc_compare/wtc13_chr10_part_small.png -p triangular -vmin 0 -vmax 10000 DNMT3Awtc13.allValidPairs.hic@500000@VC
```

3. Calculate the AB compartments and get the eigenvector values. Test output stored at `/data/humangen_external/test_area/varun/2025.07_fanc_ab_compartments`

```bash
srun -p shortterm -c 1 --mem 100GB --pty bash
cd $SCRATCH && \
cp /data/humangen_external/HiC/steinhaeuser_hic/data_n_results/DNMT3A_clones/2023.10_fastq2HiC/juicer_1.19.02/DNMT3Amut55.allValidPairs.hic ./

conda activate fanc
fanc compartments -v DNMT3Amut55_25kb.ev.txt DNMT3Amut55.allValidPairs.hic@25kb DNMT3Amut55_25kb.ab
cp DNMT3Amut55_* /data/humangen_external/test_area/varun/2025.07_fanc_ab_compartments/

# Now plot the AB compartments
fancplot 1:500kb-100mb -o /data/humangen_external/test_area/varun/2025.07_fanc_ab_compartments/DNMT3Amut55_1mb_AB_chr1_part_small.png -p square -c RdBu_r DNMT3Amut55_1mb.ab -p line DNMT3Amut55_1mb.ev.txt
fancplot 10:500kb-100mb -o /data/humangen_external/test_area/varun/2025.07_fanc_ab_compartments/DNMT3Amut55_1mb_AB_chr10_part_small.png -p square -c RdBu_r DNMT3Amut55_1mb.ab -p line DNMT3Amut55_1mb.ev.txt
fancplot 2:500kb-100mb -o /data/humangen_external/test_area/varun/2025.07_fanc_ab_compartments/DNMT3Amut55_1mb_AB_chr2_part_small.png -p square -c RdBu_r DNMT3Amut55_1mb.ab -p line DNMT3Amut55_1mb.ev.txt
fancplot 2:1b-100mb -o /data/humangen_external/test_area/varun/2025.07_fanc_ab_compartments/DNMT3Amut55_25kb_AB_chr2_part_small.png -p square -vmin -0.1 -vmax 0.1 -c RdBu_r DNMT3Amut55_25kb.ab -p line DNMT3Amut55_25kb.ev.txt

```
