# conda activate MegaBundle
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))

args = commandArgs(trailingOnly = TRUE)

deseq_result <- read.table("/Users/sreenivasan/Documents/Works/steinhaeuser_hic/data_n_results/2023.09 DESeq data from DavidBaden/DNMT3A_wo_outlier/DESeq2_DNMT3A_wo_C14.tsv", sep="\t", header = TRUE)
group1_name = "DNMT3A"
group1_genes <- read.table("/Users/sreenivasan/Documents/Works/steinhaeuser_hic/data_n_results/2023.09_diffPeakachu/diffpeakachu/DNMT3A.unique.genes.txt") %>%
    .[,1]
group2_name = "WT"
group2_genes <- read.table("/Users/sreenivasan/Documents/Works/steinhaeuser_hic/data_n_results/2023.09_diffPeakachu/diffpeakachu/WT.unique.genes.txt") %>%
    .[,1]

deseq_hic_comparison <- deseq_result %>%
    filter(genes %in% c(group1_genes, group2_genes)) %>%
    # filter(!(genes %in% group1_genes & genes %in% group2_genes)) %>% # Remove genes that are common in both
    mutate(unique_in = case_when( ((genes %in% group1_genes) & (genes %in% group2_genes)) ~ paste0(group1_name, "&", group2_name),
				 genes %in% group1_genes ~ group1_name,
				 genes %in% group2_genes ~ group2_name)) %>%
    pivot_longer(cols = starts_with("meanCounts"),
		 names_to = "Group",
		 values_to = "Expression")

write.table(file = "Mean_expression_of_uniq_genes.tsv", arrange(deseq_hic_comparison, padj), sep="\t")

# custom facet labels
#### !!!! Warning - this has not been updated for plotting when group1&group2 is included
unique_in.labs <- paste0("Genes in ", c(group1_name, group2_name), " specfic loop")
names(unique_in.labs) <- c(group1_name, group2_name)

ggplot(deseq_hic_comparison, aes(x = Group, y = Expression, color = Group)) + 
    geom_boxplot() +
    geom_jitter(alpha = 0.5) +
    facet_wrap(.~unique_in,
        labeller = labeller(unique_in = unique_in.labs)) +
    scale_y_log10() +
    scale_x_discrete(breaks = c("meanCounts_mut", "meanCounts_wt"), labels = c("DNMT3A", "WT")) +
    theme(legend.position = "none")
ggsave(filename = "Mean_expression_of_unique_genes.png", width = 10, height = 7)

ggplot(deseq_hic_comparison, aes(y = log2FoldChange_.normal., x = unique_in)) + 
    geom_boxplot() +
    geom_jitter(alpha = 0.5) +
ggsave(filename = "test.pdf", width =7 , height = 7)
