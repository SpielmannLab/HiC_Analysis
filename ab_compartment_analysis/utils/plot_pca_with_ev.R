#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(yaml)
library(tibble)
library(purrr)

params = read_yaml("input_files_for_pca.yaml")
samples = data.frame()
ev_entry <- list()
for (sample in params[["samples"]]) {
  entry <- data.frame(sample)
  samples <- rbind(samples, entry)
  ev_entry[[entry$name]] <- read.table(entry$ev_file, sep = "\t", header = FALSE) |>
    setNames(c("chr", "start", "end", "domain", entry$name , "x")) |>
    mutate(coord = paste0(chr,":",start,"-",end)) |>
    select("coord", entry$name)
}
samples <- samples |> select(-ev_file)
ev <- reduce(.x =ev_entry, .f = left_join, by = "coord") |>
  select(-coord) |>
  as.matrix() |>
  t()

pca <- prcomp(ev, center = TRUE, scale = FALSE)
dimred <- pca$x |>
  data.frame() |>
  rownames_to_column("name")
variances <- pca$sdev
combined <- left_join(samples, dimred)

plot <- ggplot(combined, aes(x = PC1, y = PC2, color = condition1, shape = condition2)) +
  geom_point(size = 4, stroke = 1) +
  scale_shape_manual(breaks = params$condition2_breaks, values = params$condition2_shapes)
ggsave(plot, filename = "PCA_DNMT3A_IDH1_AB_eigen_vector.png")


