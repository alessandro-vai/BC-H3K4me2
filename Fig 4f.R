library(tidyverse) |> suppressPackageStartupMessages()

TN_only_prom <- read.table("TN_only_prom_genes.txt")$V1
TN_only_enh <- read.table("TN_only_Enhancer.txt")$V1
TN_only <- union(TN_only_prom, TN_only_enh)
LuA_only_prom <- read.table("LuA_only_prom_genes.txt")$V1
LuA_only_enh <- read.table("LuA_only_Enhancer.txt")$V1
LuA_only <- union(LuA_only_prom, LuA_only_enh)
rm(TN_only_prom, TN_only_enh, LuA_only_prom, LuA_only_enh)
tmp <- intersect(TN_only, LuA_only) 
TN_only <- setdiff(TN_only, tmp)
LuA_only <- setdiff(LuA_only, tmp)
rm(tmp)

# From 'TCGA Analysis.Rmd'
batch_dds <- readRDS("batch_dds.rda")
res <- results(batch_dds, alpha=0.05)

p <- res |> as.data.frame() |>
    merge(rowData(brca_se)[, c("gene_id", "gene_name")], by=0) |> 
    as.data.frame() |> 
    select(-Row.names) |> 
    mutate(Legend = case_when(gene_name %in% TN_only ~ "TN only",
                              gene_name %in% LuA_only ~ "LuA only",
                              TRUE ~ "Other"),
           isReg = if_else(abs(log2FoldChange) >= 1 & padj <= 0.05, "Reg", "NR"),
           Direction = case_when(isReg == "Reg" & log2FoldChange < 0 ~ "Down",
                                 isReg == "Reg" & log2FoldChange > 0 ~ "Up",
                                 TRUE ~ "NR"),
           Legend = case_when(Legend == "TN only" & Direction == "Up" ~ "TN only",
                              Legend == "LuA only" & Direction == "Down" ~ "LuA only",
                              isReg == "Reg" ~ "Reg",
                              TRUE ~ "NR")) |> 
    ggplot(aes(x = log2FoldChange, y = -log10(padj), color=Legend, alpha = Legend, size = Legend)) +
    geom_point(show.legend = FALSE) +
    scale_alpha_manual(values = c(1, 0.5, 0.5, 1)) + 
    scale_size_manual(values = c(2.5, 1.5, 1.5, 2.5)) + 
    theme(text = element_text(family = "ArialMT", size = 14)) +
    theme_classic() +
    xlim(-11, 11) +
    scale_color_manual(values = c('#80ffff', "gray", "gray40","black"))
