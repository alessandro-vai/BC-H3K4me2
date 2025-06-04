library(tidyverse) |> suppressPackageStartupMessages()
library(clusterProfiler)

# DEGs from 'TCGA Analysis.Rmd'
TCGA_UP <- read.delim("TCGA_UP.txt",
                         header = F) |> pull(1)
TN_only_prom <- read.delim("TN_only_prom_genes.txt",
                           header = F) |> pull(1)
TN_only_enh <- read.delim("peaks_to_gene/ATAC_Enhancer/TN_only_ATAC_Enhancer.txt",
                           header = F) |> pull(1)
TN_only <- union(TN_only_prom, TN_only_ATAC)
rm(TN_only_prom, TN_only_ATAC)
TN_only_tcga <- TN_only |> intersect(RNA_seq_UP)
# DEPs from 'TCGA Analysis.Rmd'
CPTAC_UP <- read.table("CPTAC_UP.txt")$V1
TN_only_tcga_cptac <- intersect(TN_only_tcga, CPTAC_UP)

#MSigDb
C2_t2g <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2") |>  
    dplyr::select(gs_name, gene_symbol)

ck <- compareCluster(geneCluster = list(TCGA = TN_only_tcga,
                                        `TCGA+CPTAC` = TN_only_tcga_cptac),
                     fun = enricher,
                     TERM2GENE=C2_t2g)

dotplot(ck, showCategory=10) + coord_flip() + scale_x_discrete(limits=rev) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          text = element_text(family = "Arial", size = 14)) +
    xlab("")



