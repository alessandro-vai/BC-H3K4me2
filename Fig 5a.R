library(ComplexHeatmap) |> suppressPackageStartupMessages()


Zrna <- readRDS("Zscore_combat.rds") # From 'TCGA Analysis.Rmd'
candidate_genes <- c(#"SLAMF7",
    "AIM2",
    "IGF2BP3", "PLA2G4A",
    "GLS", "PLCG2", "GLIPR2")
Zrna <- Zrna[rownames(Zrna) %in% candidate_genes,]
annotRNA <- readRDS("BRCA_SE_filtered.rds")$Subtype # From 'TCGA Analysis.Rmd'
Zprot <- readRDS("ZscoreCPTAC.rds") # From 'CPTAC Analysis.R'
Zprot <- Zprot[rownames(Zprot) %in% candidate_genes,]
Zprot <- rbind(Zprot, AIM2 = rep(NA, ncol(Zprot)))
Zprot <- Zprot[rownames(Zrna),]
annotPROT <- ifelse(grepl("LuA", colnames(Zprot)), "LumA", "TN")

makeHA <- function(annot_vec){
    ha <-  HeatmapAnnotation(subtype = annot_vec,
                             col = list(subtype = c("LumA" = "#80ffff",
                                                    "TN" = "black")))
    return(ha)
}

makeHeatmap <- function(Z, ha, width_num, title) {
    h <- Heatmap(Z,
                 name = "z-score",
                 cluster_rows = F,
                 cluster_columns = T,
                 clustering_method_columns = "complete",
                 show_column_names = F,
                 row_names_side = "left",
                 top_annotation = ha,
                 na_col = "gray90",
                 col = circlize::colorRamp2(c(-4, 0, 4), c("blue", "white", "red"), space = "RGB"),
                 width = unit(width_num, "cm"),
                 column_title = title)
    return(h)
    
}

hrna <- makeHeatmap(Zrna, makeHA(annotRNA), 12, "TCGA")
hprot <- makeHeatmap(Zprot, makeHA(annotPROT), 6, "CPTAC")
hrna+hprot





