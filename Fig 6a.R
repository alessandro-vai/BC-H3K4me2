library(TCGAbiolinks)

query.hg38 <- GDCquery(
    project= "TCGA-BRCA",
    data.category = "Transcriptome Profiling", 
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts", 
    sample.type = "Primary Tumor"
)
GDCdownload(query.hg38)
data.hg38 <- GDCprepare(query.hg38, save = TRUE,
                        save.filename = "TCGA-BRCA-STAR.RData",
                        remove.files.prepared = TRUE)

load("TCGA-BRCA-STAR.RData")
data <- data[,! is.na(data$paper_BRCA_Subtype_PAM50)]
data <- data[data@rowRanges@seqnames@values != "chrM",]
data <- data[!rowData(data)$gene_type %in% c("rRNA", "rRNA_pseudogene",
                                             "scaRNA", "scRNA",
                                             "snoRNA", "snRNA",
                                             "ribozyme", "miRNA",
                                             "Mt_rRNA", "Mt_tRNA"),]
data <- data[, data$paper_BRCA_Subtype_PAM50 != "Normal"]

dds <- DESeqDataSetFromMatrix(assay(data), 
                              colData=colData(data), 
                              ~paper_BRCA_Subtype_PAM50)
smallestGroupSize <- data$paper_BRCA_Subtype_PAM50 |> table() |> min()
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
rm(keep, smallestGroupSize)
dds <- DESeq(dds, parallel = T,
             BPPARAM = BiocParallel::MulticoreParam(8))

rowdata <- rowData(data)[, c("gene_id", "gene_name")] |> as_tibble()
# This table was used as input for GraphPad to generate Fig 6a
log2(counts(dds, norm=T) + 1) |> as_tibble(rownames="gene_id") |>
    inner_join(rowdata) |> filter(gene_name == "KMT2B") |> 
    select(-gene_id,-gene_name) |> t() |> 
    merge(colData(data)[, c("barcode", "paper_BRCA_Subtype_PAM50")], by=0) |> 
    select(-Row.names)

# Stats
getKMT2BTable <- function(subtype){
    res <- results(dds, c("paper_BRCA_Subtype_PAM50", "Basal", subtype),
                   parallel = T, BPPARAM = BiocParallel::MulticoreParam(8))
    res |> as.data.frame() |> as_tibble(rownames="gene_id") |>
        inner_join(rowdata) |> filter(gene_name == "KMT2B") |> 
        mutate(Comparison = paste0("Basalvs", subtype))
}

lapply(c("LumA", "LumB", "Her2"), getKMT2BTable) |> 
    bind_rows()




