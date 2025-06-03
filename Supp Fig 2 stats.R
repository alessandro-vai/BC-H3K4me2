# se object was generated as first chunk (make se) in Fig 1b.Rmd 

mm0 <- model.matrix(~0+Subtype, data = colData(se))
colnames(mm0) <- str_remove(colnames(mm0), "Subtype")
contrasts_ <- apply(combn(colnames(mm0), 2), 2, paste0, collapse="-")
contrasts.matrix <- makeContrasts(contrasts = contrasts_, levels = colnames(mm0))
limma.res <- lmFit(assay(se), mm0) |>
    contrasts.fit(contrasts.matrix) |> 
    eBayes(robust = T)
getTT <- function(contrast) {
    topTable(limma.res, coef = contrast, p.value = 1, adjust.method = "none",
             n = Inf) |> select(P.Value) |> 
        mutate(Comparison = contrast) |> 
        rownames_to_column("PTM")
}

df1 <- lapply(contrasts_, getTT) |> bind_rows() |> filter(P.Value <= 0.05)
df2 <- df1 |> pivot_wider(id_cols = PTM, names_from = Comparison, values_from = P.Value)