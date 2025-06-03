# se object was generated as first chunk (make se) in Fig 1b.Rmd 

se_rel <- se[, !is.na(se$relapse) & se$relapse != "yes (>3 years)"]
se_rel$relapse <- ifelse(se_rel$relapse == "NO", "no", se_rel$relapse)
mm <- model.matrix(~0+relapse, data = colData(se_rel))
colnames(mm) <- gsub("relapse", "", colnames(mm))
contrasts.matrix <- makeContrasts(contrasts = "yes-no", levels = colnames(mm))
lmFit(assay(se_rel), mm) |>
    contrasts.fit(contrasts.matrix) |> 
    eBayes(robust = T) |> 
    topTable(coef = "yes-no", p.value = 1, n = Inf) |> 
    as_tibble(rownames = "PTMs") |> dplyr::select(1,2,5)