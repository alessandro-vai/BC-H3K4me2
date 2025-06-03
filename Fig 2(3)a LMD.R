library(limma)

ratios <- readxl::read_xlsx("Supplementary Dataset 2.xlsx",
                            skip = 1, sheet = "LMD samples", col_names = F) |>
    tibble::column_to_rownames("...1")
group <- c(rep("LuA", 4), rep("TNnoRel", 6), rep("TNRel", 6))
mm0 <- model.matrix(~0+group, data = ratios)
colnames(mm0) <- stringr::str_remove(colnames(mm0), "group")
contrast.matrix <- makeContrasts(contrasts = c("TNnoRel-LuA",
                                               "TNRel-LuA",
                                               "TNRel-TNnoRel",
                                               "(TNnoRel+TNRel)/2-LuA"),
                                 levels = colnames(mm0))
limma.res <- lmFit(ratios, mm0) |>
    contrasts.fit(contrast.matrix) |> 
    eBayes(robust = T)

getTopTable <- function(contrast){
    topTable(limma.res, adjust.method = "none",
             p.value = 1, n=Inf,
             coef = contrast) |> 
        dplyr::mutate(Contrast = contrast) |> 
        tibble::rownames_to_column("PTM")
}

lapply(colnames(contrast.matrix), getTopTable) |> 
    purrr::reduce(dplyr::bind_rows) |> 
    dplyr::mutate(Contrast = ifelse(Contrast == "(TNnoRel+TNRel)/2-LuA", "TN-LuA", Contrast)) |> 
    writexl::write_xlsx("stats TN LuA LMD.xlsx", format_headers = F)

group <- c(rep("LuA", 4), rep("TN", 12))
mm0 <- model.matrix(~group, data = ratios)
colnames(mm0) <- stringr::str_remove(colnames(mm0), "group")
limma.res <- lmFit(ratios, mm0) |>
    eBayes(robust = T)
getTopTable(2)
