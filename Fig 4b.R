library(edgeR)
library(tidyverse) |> suppressPackageStartupMessages()

count_table <- read.delim("promoters_count_table.txt", comment.char = "#")
count_table <- read.delim("enhancers_count_table.txt", comment.char = "#")
# Same code starting from different count tables
count_table <- count_table[, ! colnames(count_table) %in% c("Chr", "Start", "End", "Length", "Strand")] %>% 
    select(-starts_with("INPUT")) %>%
    tibble::column_to_rownames("Geneid") %>% as.matrix()

colnames(count_table) <- str_remove(colnames(count_table), ".sorted.bam")
y <- DGEList(count_table, group = ifelse(grepl("LuA", colnames(count_table)), "LuA", "TN"))
keep.exprs <- filterByExpr(y)
y <- y[keep.exprs,, keep.lib.sizes=FALSE]
logcpm <- log2(cpm(y) + 1)

temp <- plotMDS(logcpm, plot=F)
plotdata <- data.frame(x = temp$x, y = temp$y, subtype = y$samples$group)

plotdata |> ggplot(aes(x = x, y = y, color = subtype)) +
    geom_point(size=4) + theme_bw() + xlab(paste0(temp$axislabel, "1 (", round(temp$var.explained[1] * 100), "%)")) +
    ylab(paste0(temp$axislabel, "2 (",round(temp$var.explained[2] * 100), "%)")) +
    scale_color_manual(values=c("#80ffff", "black"))