library(tidyverse) |> suppressPackageStartupMessages()
library(extrafont)
library(FactoMineR)
library(factoextra)

# Samples were filtered for 70% of valid values,
# ratios were mean subtracted and log transformed,
# Missing values were replaced with 0

# Each row correspond to a histone PTM and each column
# to a sample, called like Normal 1, LuA 1 etc.

ratios <- read.delim('ratios for PCA.txt',
                     row.names = 1) |> t() |> 
    as_tibble(rownames = 'Tissue') |> 
    mutate(Tissue = str_remove(Tissue, '\\..*'),
           Cat = if_else(grepl('Normal', Tissue), 'Normal', 'Tumor'))
df.pca <- PCA(ratios |> select(-c(Tissue, Cat)),
              graph = FALSE, scale.unit = F, ncp = 2)

points_data <- df.pca$ind$coord |> as_tibble() |>  
    mutate(Tissue = ratios$Tissue,
           Cat = ratios$Cat)

scale_var_data <- function(PCAobj, select_var=10){
    # from https://github.com/kassambara/factoextra/blob/master/R/fviz_pca.R
    
    # Data frame to be used for plotting
    var <- facto_summarize(PCAobj, element = "var", 
                           result = c("coord", "contrib", "cos2"),
                           axes = c(1,2)) # PC1 ad PC2
    colnames(var)[2:3] <-  c("x", "y")
    
    pca.ind <- get_pca_ind(PCAobj)
    ind <- data.frame(pca.ind$coord[, c(1,2), drop=FALSE],
                      stringsAsFactors = TRUE)
    colnames(ind)<- c("x", "y")
    
    # rescale variable coordinates
    scale_factor <- min(
        (max(ind[,"x"])-min(ind[,"x"])/(max(var[,"x"])-min(var[,"x"]))),
        (max(ind[,"y"])-min(ind[,"y"])/(max(var[,"y"])-min(var[,"y"])))
    ) * 0.7
    var[, c("x", "y")] <- var[, c("x", "y")]*scale_factor
    var <- var |> arrange(desc(contrib)) |> slice_head(n = select_var)
    return(var)
}

ggplot(points_data, aes(x = Dim.1, y = Dim.2)) +
    geom_point(aes(shape=Cat, color = Tissue), size = 3.5) + # size = 4
    scale_shape_manual(values = c(8, 16)) + theme_classic() +
    scale_color_manual(values = c('#ff80c0', # Her
                                  '#80ffff', # LuA
                                  '#0080ff', # LuB
                                  'black',
                                  'black')) +
    xlab(paste0('Component 1 (',
                round(df.pca$eig['comp 1', 'percentage of variance'], 1),
                '%)')) +
    ylab(paste0('Component 2 (',
                round(df.pca$eig['comp 2', 'percentage of variance'], 1),
                '%)')) + 
    theme(#text = element_text(family = 'ArialMT'),
        legend.position="none") +
    geom_segment(data = scale_var_data(df.pca) |> 
                     mutate(xstart = 0, ystart = 0), 
                 aes(x = xstart, y = ystart,
                     xend = x, yend = y),
                 arrow = grid::arrow(length = grid::unit(0.2, 'cm')),
                 color = 'black', alpha = .8, linewidth = .5) +
    ggrepel::geom_text_repel(data = scale_var_data(df.pca),
                             aes(x = x, y = y, label = name))


