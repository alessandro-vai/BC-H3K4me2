library(tidyverse) |> suppressPackageStartupMessages()

annot_table <- read.delim("annot_table.txt") |> 
    mutate(
        Feature = Genomic_feature |> as_factor() |> fct_inorder() |> fct_rev(),
        Legend = factor(Legend, levels = c("TN only", "LuA only", "Common")),
        .by = Legend
    )

pastel_palette <- c(
    # Reds
    "#FA8072",
    "#FFDAB9",
    "#FFB6C1",
    
    # Greens
    "#90EE90",
    "#FAFAD2",
    "#66CDAA",
    "#AFEEEE",
    
    # Blues
    "#87CEEB",
    "#F0F8FF",
    
    # Additional color
    "#D8BFD8"
)

for (leg in c("TN only", "LuA only", "Common")) {
    p <- ggplot(annot_table |> filter(Legend == leg),
                aes(x = perc, y = factor(1), fill=Feature)) +
        geom_bar(stat = "identity", width = 0.5) +
        scale_fill_manual(values = rev(pastel_palette)) +
        theme_void() +
        ggtitle(leg) +
        theme(legend.position = "none",
              title = element_text(family = "ArialMT", size = 14),
              plot.title = element_text(hjust = 0.5, vjust = -5))
}

