############################################################
# Boxplot Analysis of Candidate Genes
#
# Study:
# Integrated Transcriptomic and Functional Analyses Identify
# LAIR1 as a Key Immunoregulatory Molecule in Postoperative Sepsis
#
# This script compares the expression levels of four
# candidate genes between SIRS and sepsis groups.
############################################################


############################################################
# 1. Load required packages
############################################################

library(ggpubr)
library(ggplot2)
library(cowplot)


############################################################
# 2. Define input paths
############################################################

expression_file <- "data/gene_boxplot.csv"

group_file <- "data/group.csv"


############################################################
# 3. Load expression data
############################################################

Exp <- read.csv(
  expression_file,
  header = T
)

class(Exp)

row.names(Exp) <- make.names(
  Exp[, 1],
  TRUE
)

Exp <- Exp[, -1]

Exp = t(Exp)

Exp = as.data.frame(Exp)


############################################################
# 4. Select candidate genes
############################################################

gene <- c(
  "LAIR1",
  "PLAC8",
  "ACSS2",
  "PCOLCE2"
)

gene <- as.vector(gene)

Exp_plot <- Exp[, gene]


############################################################
# 5. Load sample group information
############################################################

info <- read.csv(
  group_file,
  header = T
)

Exp_plot <- Exp_plot[
  info$Sample,
]

Exp_plot$sam = info$Type

Exp_plot$sam <- factor(
  Exp_plot$sam,
  levels = c(
    "non_sepsis",
    "sepsis"
  )
)


############################################################
# 6. Define plot colors
############################################################

col <- c(
  "#337AB7",
  "#F0AD4E"
)


############################################################
# 7. Generate boxplots for the four candidate genes
############################################################

plist2 <- list()

for (i in 1:length(gene)) {
  
  bar_tmp <- Exp_plot[
    ,
    c(
      gene[i],
      "sam"
    )
  ]
  
  colnames(bar_tmp) <- c(
    "Expression",
    "sam"
  )
  
  my_comparisons1 <- list(
    c(
      "SIRS",
      "sepsis"
    )
  )
  
  pb1 <- ggboxplot(
    bar_tmp,
    x = "sam",
    y = "Expression",
    color = "sam",
    fill = NULL,
    add = "jitter",
    bxp.errorbar.width = 0.6,
    width = 0.6,
    size = 0.01,
    font.label = list(
      size = 30
    ),
    palette = col
  ) +
    theme(
      panel.background = element_blank()
    )
  
  pb1 <- pb1 +
    theme(
      axis.line = element_line(
        colour = "black"
      )
    ) +
    theme(
      axis.title.x = element_blank()
    )
  
  pb1 <- pb1 +
    theme(
      axis.title.y = element_blank()
    ) +
    theme(
      axis.text.x = element_text(
        size = 15,
        angle = 45,
        vjust = 1,
        hjust = 1
      )
    )
  
  pb1 <- pb1 +
    theme(
      axis.text.y = element_text(
        size = 15
      )
    ) +
    ggtitle(
      gene[i]
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 15,
        face = "bold"
      )
    )
  
  pb1 <- pb1 +
    theme(
      legend.position = "NA"
    )
  
  pb1 <- pb1 +
    stat_compare_means(
      method = "t.test",
      hide.ns = F,
      comparisons = my_comparisons1,
      label = "p.signif"
    )
  
  plist2[[i]] <- pb1
}


############################################################
# 8. Combine the four boxplots
############################################################

plot_grid(
  plist2[[1]],
  plist2[[2]],
  plist2[[3]],
  plist2[[4]],
  ncol = 2
)
