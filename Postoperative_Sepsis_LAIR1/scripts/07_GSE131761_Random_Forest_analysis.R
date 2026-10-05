############################################################
# Random Forest Feature Selection Analysis
#
# Study:
# Integrated Transcriptomic and Functional Analyses Identify
# LAIR1 as a Key Immunoregulatory Molecule in Postoperative Sepsis
#
# This script performs:
# 1. Random forest classification
# 2. Optimization of the number of trees
# 3. Calculation of gene importance
# 4. Selection and visualization of the top-ranked genes
############################################################


############################################################
# 1. Load required packages
############################################################

library(randomForest)
library(ggpubr)


############################################################
# 2. Define input and output paths
############################################################

expression_file <- "data/gse131761_GENE.csv"

group_file <- "data/group.csv"

output_dir <- "results"

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


############################################################
# 3. Set random seed
############################################################

set.seed(123)


############################################################
# 4. Load gene expression and group information
############################################################

# Gene expression matrix

data <- read.csv(
  expression_file,
  row.names = 1
)

# Sample classification information

genelist <- read.csv(
  group_file,
  row.names = 1
)


############################################################
# 5. Prepare data for Random Forest analysis
############################################################

data = t(data)

# Use the "lasso" column as the sample group label

group <- as.character(
  genelist$lasso
)


############################################################
# 6. Construct the initial Random Forest model
############################################################

rf = randomForest(
  as.factor(group) ~ .,
  data = data,
  ntree = 500
)


############################################################
# 7. Plot Random Forest classification error
############################################################

pdf(
  file = file.path(
    output_dir,
    "Random_Forest.pdf"
  ),
  width = 6,
  height = 6
)

plot(
  rf,
  main = "Random forest",
  lwd = 2
)

dev.off()


############################################################
# 8. Determine the optimal number of trees
############################################################

optionTrees = which.min(
  rf$err.rate[, 1]
)

optionTrees


############################################################
# 9. Reconstruct Random Forest using optimal tree number
############################################################

rf2 = randomForest(
  as.factor(group) ~ .,
  data = data,
  ntree = optionTrees
)

plot(
  rf2,
  main = "Random forest",
  lwd = 2
)


############################################################
# 10. Extract gene importance
############################################################

importance = importance(
  x = rf2
)

importance = as.data.frame(
  importance
)

importance$size = rownames(
  importance
)

importance = importance[
  ,
  c(2, 1)
]

names(importance) = c(
  "Gene",
  "importance"
)


############################################################
# 11. Select the top 10 genes by importance
############################################################

af = importance[
  order(
    importance$importance,
    decreasing = T
  ),
]

af = af[
  1:10,
]


############################################################
# 12. Bar plot of the top 10 genes
############################################################

p1 = ggplot(
  af,
  aes(
    x = reorder(
      Gene,
      importance
    ),
    y = importance,
    fill = importance
  )
) +
  geom_bar(
    stat = "identity",
    width = 0.7
  ) +
  coord_flip() +
  scale_fill_gradient(
    low = ggsci::pal_npg()(2)[1],
    high = ggsci::pal_npg()(2)[2]
  ) +
  labs(
    x = "Gene",
    y = "Importance",
    title = "Top 10 Genes by Importance"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 1
    ),
    axis.text.y = element_text(
      size = 12
    ),
    plot.title = element_text(
      hjust = 0.5
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


############################################################
# 13. Save top-gene importance bar plot
############################################################

pdf(
  file = file.path(
    output_dir,
    "Random_Forest_importance_1.pdf"
  ),
  width = 6,
  height = 6
)

print(p1)

dev.off()


############################################################
# 14. Dot plot of the top 10 genes
############################################################

p2 = ggdotchart(
  af,
  x = "Gene",
  y = "importance",
  color = "importance",
  sorting = "descending",
  add = "segments",
  add.params = list(
    color = "lightgray",
    size = 2
  ),
  dot.size = 6,
  font.label = list(
    color = "white",
    size = 9,
    vjust = 0.5
  ),
  ggtheme = theme_bw(),
  rotate = TRUE
)

p3 = p2 +
  geom_hline(
    yintercept = 0,
    linetype = 2,
    color = "lightgray"
  ) +
  gradient_color(
    palette = c(
      ggsci::pal_npg()(2)[2],
      ggsci::pal_npg()(2)[1]
    )
  )

p3


############################################################
# 15. Save top-gene importance dot plot
############################################################

pdf(
  file = file.path(
    output_dir,
    "Random_Forest_importance_2.pdf"
  ),
  width = 6,
  height = 6
)

print(p3)


############################################################
# 16. Export ranked Random Forest feature genes
############################################################

rfGenes = importance[
  order(
    importance[, "importance"],
    decreasing = TRUE
  ),
]

write.csv(
  rfGenes,
  file.path(
    output_dir,
    "Random_Forest_feature_genes.csv"
  )
)