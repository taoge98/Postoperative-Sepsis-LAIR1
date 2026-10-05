############################################################
# LASSO Regression Feature Selection Analysis
#
# Study:
# Integrated Transcriptomic and Functional Analyses Identify
# LAIR1 as a Key Immunoregulatory Molecule in Postoperative Sepsis
#
# This script performs LASSO regression analysis for
# candidate gene feature selection using GSE131761.
############################################################


############################################################
# 1. Load required package
############################################################

library(glmnet)


############################################################
# 2. Define input and output paths
############################################################

expression_file <- "data/gse131761_GENE.csv"

candidate_gene_file <- "data/gene16.csv"

group_file <- "data/fen.csv"

output_dir <- "results"

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


############################################################
# 3. Load gene expression data
############################################################

k <- read.csv(
  expression_file
)

row.names(k) <- make.names(
  k[, 1],
  TRUE
)

k <- k[, -1]


############################################################
# 4. Load candidate genes and retain selected genes
############################################################

gene <- read.csv(
  candidate_gene_file
)

k <- k[
  rownames(k) %in% gene$gene,
]


############################################################
# 5. Prepare sample-level expression matrix
############################################################

k <- as.data.frame(
  t(k)
)

fen <- read.csv(
  group_file,
  row.names = 1
)

k <- k[
  rownames(fen),
]

k <- cbind(
  fen$lasso,
  k
)


############################################################
# 6. Prepare predictor matrix and outcome variable
############################################################

states <- as.matrix(k)

x <- states[, -1]

y <- states[, 1]


############################################################
# 7. Five-fold cross-validation for LASSO
############################################################

cvfit = cv.glmnet(
  x,
  y,
  family = "binomial",
  type.measure = "mse",
  nfolds = 5,
  alpha = 1
)

plot(cvfit)

cvfit$lambda.min

c(
  cvfit$lambda.min,
  cvfit$lambda.1se
)


############################################################
# 8. Fit LASSO regression model
############################################################

lasso <- glmnet(
  x,
  y,
  family = "binomial",
  alpha = 1,
  nlambda = 100
)


############################################################
# 9. Extract coefficients at specified lambda values
############################################################

coef(
  lasso,
  s = c(
    0.01732064,
    0.05805191
  )
)


############################################################
# 10. Visualize LASSO coefficient paths
############################################################

win.graph(
  width = 6,
  height = 5,
  pointsize = 9
)

print(lasso)

plot(
  lasso,
  label = T
)

plot(
  lasso,
  xvar = "lambda",
  label = T
)


############################################################
# 11. Extract regression coefficients
############################################################

lasso.coef <- predict(
  lasso,
  s = 0.4,
  type = "coefficients"
)


############################################################
# 12. Plot explained deviance
############################################################

plot(
  lasso,
  xvar = "dev",
  label = T
)


############################################################
# 13. Generate fitted probabilities
############################################################

lasso.y <- predict(
  lasso,
  newx = x,
  type = "response",
  s = 0.4
)