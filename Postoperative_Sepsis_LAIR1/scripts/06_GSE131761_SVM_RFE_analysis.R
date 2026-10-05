############################################################
# SVM-RFE Feature Selection Analysis
#
# Study:
# Integrated Transcriptomic and Functional Analyses Identify
# LAIR1 as a Key Immunoregulatory Molecule in Postoperative Sepsis
#
# This script performs SVM-RFE feature selection using
# five-fold cross-validation.
############################################################


############################################################
# 1. Load required package and helper functions
############################################################

library(e1071)

source("msvmRFE.R")


############################################################
# 2. Define input and output paths
############################################################

input_file <- "data/k.csv"

output_dir <- "results"

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


############################################################
# 3. Load input data
############################################################

input <- read.csv(
  input_file,
  row.names = 1,
  check.names = FALSE
)


############################################################
# 4. Define class labels
############################################################

input$lasso <- factor(
  input$lasso,
  levels = c(0, 1),
  labels = c("non_sepsis", "sepsis")
)


############################################################
# 5. Set random seed
############################################################

set.seed(2026)


############################################################
# 6. Generate five cross-validation folds
############################################################

nfold <- 5

fold_id <- sample(
  rep(
    seq_len(nfold),
    length.out = nrow(input)
  )
)

folds <- lapply(
  seq_len(nfold),
  function(i) which(fold_id == i)
)


############################################################
# 7. Perform SVM-RFE
############################################################

results <- lapply(
  folds,
  svmRFE.wrap,
  input,
  k = 5,
  halve.above = 100
)


############################################################
# 8. Rank selected features
############################################################

top.features <- WriteFeatures(
  results,
  input,
  save = FALSE
)


############################################################
# 9. Evaluate the top 1-15 ranked features
############################################################

featsweep <- lapply(
  seq_len(15),
  FeatSweep.wrap,
  results,
  input
)

errors <- vapply(
  featsweep,
  function(x) x$error,
  numeric(1)
)

best_n <- which.min(
  errors
)


############################################################
# 10. Export feature ranking
############################################################

write.csv(
  top.features,
  file.path(
    output_dir,
    "feature_svm_corrected.csv"
  ),
  row.names = FALSE
)


############################################################
# 11. Export optimal features
############################################################

write.csv(
  head(
    top.features$FeatureName,
    best_n
  ),
  file.path(
    output_dir,
    "top_corrected.csv"
  ),
  row.names = FALSE
)


############################################################
# 12. Calculate baseline classification error
############################################################

baseline_error <- min(
  prop.table(
    table(input$lasso)
  )
)


############################################################
# 13. Plot cross-validation error
############################################################

PlotErrors(
  errors,
  no.info = baseline_error
)


############################################################
# 14. Plot cross-validation accuracy
############################################################

Plotaccuracy(
  1 - errors,
  no.info = 1 - baseline_error
)
