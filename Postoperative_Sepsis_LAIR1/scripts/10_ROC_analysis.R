############################################################
# ROC Analysis of Candidate Genes and Four-Gene Model
#
# Study:
# Integrated Transcriptomic and Functional Analyses Identify
# LAIR1 as a Key Immunoregulatory Molecule in Postoperative Sepsis
#
# This script evaluates the diagnostic performance of
# LAIR1, PLAC8, ACSS2, and PCOLCE2 individually and as a
# combined four-gene logistic regression model.
############################################################


############################################################
# 1. Load required package
############################################################

library(pROC)


############################################################
# 2. Define input path
############################################################

input_file <- "data/gse131761_GENE_noncontrol.csv"


############################################################
# 4. Load data
############################################################

roc_data = read.csv(
  input_file,
  header = T
)

rownames(roc_data) = make.names(
  roc_data[, 1]
)

roc_data = roc_data[, -1]


############################################################
# 5. ROC analysis for LAIR1
############################################################

roc1 = plot.roc(
  roc_data$type,
  roc_data$`LAIR1`,
  col = "#A31621",
  percent = TRUE,
  lwd = 2,
  print.auc = TRUE,
  print.auc.cex = 1,
  print.auc.pattern = "AUC: %.1f%%",
  print.auc.y = 20
)

ci.auc(roc1)


############################################################
# 6. ROC analysis for PLAC8
############################################################

roc1 = plot.roc(
  roc_data$type,
  roc_data$`PLAC8`,
  col = "#A31621",
  percent = TRUE,
  lwd = 2,
  print.auc = TRUE,
  print.auc.cex = 1,
  print.auc.pattern = "AUC: %.1f%%",
  print.auc.y = 20
)

ci.auc(roc1)


############################################################
# 7. ROC analysis for ACSS2
############################################################

roc1 = plot.roc(
  roc_data$type,
  roc_data$`ACSS2`,
  col = "#A31621",
  percent = TRUE,
  lwd = 2,
  print.auc = TRUE,
  print.auc.cex = 1,
  print.auc.pattern = "AUC: %.1f%%",
  print.auc.y = 20
)

ci.auc(roc1)


############################################################
# 8. ROC analysis for PCOLCE2
############################################################

roc1 = plot.roc(
  roc_data$type,
  roc_data$`PCOLCE2`,
  col = "#A31621",
  percent = TRUE,
  lwd = 2,
  print.auc = TRUE,
  print.auc.cex = 1,
  print.auc.pattern = "AUC: %.1f%%",
  print.auc.y = 20
)

ci.auc(roc1)


############################################################
# 9. Construct four-gene logistic regression model
############################################################

model_1 <- glm(
  type ~ LAIR1 + PLAC8 + ACSS2 + PCOLCE2,
  data = roc_data,
  family = binomial(
    link = "logit"
  )
)

summary(model_1)


############################################################
# 10. Generate predicted probabilities
############################################################

fitted.prob <- predict(
  model_1,
  newdata = roc_data,
  type = "response"
)

roc_data$pred <- model_1$fitted.values


############################################################
# 11. ROC analysis of the four-gene model
############################################################

roc_multivar_1 <- roc(
  roc_data$type,
  roc_data[, "pred"]
)

plot.roc(
  roc_multivar_1,
  col = "red"
)


############################################################
# 12. Plot ROC curve and calculate confidence interval
#     for the four-gene model
############################################################

roc1 = plot.roc(
  roc_data$type,
  roc_data$pred,
  col = "#A31621",
  percent = TRUE,
  lwd = 2,
  print.auc = TRUE,
  print.auc.cex = 1,
  print.auc.pattern = "AUC: %.1f%%",
  print.auc.y = 20
)

ci.auc(roc1)
