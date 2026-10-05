############################################################
# Nomogram, Calibration, Decision Curve and Clinical Impact
# Curve Analyses
#
# Study:
# Integrated Transcriptomic and Functional Analyses Identify
# LAIR1 as a Key Immunoregulatory Molecule in Postoperative Sepsis
#
# This script performs:
# 1. Construction of the four-gene logistic regression model
# 2. Standard nomogram construction
# 3. Calculation of individual nomogram scores
# 4. Graphical nomogram using regplot
# 5. Bootstrap calibration analysis
# 6. Decision curve analysis (DCA)
# 7. Clinical impact curve (CIC) of the four-gene model
# 8. Clinical impact curves for individual diagnostic genes
#
# Outcome coding:
# type = 0: control
# type = 1: sepsis
############################################################


############################################################
# 1. Load packages
############################################################

library(rms)
library(rmda)
library(nomogramFormula)
library(regplot)


############################################################
# 2. Define input and output paths
############################################################

input_file <- "data/GSE131761_model_input.csv"

output_dir <- "results"

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


############################################################
# 3. Load and prepare data
############################################################

mydata <- read.csv(
  input_file,
  row.names = 1,
  check.names = FALSE
)

# Variables required for the diagnostic model
required_variables <- c(
  "LAIR1",
  "PLAC8",
  "ACSS2",
  "PCOLCE2",
  "type"
)

# Check whether all required variables are available
stopifnot(
  all(required_variables %in% colnames(mydata))
)

# Retain only variables used in the present analysis
mydata <- mydata[
  ,
  required_variables,
  drop = FALSE
]

# Remove samples with missing values
mydata <- na.omit(mydata)

# Outcome coding:
# 0 = control
# 1 = sepsis
mydata$type <- as.numeric(
  as.character(mydata$type)
)

# Confirm correct binary coding
stopifnot(
  all(mydata$type %in% c(0, 1))
)

# Display group distribution
print(
  table(mydata$type)
)

# Keep a clean copy for DCA/CIC analyses
dca_data <- mydata


############################################################
# 4. Construct the four-gene logistic regression model
############################################################

dd <- datadist(mydata)

old_options <- options(
  datadist = "dd"
)

formula_4gene <- type ~
  LAIR1 +
  PLAC8 +
  ACSS2 +
  PCOLCE2

model_lrm <- lrm(
  formula_4gene,
  data = mydata,
  x = TRUE,
  y = TRUE
)

print(model_lrm)


############################################################
# 5. Construct standard nomogram
############################################################

nomogram_4gene <- nomogram(
  model_lrm,
  fun = function(x) {
    1 / (1 + exp(-x))
  },
  lp = FALSE,
  fun.at = c(
    0.1,
    0.3,
    0.5,
    0.7,
    0.9
  ),
  funlabel = "Risk"
)

pdf(
  file.path(
    output_dir,
    "Four_gene_nomogram_standard.pdf"
  ),
  width = 12,
  height = 7
)

plot(
  nomogram_4gene,
  xfrac = 0.35,
  cex.var = 1.6,
  cex.axis = 1.4,
  tcl = -0.5,
  lmgp = 0.3,
  label.every = 1,
  col.grid = gray(
    c(0.8, 0.95)
  )
)

dev.off()


############################################################
# 6. Calculate nomogram scores for all samples
############################################################

nomogram_formula <- formula_rd(
  nomogram = nomogram_4gene
)

mydata$Nomogram_score <- points_cal(
  formula = nomogram_formula$formula,
  rd = mydata
)

# Display first several scores
print(
  head(
    mydata[
      ,
      c(
        "LAIR1",
        "PLAC8",
        "ACSS2",
        "PCOLCE2",
        "type",
        "Nomogram_score"
      )
    ]
  )
)

# Save nomogram scores
write.csv(
  mydata,
  file = file.path(
    output_dir,
    "GSE131761_nomogram_scores.csv"
  ),
  row.names = TRUE
)


############################################################
# 7. Graphical nomogram using regplot
############################################################

# The second sample is displayed here to remain consistent
# with the original analysis workflow.
#
# If another sample was used in the final published figure,
# replace 2 with the corresponding row number or sample ID.

selected_observation <- mydata[
  2,
  ,
  drop = FALSE
]

pdf(
  file.path(
    output_dir,
    "Four_gene_nomogram_regplot.pdf"
  ),
  width = 16,
  height = 6
)

regplot(
  model_lrm,
  plots = c(
    "violin",
    "boxes"
  ),
  observation = selected_observation,
  center = TRUE,
  subticks = TRUE,
  droplines = TRUE,
  title = "Nomogram",
  points = TRUE,
  odds = TRUE,
  showP = TRUE,
  rank = "sd",
  interval = "confidence",
  clickable = FALSE
)

dev.off()


############################################################
# 8. Bootstrap calibration analysis
############################################################

set.seed(123)

calibration <- calibrate(
  model_lrm,
  method = "boot",
  B = 1000
)

pdf(
  file.path(
    output_dir,
    "Four_gene_calibration_curve.pdf"
  ),
  width = 7,
  height = 7
)

plot(
  calibration,
  xlim = c(0, 1),
  ylim = c(0, 1),
  xlab = "Predicted Probability",
  ylab = "Observed Probability",
  subtitles = FALSE
)

dev.off()


############################################################
# 9. Decision curve analysis of the four-gene model
############################################################

dca_4gene <- decision_curve(
  formula_4gene,
  data = dca_data,
  family = binomial(
    link = "logit"
  ),
  thresholds = seq(
    0,
    1,
    by = 0.01
  ),
  confidence.intervals = 0.95,
  study.design = "case-control",
  population.prevalence = 0.30
)


############################################################
# 10. Plot four-gene decision curve
############################################################

pdf(
  file.path(
    output_dir,
    "Four_gene_DCA.pdf"
  ),
  width = 7,
  height = 7
)

plot_decision_curve(
  dca_4gene,
  curve.names = "Four-gene model",
  xlim = c(0, 1),
  cost.benefit.axis = FALSE,
  col = "red",
  confidence.intervals = FALSE,
  standardize = FALSE
)

dev.off()


############################################################
# 11. Clinical impact curve of the four-gene model
############################################################

pdf(
  file.path(
    output_dir,
    "Four_gene_CIC.pdf"
  ),
  width = 7,
  height = 7
)

plot_clinical_impact(
  dca_4gene,
  population.size = 1000,
  cost.benefit.axis = TRUE,
  n.cost.benefits = 8,
  col = c(
    "red",
    "blue"
  ),
  confidence.intervals = TRUE
)

dev.off()


############################################################
# 12. Define individual-gene models
############################################################

formula_LAIR1 <- type ~ LAIR1

formula_PLAC8 <- type ~ PLAC8

formula_ACSS2 <- type ~ ACSS2

formula_PCOLCE2 <- type ~ PCOLCE2


############################################################
# 13. Generate DCA objects for individual genes
############################################################

dca_LAIR1 <- decision_curve(
  formula_LAIR1,
  data = dca_data,
  family = binomial(
    link = "logit"
  ),
  thresholds = seq(
    0,
    1,
    by = 0.01
  ),
  confidence.intervals = 0.95,
  study.design = "case-control",
  population.prevalence = 0.30
)


dca_PLAC8 <- decision_curve(
  formula_PLAC8,
  data = dca_data,
  family = binomial(
    link = "logit"
  ),
  thresholds = seq(
    0,
    1,
    by = 0.01
  ),
  confidence.intervals = 0.95,
  study.design = "case-control",
  population.prevalence = 0.30
)


dca_ACSS2 <- decision_curve(
  formula_ACSS2,
  data = dca_data,
  family = binomial(
    link = "logit"
  ),
  thresholds = seq(
    0,
    1,
    by = 0.01
  ),
  confidence.intervals = 0.95,
  study.design = "case-control",
  population.prevalence = 0.30
)


dca_PCOLCE2 <- decision_curve(
  formula_PCOLCE2,
  data = dca_data,
  family = binomial(
    link = "logit"
  ),
  thresholds = seq(
    0,
    1,
    by = 0.01
  ),
  confidence.intervals = 0.95,
  study.design = "case-control",
  population.prevalence = 0.30
)


############################################################
# 14. Plot individual-gene clinical impact curves
#     Corresponding to Supplementary Figure 3
############################################################

pdf(
  file.path(
    output_dir,
    "Supplementary_Figure_3_individual_gene_CIC.pdf"
  ),
  width = 12,
  height = 10
)

par(
  mfrow = c(2, 2),
  mar = c(
    5,
    5,
    4,
    2
  )
)


############################################################
# Panel A: LAIR1
############################################################

plot_clinical_impact(
  dca_LAIR1,
  population.size = 1000,
  cost.benefit.axis = TRUE,
  n.cost.benefits = 8,
  col = c(
    "red",
    "blue"
  ),
  confidence.intervals = TRUE
)

title(
  main = "LAIR1",
  cex.main = 1.4
)

mtext(
  "A",
  side = 3,
  adj = -0.12,
  line = 1,
  cex = 1.5
)


############################################################
# Panel B: PLAC8
############################################################

plot_clinical_impact(
  dca_PLAC8,
  population.size = 1000,
  cost.benefit.axis = TRUE,
  n.cost.benefits = 8,
  col = c(
    "red",
    "blue"
  ),
  confidence.intervals = TRUE
)

title(
  main = "PLAC8",
  cex.main = 1.4
)

mtext(
  "B",
  side = 3,
  adj = -0.12,
  line = 1,
  cex = 1.5
)


############################################################
# Panel C: ACSS2
############################################################

plot_clinical_impact(
  dca_ACSS2,
  population.size = 1000,
  cost.benefit.axis = TRUE,
  n.cost.benefits = 8,
  col = c(
    "red",
    "blue"
  ),
  confidence.intervals = TRUE
)

title(
  main = "ACSS2",
  cex.main = 1.4
)

mtext(
  "C",
  side = 3,
  adj = -0.12,
  line = 1,
  cex = 1.5
)


############################################################
# Panel D: PCOLCE2
############################################################

plot_clinical_impact(
  dca_PCOLCE2,
  population.size = 1000,
  cost.benefit.axis = TRUE,
  n.cost.benefits = 8,
  col = c(
    "red",
    "blue"
  ),
  confidence.intervals = TRUE
)

title(
  main = "PCOLCE2",
  cex.main = 1.4
)

mtext(
  "D",
  side = 3,
  adj = -0.12,
  line = 1,
  cex = 1.5
)

dev.off()
