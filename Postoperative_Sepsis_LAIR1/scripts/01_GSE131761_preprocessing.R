############################################################
# GSE131761 Data Preprocessing
#
# Study:
# Integrated Transcriptomic and Functional Analyses Identify
# LAIR1 as a Key Immunoregulatory Molecule in Postoperative Sepsis
#
# This script performs probe-to-gene annotation and generates
# a gene-level expression matrix for GSE131761.
############################################################


############################################################
# 1. Load required package
############################################################

library(limma)


############################################################
# 2. Define input and output paths
############################################################

expression_file <- "data/GSE131761_series_matrix.txt"

platform_file <- "data/GPL13497-9755.txt"

output_file <- "results/GSE131761_gene.csv"


############################################################
# 3. Load expression matrix
############################################################

# Rows = probes
# Columns = samples
# The first column contains probe IDs and is used as row names.

expr <- read.table(
  expression_file,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  comment.char = "!",
  check.names = FALSE
)


############################################################
# 4. Load platform annotation
############################################################

gpl <- read.table(
  platform_file,
  header = TRUE,
  sep = "\t",
  comment.char = "#",
  quote = "",
  fill = TRUE
)


############################################################
# 5. Check required annotation columns
############################################################

stopifnot(
  all(
    c("ID", "GENE_SYMBOL") %in% colnames(gpl)
  )
)


############################################################
# 6. Convert expression data to numeric matrix
############################################################

expr <- as.matrix(expr)

stopifnot(
  is.numeric(expr)
)


############################################################
# 7. Match probe IDs to gene symbols
############################################################

# Gene symbols are matched according to probe IDs while
# preserving the original row order of the expression matrix.

symbol <- gpl$GENE_SYMBOL[
  match(
    rownames(expr),
    gpl$ID
  )
]

symbol <- trimws(
  as.character(symbol)
)


############################################################
# 8. Remove probes without valid gene symbols or expression
############################################################

keep <- !is.na(symbol) &
  nzchar(symbol) &
  rowSums(
    !is.finite(expr)
  ) == 0


############################################################
# 9. Average multiple probes mapped to the same gene
############################################################

expr_gene <- avereps(
  expr[
    keep,
    ,
    drop = FALSE
  ],
  ID = symbol[
    keep
  ]
)


############################################################
# 10. Prepare output table
############################################################

result <- data.frame(
  GENE_SYMBOL = rownames(expr_gene),
  expr_gene,
  check.names = FALSE
)


############################################################
# 11. Export gene-level expression matrix
############################################################

write.csv(
  result,
  output_file,
  row.names = FALSE
)


############################################################
# 12. Display final matrix dimensions
############################################################

dim(expr_gene)

