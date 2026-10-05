############################################################
# GSE74224 Data Preprocessing
#
# Study:
# Integrated Transcriptomic and Functional Analyses Identify
# LAIR1 as a Key Immunoregulatory Molecule in Postoperative Sepsis
#
# This script performs probe annotation and gene-level
# expression matrix generation for GSE74224.
############################################################


############################################################
# 1. Define input and output paths
############################################################

mat = "data/GSE74224_series_matrix.txt"

gpl = "data/GPL5175-3188.txt"

output_file = "results/GSE74224_gene_expression.txt"


############################################################
# 2. Load required packages
############################################################

library("limma")
library("tidyverse")


############################################################
# 3. Load expression matrix and GPL annotation
############################################################

exp = read.table(
  mat,
  header = T,
  sep = "\t",
  dec = ".",
  comment.char = "!",
  na.strings = c("NA"),
  fill = T
)

GPL_file = read.table(
  gpl,
  header = T,
  quote = "",
  sep = "\t",
  dec = ".",
  comment.char = "#",
  na.strings = c("NA"),
  fill = T
)


############################################################
# 4. Extract gene annotation information
############################################################

separate(
  data = GPL_file,
  col = gene_assignment,
  into = c("a", "b", "c", "d", "e"),
  sep = "//"
)

GPL_file[c("a", "Symbol", "c")] =
  str_split_fixed(
    GPL_file$gene_assignment,
    "//",
    3
  )


############################################################
# 5. Match probe IDs to gene symbols
############################################################

gpl_file = GPL_file[, c(1, 4)]

exp = as.data.frame(exp)

colnames(exp)[1] = "ID"

exp_symbol = merge(
  exp,
  gpl_file,
  by = "ID"
)


############################################################
# 6. Remove probes without valid gene symbols
############################################################

exp_symbol[exp_symbol == ""] = NA

exp_symbol = na.omit(exp_symbol)

exp_symbol[
  ,
  grep(
    "Symbol",
    colnames(exp_symbol)
  )
] = trimws(
  exp_symbol[
    ,
    grep(
      "Symbol",
      colnames(exp_symbol)
    )
  ]
)


############################################################
# 7. Check duplicated gene symbols
############################################################

table(
  duplicated(
    exp_symbol[
      ,
      ncol(exp_symbol)
    ]
  )
)

d = data.frame(
  duplicated(
    exp_symbol[
      ,
      ncol(exp_symbol)
    ]
  )
)


############################################################
# 8. Average multiple probes mapped to the same gene
############################################################

exp_symbol = avereps(
  exp_symbol[
    ,
    -c(
      1,
      ncol(exp_symbol)
    )
  ],
  ID = exp_symbol$Symbol
)

table(
  duplicated(
    rownames(exp_symbol)
  )
)


############################################################
# 9. Export processed expression matrix
############################################################

write.table(
  exp_symbol,
  file = output_file.txt,
  sep = "\t",
  quote = F,
  col.names = T
)

