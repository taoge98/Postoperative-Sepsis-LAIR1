############################################################
# WGCNA Analysis of GSE131761
#
# Study:
# Integrated Transcriptomic and Functional Analyses Identify
# LAIR1 as a Key Immunoregulatory Molecule in Postoperative Sepsis
#
# This script performs:
# 1. Gene variability filtering
# 2. Sample clustering and outlier removal
# 3. Soft-threshold power selection
# 4. Weighted gene co-expression network construction
# 5. Module-trait correlation analysis
# 6. Module membership and gene significance analysis
# 7. Identification of hub genes in the brown module
############################################################


############################################################
# 1. Clear workspace and load required package
############################################################

rm(list = ls())

library(WGCNA)

enableWGCNAThreads()


############################################################
# 2. Define input and output paths
############################################################

expression_file <- "data/tpm.csv"

trait_file <- "metadata/GSE131761_WGCNA_traits.csv"

output_dir <- "results"

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


############################################################
# 3. Load expression data
############################################################

data <- read.csv(
  expression_file,
  header = T
)

row.names(data) <- make.names(
  data[, 1],
  TRUE
)

data <- data[, -1]


############################################################
# 4. Gene variability filtering
############################################################

m.mad <- apply(data, 1, mad)

y = as.data.frame(m.mad)

dataExprVar <- data[
  which(
    m.mad >
      max(
        quantile(
          m.mad,
          probs = seq(0, 1, 0.25)
        )[2],
        0.5
      )
  ),
]

dim(dataExprVar)

dat = t(dataExprVar)


# Retain the top 25% most variable genes

data <- as.data.frame(
  t(data)
)

m.vars <- apply(
  data,
  2,
  var
)

expro.upper <- data[
  ,
  which(
    m.vars >
      quantile(
        m.vars,
        probs = seq(0, 1, 0.25)
      )[4]
  )
]

data <- as.data.frame(
  t(expro.upper)
)

dat = t(data)


############################################################
# 5. Check data quality
############################################################

gsg <- goodSamplesGenes(
  dat,
  verbose = 3
)

gsg$allOK


############################################################
# 6. Sample clustering and outlier detection
############################################################

sampleTree = hclust(
  dist(dat),
  method = "average"
)

sizeGrWindow(12, 9)

par(
  cex = 0.6
)

par(
  mar = c(0, 4, 2, 0)
)

plot(
  sampleTree,
  main = "Sample clustering to detectoutliers",
  sub = "",
  xlab = "",
  cex.lab = 1.5,
  cex.axis = 1.5,
  cex.main = 2
)

abline(
  h = 80,
  col = "red"
)


############################################################
# 7. Remove sample outliers
############################################################

clust = cutreeStatic(
  sampleTree,
  cutHeight = 80,
  minSize = 10
)

table(clust)

keepSamples = (clust == 1)

datExpr = dat[
  keepSamples,
]

nGenes = ncol(datExpr)

nSamples = nrow(datExpr)


############################################################
# 8. Load clinical trait data
############################################################

traitData = read.csv(
  trait_file,
  header = T
)

dim(traitData)

names(traitData)

femaleSamples = rownames(datExpr)

traitRows = match(
  femaleSamples,
  traitData$ID
)

datTraits = traitData[
  traitRows,
  -1
]

rownames(datTraits) = traitData[
  traitRows,
  1
]

collectGarbage()


############################################################
# 9. Visualize sample clustering and clinical traits
############################################################

sampleTree2 = hclust(
  dist(datExpr),
  method = "average"
)

traitColors = numbers2colors(
  datTraits,
  signed = FALSE
)

plotDendroAndColors(
  sampleTree2,
  traitColors,
  groupLabels = names(datTraits),
  main = "Sample dendrogram and trait heatmap"
)


############################################################
# 10. Select soft-thresholding power
############################################################

powers = c(
  c(1:10),
  seq(
    from = 12,
    to = 20,
    by = 2
  )
)

sft = pickSoftThreshold(
  datExpr,
  powerVector = powers,
  verbose = 5
)

sizeGrWindow(9, 5)

par(
  mfrow = c(1, 2)
)

cex1 = 0.85

plot(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) *
    sft$fitIndices[, 2],
  xlab = "SoftThreshold(power)",
  ylab = "ScaleFreeTopologyModelFit,signedR^2",
  type = "n",
  main = paste("Scaleindependence")
)

text(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) *
    sft$fitIndices[, 2],
  labels = powers,
  cex = cex1,
  col = "red"
)

abline(
  h = 0.9,
  col = "red"
)

plot(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  xlab = "SoftThreshold(power)",
  ylab = "MeanConnectivity",
  type = "n",
  main = paste("Meanconnectivity")
)

text(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  labels = powers,
  cex = cex1,
  col = "red"
)


############################################################
# 11. Construct weighted co-expression network
############################################################

net = blockwiseModules(
  datExpr,
  power = 12,
  TOMType = "unsigned",
  minModuleSize = 60,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = file.path(
    output_dir,
    "femaleMouseTOM"
  ),
  verbose = 3
)

table(net$colors)


############################################################
# 12. Visualize identified modules
############################################################

sizeGrWindow(12, 9)

mergedColors = labels2colors(
  net$colors
)

plotDendroAndColors(
  net$dendrograms[[1]],
  mergedColors[
    net$blockGenes[[1]]
  ],
  "Modulecolors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)


############################################################
# 13. Save network construction objects
############################################################

moduleLabels = net$colors

moduleColors = labels2colors(
  net$colors
)

MEs = net$MEs

geneTree = net$dendrograms[[1]]

save(
  MEs,
  moduleLabels,
  moduleColors,
  geneTree,
  file = file.path(
    output_dir,
    "FemaleLiver-02-networkConstruction-auto.RData"
  )
)


############################################################
# 14. Load clinical information for module-trait analysis
############################################################

data_clin <- read.csv(
  trait_file,
  header = T
)

rownames(data_clin) <- data_clin$ID

data_clin <- data_clin[, -1]


############################################################
# 15. Module-trait correlation analysis
############################################################

nGenes <- ncol(datExpr)

nSamples <- nrow(datExpr)

MEs0 <- moduleEigengenes(
  datExpr,
  moduleColors
)$eigengenes

MEs <- orderMEs(
  MEs0
)

moduleTraitCor <- cor(
  MEs,
  data_clin,
  use = "p"
)

moduleTraitPvalue <- corPvalueStudent(
  moduleTraitCor,
  nSamples
)


############################################################
# 16. Visualize module-trait relationships
############################################################

textMatrix <- paste(
  signif(
    moduleTraitCor,
    2
  ),
  "\n(",
  signif(
    moduleTraitPvalue,
    1
  ),
  ")",
  sep = ""
)

dim(textMatrix) <- dim(
  moduleTraitCor
)

par(
  mar = c(
    6,
    8.5,
    3,
    3
  )
)

labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(data_clin),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.5,
  zlim = c(-1, 1),
  main = paste(
    "Module-trait relationships"
  )
)

dev.off()


############################################################
# 17. Calculate module membership and gene significance
############################################################

nSamples <- nrow(datExpr)

modNames <- substring(
  names(MEs),
  3
)

geneModuleMembership <- as.data.frame(
  cor(
    datExpr,
    MEs,
    use = "p"
  )
)

MMPvalue <- as.data.frame(
  corPvalueStudent(
    as.matrix(
      geneModuleMembership
    ),
    nSamples
  )
)

names(geneModuleMembership) <- paste(
  "MM",
  modNames,
  sep = ""
)

names(MMPvalue) <- paste(
  "p.MM",
  modNames,
  sep = ""
)

geneTraitSignificance <- as.data.frame(
  cor(
    datExpr,
    data_clin,
    use = "p"
  )
)

GSPvalue <- as.data.frame(
  corPvalueStudent(
    as.matrix(
      geneTraitSignificance
    ),
    nSamples
  )
)

names(geneTraitSignificance) <- paste(
  "GS.",
  colnames(data_clin),
  sep = ""
)

names(GSPvalue) <- paste(
  "p.GS.",
  colnames(data_clin),
  sep = ""
)


############################################################
# 18. Export gene module membership
############################################################

write.csv(
  geneModuleMembership,
  file.path(
    output_dir,
    "geneModuleMembership.csv"
  )
)


############################################################
# 19. Analyze the brown module
############################################################

module <- "brown"

column <- match(
  module,
  modNames
)

brown_moduleGenes <- names(
  net$colors
)[
  which(
    moduleColors == module
  )
]

MM <- abs(
  geneModuleMembership[
    brown_moduleGenes,
    column
  ]
)

GS <- abs(
  geneTraitSignificance[
    brown_moduleGenes,
    1
  ]
)


############################################################
# 20. Plot module membership vs gene significance
############################################################

png(
  file.path(
    output_dir,
    "brown_membership_gene_significance.png"
  ),
  width = 800,
  height = 600
)

par(
  mfrow = c(1, 1)
)

verboseScatterplot(
  MM,
  GS,
  xlab = paste(
    "Module Membership in",
    module,
    "module"
  ),
  ylab = "Gene significance for Basal",
  main = paste(
    "Module membership vs. gene significance\n"
  ),
  cex.main = 1.2,
  cex.lab = 1.2,
  cex.axis = 1.2,
  col = module
)

dev.off()


############################################################
# 21. Identify hub genes in the brown module
############################################################

# Hub gene criteria:
# |MM| > 0.8 and |GS| > 0.2

brown_hub <- brown_moduleGenes[
  (
    GS > 0.2 &
      MM > 0.8
  )
]

length(brown_hub)


############################################################
# 22. Export brown-module hub genes
############################################################

write.csv(
  brown_hub,
  file.path(
    output_dir,
    "brown_hub_gene.csv"
  )
)


############################################################
# 23. Export gene-module assignments
############################################################

moduleGenes <- names(
  net$colors
)

a = data.frame(
  moduleGenes,
  moduleColors
)

write.csv(
  a,
  file.path(
    output_dir,
    "module_gene.csv"
  )
)

text = unique(
  moduleColors
)


