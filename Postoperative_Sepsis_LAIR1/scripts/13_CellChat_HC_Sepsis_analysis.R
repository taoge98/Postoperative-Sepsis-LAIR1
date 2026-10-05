############################################################
# CellChat Analysis of HC and Sepsis Groups
#
# Study:
# Integrated Transcriptomic and Functional Analyses Identify
# LAIR1 as a Key Immunoregulatory Molecule in Postoperative Sepsis
#
# This script performs:
# 1. CellChat analysis in HC and Sepsis groups
# 2. Comparison of intercellular communication networks
# 3. Signaling-role analysis
# 4. Monocyte signaling-change analysis
# 5. MIF pathway-specific analysis
############################################################


############################################################
# 1. Load required packages
############################################################

library(Seurat)
library(CellChat)
library(ggplot2)


############################################################
# 2. Define input and output paths
############################################################

input_file <- "data/seurat_object.Rda"

output_dir <- "results/CellChat"

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


############################################################
# 3. Load annotated Seurat object
############################################################

load(input_file)

head(seurat_object@meta.data)

unique(seurat_object$orig.ident)


############################################################
# 4. Seurat v5 compatibility
############################################################

if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
  seurat_object <- JoinLayers(
    seurat_object
  )
}


############################################################
# 5. Define cell-type annotation
############################################################

# If cell_type is not already stored in metadata,
# use the active Seurat identities.

if (!"cell_type" %in% colnames(seurat_object@meta.data)) {
  seurat_object$cell_type <- Idents(
    seurat_object
  )
}


############################################################
# 6. Define HC and Sepsis groups
############################################################

hc_samples <- c(
  "HC1",
  "HC2"
)

seurat_object$group <- ifelse(
  seurat_object$orig.ident %in% hc_samples,
  "HC",
  "Sepsis"
)

print(
  table(
    seurat_object$orig.ident,
    seurat_object$group
  )
)

print(
  table(
    seurat_object$cell_type
  )
)


############################################################
# 7. Define CellChat analysis function
############################################################

run_cellchat <- function(
    seurat_subset,
    group_name
) {
  
  expression_data <- GetAssayData(
    seurat_subset,
    assay = "RNA",
    layer = "data"
  )
  
  metadata <- seurat_subset@meta.data[
    ,
    c(
      "cell_type",
      "group"
    )
  ]
  
  cellchat <- createCellChat(
    object = expression_data
  )
  
  cellchat <- addMeta(
    cellchat,
    meta = metadata
  )
  
  cellchat <- setIdent(
    cellchat,
    ident.use = "cell_type"
  )
  
  message(
    group_name,
    " cell types: ",
    paste(
      levels(cellchat@idents),
      collapse = ", "
    )
  )
  
  cellchat@DB <- CellChatDB.human
  
  cellchat <- subsetData(
    cellchat
  )
  
  cellchat <- identifyOverExpressedGenes(
    cellchat
  )
  
  cellchat <- identifyOverExpressedInteractions(
    cellchat
  )
  
  cellchat <- projectData(
    cellchat,
    PPI.human
  )
  
  cellchat <- computeCommunProb(
    cellchat,
    raw.use = TRUE
  )
  
  cellchat <- filterCommunication(
    cellchat,
    min.cells = 10
  )
  
  cellchat <- computeCommunProbPathway(
    cellchat
  )
  
  cellchat <- aggregateNet(
    cellchat
  )
  
  cellchat <- netAnalysis_computeCentrality(
    cellchat,
    slot.name = "netP"
  )
  
  return(
    cellchat
  )
}


############################################################
# 8. Run CellChat separately for HC and Sepsis
############################################################

seurat_HC <- subset(
  seurat_object,
  subset = group == "HC"
)

seurat_Sepsis <- subset(
  seurat_object,
  subset = group == "Sepsis"
)


cellchat.HC <- run_cellchat(
  seurat_HC,
  group_name = "HC"
)

cellchat.Sepsis <- run_cellchat(
  seurat_Sepsis,
  group_name = "Sepsis"
)


############################################################
# 9. Export communication tables
############################################################

HC.net <- subsetCommunication(
  cellchat.HC
)

Sepsis.net <- subsetCommunication(
  cellchat.Sepsis
)

write.csv(
  HC.net,
  file.path(
    output_dir,
    "HC_CellChat_interactions.csv"
  ),
  row.names = FALSE
)

write.csv(
  Sepsis.net,
  file.path(
    output_dir,
    "Sepsis_CellChat_interactions.csv"
  ),
  row.names = FALSE
)


############################################################
# 10. Merge CellChat objects
############################################################

object.list <- list(
  HC = cellchat.HC,
  Sepsis = cellchat.Sepsis
)

cellchat.merged <- mergeCellChat(
  object.list,
  add.names = names(object.list)
)


############################################################
# 11. Figure H: Number of interactions
############################################################

groupSize.HC <- as.numeric(
  table(
    cellchat.HC@idents
  )
)

groupSize.Sepsis <- as.numeric(
  table(
    cellchat.Sepsis@idents
  )
)


pdf(
  file.path(
    output_dir,
    "Figure_H_Number_of_interactions_HC.pdf"
  ),
  width = 6,
  height = 6
)

netVisual_circle(
  cellchat.HC@net$count,
  vertex.weight = groupSize.HC,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Number of interactions - HC"
)

dev.off()


pdf(
  file.path(
    output_dir,
    "Figure_H_Number_of_interactions_Sepsis.pdf"
  ),
  width = 6,
  height = 6
)

netVisual_circle(
  cellchat.Sepsis@net$count,
  vertex.weight = groupSize.Sepsis,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Number of interactions - Sepsis"
)

dev.off()


############################################################
# 12. Figure I: Incoming/outgoing interaction strength
############################################################

pdf(
  file.path(
    output_dir,
    "Figure_I_Signaling_role_scatter_HC.pdf"
  ),
  width = 5,
  height = 5
)

netAnalysis_signalingRole_scatter(
  cellchat.HC,
  title = "HC"
)

dev.off()


pdf(
  file.path(
    output_dir,
    "Figure_I_Signaling_role_scatter_Sepsis.pdf"
  ),
  width = 5,
  height = 5
)

netAnalysis_signalingRole_scatter(
  cellchat.Sepsis,
  title = "Sepsis"
)

dev.off()


############################################################
# 13. Figure J: Overall signaling patterns
############################################################

pdf(
  file.path(
    output_dir,
    "Figure_J_Overall_signaling_patterns_HC.pdf"
  ),
  width = 7,
  height = 9
)

netAnalysis_signalingRole_heatmap(
  cellchat.HC,
  pattern = "all",
  title = "overall signaling patterns - HC"
)

dev.off()


pdf(
  file.path(
    output_dir,
    "Figure_J_Overall_signaling_patterns_Sepsis.pdf"
  ),
  width = 7,
  height = 9
)

netAnalysis_signalingRole_heatmap(
  cellchat.Sepsis,
  pattern = "all",
  title = "overall signaling patterns - Sepsis"
)

dev.off()


############################################################
# 14. Figure K: Signaling changes of Monocytes
############################################################

p.monocytes <- netAnalysis_signalingChanges_scatter(
  cellchat.merged,
  idents.use = "Monocytes"
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure_K_Monocytes_signaling_changes_Sepsis_vs_HC.pdf"
  ),
  plot = p.monocytes,
  width = 7,
  height = 6
)


############################################################
# 15. Additional comparative analyses
############################################################

p.number <- compareInteractions(
  cellchat.merged,
  show.legend = FALSE,
  group = c(1, 2),
  measure = "count"
)

ggsave(
  filename = file.path(
    output_dir,
    "CellChat_compare_number_of_interactions.pdf"
  ),
  plot = p.number,
  width = 5,
  height = 5
)


p.weight <- compareInteractions(
  cellchat.merged,
  show.legend = FALSE,
  group = c(1, 2),
  measure = "weight"
)

ggsave(
  filename = file.path(
    output_dir,
    "CellChat_compare_interaction_strength.pdf"
  ),
  plot = p.weight,
  width = 5,
  height = 5
)


pdf(
  file.path(
    output_dir,
    "CellChat_differential_number_of_interactions.pdf"
  ),
  width = 7,
  height = 7
)

netVisual_diffInteraction(
  cellchat.merged,
  weight.scale = TRUE,
  measure = "count"
)

dev.off()


pdf(
  file.path(
    output_dir,
    "CellChat_differential_interaction_strength.pdf"
  ),
  width = 7,
  height = 7
)

netVisual_diffInteraction(
  cellchat.merged,
  weight.scale = TRUE,
  measure = "weight"
)

dev.off()


############################################################
# 16. Optional pathway-specific analysis: MIF
############################################################

pathways.show <- c(
  "MIF"
)

case.cellchat <- cellchat.Sepsis

groupSize <- as.numeric(
  table(
    case.cellchat@idents
  )
)

groupSize.scaled <- groupSize / 100


pdf(
  file.path(
    output_dir,
    "MIF_Sepsis_hierarchy.pdf"
  ),
  width = 7,
  height = 6
)

netVisual_aggregate(
  case.cellchat,
  signaling = pathways.show,
  vertex.receiver = c(
    1,
    2,
    3,
    4
  ),
  layout = "hierarchy",
  vertex.size = groupSize.scaled
)

dev.off()


pdf(
  file.path(
    output_dir,
    "MIF_Sepsis_circle.pdf"
  ),
  width = 7,
  height = 7
)

netVisual_aggregate(
  case.cellchat,
  signaling = pathways.show,
  layout = "circle",
  vertex.size = groupSize.scaled
)

dev.off()


pdf(
  file.path(
    output_dir,
    "MIF_Sepsis_chord.pdf"
  ),
  width = 8,
  height = 8
)

netVisual_aggregate(
  case.cellchat,
  signaling = pathways.show,
  layout = "chord",
  vertex.size = groupSize.scaled
)

dev.off()


pdf(
  file.path(
    output_dir,
    "MIF_Sepsis_heatmap.pdf"
  ),
  width = 7,
  height = 6
)

netVisual_heatmap(
  case.cellchat,
  signaling = pathways.show,
  color.heatmap = "Blues"
)

dev.off()


p.MIF.expression <- plotGeneExpression(
  case.cellchat,
  signaling = "MIF"
)

ggsave(
  filename = file.path(
    output_dir,
    "MIF_Sepsis_gene_expression.pdf"
  ),
  plot = p.MIF.expression,
  width = 8,
  height = 5
)


pdf(
  file.path(
    output_dir,
    "MIF_Sepsis_LR_contribution.pdf"
  ),
  width = 6,
  height = 5
)

netAnalysis_contribution(
  case.cellchat,
  signaling = pathways.show
)

dev.off()


############################################################
# 17. Optional NMF-based communication-pattern analysis
############################################################

# The number of communication patterns should be selected
# after inspecting selectK() results.
#
# Example:
#
# selectK(cellchat.Sepsis, pattern = "outgoing")
#
# nPatterns.outgoing <- 5
#
# cellchat.Sepsis <- identifyCommunicationPatterns(
#   cellchat.Sepsis,
#   pattern = "outgoing",
#   k = nPatterns.outgoing
# )
#
# netAnalysis_river(
#   cellchat.Sepsis,
#   pattern = "outgoing"
# )
#
# netAnalysis_dot(
#   cellchat.Sepsis,
#   pattern = "outgoing"
# )
#
# nPatterns.incoming <- 5
#
# cellchat.Sepsis <- identifyCommunicationPatterns(
#   cellchat.Sepsis,
#   pattern = "incoming",
#   k = nPatterns.incoming
# )
#
# netAnalysis_river(
#   cellchat.Sepsis,
#   pattern = "incoming"
# )
#
# netAnalysis_dot(
#   cellchat.Sepsis,
#   pattern = "incoming"
# )


############################################################
# 18. Save CellChat objects
############################################################

saveRDS(
  cellchat.HC,
  file = file.path(
    output_dir,
    "cellchat_HC.rds"
  )
)

saveRDS(
  cellchat.Sepsis,
  file = file.path(
    output_dir,
    "cellchat_Sepsis.rds"
  )
)

saveRDS(
  cellchat.merged,
  file = file.path(
    output_dir,
    "cellchat_merged_HC_Sepsis.rds"
  )
)


############################################################
# 19. Save session information
############################################################

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "CellChat_sessionInfo.txt"
  )
)
