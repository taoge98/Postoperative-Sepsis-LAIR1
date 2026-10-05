############################################################
# GSE167363 Single-Cell RNA-seq Analysis
#
# Study:
# Integrated Transcriptomic and Functional Analyses Identify
# LAIR1 as a Key Immunoregulatory Molecule in Postoperative Sepsis
############################################################

# =========================
# 1. Load required packages
# =========================
library(Seurat)
library(harmony)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(scales)


# =========================
# 2. Define input and output paths
# =========================
data_root <- "data/GSE167363"
output_dir <- "results/GSE167363"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =========================
# 3. Define sample information
# =========================
sample_ids <- c(
  "HC1", "HC2",
  "NS_ES_0h", "NS_ES_6h", "NS_LS_0h", "NS_LS_6h",
  "S1_0h", "S1_6h", "S2_0h", "S2_6h", "S3_0h", "S3_6h"
)

sample_group <- c(
  "HC", "HC",
  "NS", "NS", "NS", "NS",
  "S", "S", "S", "S", "S", "S"
)

sample_info <- data.frame(
  sample = sample_ids,
  group = sample_group,
  stringsAsFactors = FALSE
)

# =========================
# 4. Import 10X data
# =========================

data_dir <- file.path(data_root, "HC1")  ##指定数据所在目录
list.files(data_dir)  ##列出文件名
HC1 <- Read10X(data.dir = data_dir) ##读取数据 
colnames(HC1) <- paste('HC1',colnames(HC1),sep='_')##为细胞ID加上一个前缀，防止数据合并时出现重复名称
head(colnames(HC1))

data_dir <- file.path(data_root, "HC2")  ##指定数据所在目录
list.files(data_dir)  ##列出文件名
HC2 <- Read10X(data.dir = data_dir) ##读取数据 
colnames(HC2) <- paste('HC2',colnames(HC2),sep='_')##为细胞ID加上一个前缀，防止数据合并时出现重复名称
head(colnames(HC2))

data_dir <- file.path(data_root, "NS ES 0h")  ##指定数据所在目录
list.files(data_dir)  ##列出文件名
NS_ES_0h <- Read10X(data.dir = data_dir) ##读取数据 
colnames(NS_ES_0h) <- paste('NS_ES_0h',colnames(NS_ES_0h),sep='_')##为细胞ID加上一个前缀，防止数据合并时出现重复名称
head(colnames(NS_ES_0h))

data_dir <- file.path(data_root, "NS ES 6h")  ##指定数据所在目录
list.files(data_dir)  ##列出文件名
NS_ES_6h <- Read10X(data.dir = data_dir) ##读取数据 
colnames(NS_ES_6h) <- paste('NS_ES_6h',colnames(NS_ES_6h),sep='_')##为细胞ID加上一个前缀，防止数据合并时出现重复名称
head(colnames(NS_ES_6h))

data_dir <- file.path(data_root, "NS LS 0h")  ##指定数据所在目录
list.files(data_dir)  ##列出文件名
NS_LS_0h <- Read10X(data.dir = data_dir) ##读取数据 
colnames(NS_LS_0h) <- paste('NS_LS_0h',colnames(NS_LS_0h),sep='_')##为细胞ID加上一个前缀，防止数据合并时出现重复名称
head(colnames(NS_LS_0h))

data_dir <- file.path(data_root, "NS LS 6h")  ##指定数据所在目录
list.files(data_dir)  ##列出文件名
NS_LS_6h <- Read10X(data.dir = data_dir) ##读取数据 
colnames(NS_LS_6h) <- paste('NS_LS_6h',colnames(NS_LS_6h),sep='_')##为细胞ID加上一个前缀，防止数据合并时出现重复名称
head(colnames(NS_LS_6h))

data_dir <- file.path(data_root, "S1 0h")  ##指定数据所在目录
list.files(data_dir)  ##列出文件名
S1_0h <- Read10X(data.dir = data_dir) ##读取数据 
colnames(S1_0h) <- paste('S1_0h',colnames(S1_0h),sep='_')##为细胞ID加上一个前缀，防止数据合并时出现重复名称
head(colnames(S1_0h))

data_dir <- file.path(data_root, "S1 6h")  ##指定数据所在目录
list.files(data_dir)  ##列出文件名
S1_6h <- Read10X(data.dir = data_dir) ##读取数据 
colnames(S1_6h) <- paste('S1_6h',colnames(S1_6h),sep='_')##为细胞ID加上一个前缀，防止数据合并时出现重复名称
head(colnames(S1_6h))

data_dir <- file.path(data_root, "S2 0h")  ##指定数据所在目录
list.files(data_dir)  ##列出文件名
S2_0h <- Read10X(data.dir = data_dir) ##读取数据 
colnames(S2_0h) <- paste('S2_0h',colnames(S2_0h),sep='_')##为细胞ID加上一个前缀，防止数据合并时出现重复名称
head(colnames(S2_0h))

data_dir <- file.path(data_root, "S2 6h")  ##指定数据所在目录
list.files(data_dir)  ##列出文件名
S2_6h <- Read10X(data.dir = data_dir) ##读取数据 
colnames(S2_6h) <- paste('S2_6h',colnames(S2_6h),sep='_')##为细胞ID加上一个前缀，防止数据合并时出现重复名称
head(colnames(S2_6h))

data_dir <- file.path(data_root, "S3 0h")  ##指定数据所在目录
list.files(data_dir)  ##列出文件名
S3_0h <- Read10X(data.dir = data_dir) ##读取数据 
colnames(S3_0h) <- paste('S3_0h',colnames(S3_0h),sep='_')##为细胞ID加上一个前缀，防止数据合并时出现重复名称
head(colnames(S3_0h))

data_dir <- file.path(data_root, "S3 6h")  ##指定数据所在目录
list.files(data_dir)  ##列出文件名
S3_6h <- Read10X(data.dir = data_dir) ##读取数据 
colnames(S3_6h) <- paste('S3_6h',colnames(S3_6h),sep='_')##为细胞ID加上一个前缀，防止数据合并时出现重复名称
head(colnames(S3_6h))

HC1 <- CreateSeuratObject(counts = HC1, project = "HC1", min.cells = 3, min.features = 200)
HC2 <- CreateSeuratObject(counts = HC2, project = "HC2", min.cells = 3, min.features = 200)
NS_ES_0h <- CreateSeuratObject(counts = NS_ES_0h, project = "NS_ES_0h", min.cells = 3, min.features = 200)
NS_ES_6h <- CreateSeuratObject(counts = NS_ES_6h, project = "NS_ES_6h", min.cells = 3, min.features = 200)
NS_LS_0h <- CreateSeuratObject(counts = NS_LS_0h, project = "NS_LS_0h", min.cells = 3, min.features = 200)
NS_LS_6h <- CreateSeuratObject(counts = NS_LS_6h, project = "NS_LS_6h", min.cells = 3, min.features = 200)
S1_0h <- CreateSeuratObject(counts = S1_0h, project = "S1_0h", min.cells = 3, min.features = 200)
S1_6h <- CreateSeuratObject(counts = S1_6h, project = "S1_6h", min.cells = 3, min.features = 200)
S2_0h <- CreateSeuratObject(counts = S2_0h, project = "S2_0h", min.cells = 3, min.features = 200)
S2_6h <- CreateSeuratObject(counts = S2_6h, project = "S2_6h", min.cells = 3, min.features = 200)
S3_0h <- CreateSeuratObject(counts = S3_0h , project = "S3_0h", min.cells = 3, min.features = 200)
S3_6h <- CreateSeuratObject(counts = S3_6h, project = "S3_6h", min.cells = 3, min.features = 200)
seurat_object <- merge(HC1,y=c(HC2,NS_ES_0h,NS_ES_6h,NS_LS_0h,
                               NS_LS_6h,S1_0h,S1_6h,S2_0h,S2_6h,S3_0h,S3_6h), add.cell.ids = c("HC1","HC2","NS_ES_0h","NS_ES_6h","NS_LS_0h",
                                                                                               "NS_LS_6h","S1_0h","S1_6h","S2_0h","S2_6h","S3_0h","S3_6h"), project = "Twelve_merged")

# =========================
# 5. Quality control
# =========================
seurat_object[["percent.mt"]] <- PercentageFeatureSet(
  seurat_object,
  pattern = "^MT-"
)

# Optional QC plot before filtering
p_qc_before <- VlnPlot(
  seurat_object,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)

ggsave(
  filename = file.path(output_dir, "QC_before_filtering.pdf"),
  plot = p_qc_before,
  width = 10,
  height = 4
)

# Retain cells using the original QC thresholds:
#   500 < detected genes < 5000
#   total RNA counts < 25000
#   mitochondrial transcript percentage < 15%
seurat_object <- subset(
  seurat_object,
  subset = nFeature_RNA > 500 &
    nFeature_RNA < 5000 & nCount_RNA<25000 &
    percent.mt < 15
)

# QC scatter plots after filtering
p_qc1 <- FeatureScatter(
  seurat_object,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
)

p_qc2 <- FeatureScatter(
  seurat_object,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)

ggsave(
  filename = file.path(output_dir, "QC_scatter_after_filtering.pdf"),
  plot = p_qc1 + p_qc2,
  width = 9,
  height = 4
)

# =========================
# 6. Normalization, variable-feature selection, scaling, and PCA
# =========================
seurat_object <- NormalizeData(
  seurat_object,
  normalization.method = "LogNormalize",
  scale.factor = 10000
) %>%
  FindVariableFeatures(
    selection.method = "vst",
    nfeatures = 2000
  ) %>%
  ScaleData() %>%
  RunPCA()

# =========================
# 7. Harmony batch correction
# =========================
seurat_object <- RunHarmony(
  object = seurat_object,
  group.by.vars = "orig.ident",
  plot_convergence = TRUE
)

# =========================
# 8. UMAP and clustering
# =========================
seurat_object <- seurat_object %>%
  RunUMAP(
    reduction = "harmony",
    dims = 1:15
  ) %>%
  FindNeighbors(
    reduction = "harmony",
    dims = 1:15
  ) %>%
  FindClusters(
    resolution = 0.5
  )

#降维可视化
p1 <- DimPlot(seurat_object, reduction = "umap", group.by = "orig.ident", pt.size = 0.5)
p2 <- DimPlot(seurat_object, reduction = "umap", pt.size = 0.5)
p1+p2

# =========================
# 9. Marker-gene analysis
# =========================
# Join assay layers for Seurat v5 compatibility.
if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
  seurat_object <- JoinLayers(seurat_object)
}

all_markers <- FindAllMarkers(
  seurat_object,
  only.pos = TRUE,
  logfc.threshold = 0.6,
  min.pct = 0.4
)

top5_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 5,
    with_ties = FALSE
  ) %>%
  ungroup()

write.csv(
  all_markers,
  file = file.path(output_dir, "All_cluster_markers.csv"),
  row.names = FALSE
)

write.csv(
  top5_markers,
  file = file.path(output_dir, "Top5_markers_per_cluster.csv"),
  row.names = FALSE
)

# =========================
# 10. Cell-type annotation
# =========================
cluster_map <- c(
  "0"  = "T cells",
  "1"  = "B cells",
  "2"  = "Monocytes",
  "3"  = "NK cells",
  "4"  = "T cells",
  "5"  = "Monocytes",
  "6"  = "B cells",
  "7"  = "Platelets",
  "8"  = "Neutrophils",
  "9"  = "Erythroid",
  "10" = "Dendritic Cells",
  "11" = "B cells",
  "12" = "Proliferating cells",
  "13" = "Monocytes",
  "14" = "Plasma cells",
  "15" = "Monocytes"
)

# Check that all detected clusters have an annotation
detected_clusters <- levels(Idents(seurat_object))
missing_annotation <- setdiff(detected_clusters, names(cluster_map))

if (length(missing_annotation) > 0) {
  stop(
    "The following clusters are missing from cluster_map: ",
    paste(missing_annotation, collapse = ", ")
  )
}

seurat_object <- RenameIdents(
  seurat_object,
  cluster_map
)

seurat_object$celltype <- Idents(seurat_object)
Idents(seurat_object) <- "celltype"

# =========================
# 11. Add clinical groups
# =========================
group_map <- setNames(sample_info$group, sample_info$sample)

seurat_object$group <- unname(
  group_map[as.character(seurat_object$orig.ident)]
)

if (any(is.na(seurat_object$group))) {
  stop("Some cells could not be assigned to HC, NS, or S.")
}

# Keep the same group order as used in the manuscript figure
seurat_object$group <- factor(
  seurat_object$group,
  levels = c("HC", "S", "NS")
)

# =========================
# 12. annotated UMAP
# =========================
p_umap <- DimPlot(
  seurat_object,
  reduction = "umap",
  group.by = "celltype",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.15
) +
  labs(
    x = "UMAP_1",
    y = "UMAP_2"
  ) +
  theme_classic() +
  theme(
    legend.title = element_blank()
  )

ggsave(
  filename = file.path(output_dir, "Figure_A_celltype_UMAP.pdf"),
  plot = p_umap,
  width = 8.5,
  height = 5.5
)

# =========================
# 13. cell-type composition by sample
# =========================
cell_stat <- as.data.frame(
  table(
    sample = seurat_object$orig.ident,
    celltype = seurat_object$celltype
  )
)

colnames(cell_stat) <- c("sample", "celltype", "Freq")

cell_stat <- cell_stat %>%
  group_by(sample) %>%
  mutate(
    Proportion = Freq / sum(Freq),
    Label = percent(Proportion, accuracy = 1)
  ) %>%
  ungroup()

# Preserve sample order
cell_stat$sample <- factor(
  cell_stat$sample,
  levels = sample_ids
)

cell_colors <- c(
  "T cells" = "#00468B",
  "B cells" = "#925E9F",
  "Monocytes" = "#759EDD",
  "NK cells" = "#0099B4",
  "Platelets" = "#76D1B1",
  "Neutrophils" = "#42B540",
  "Erythroid" = "#B8D24D",
  "Dendritic Cells" = "#EDE447",
  "Proliferating cells" = "#FAB158",
  "Plasma cells" = "#FF7777"
)

p_composition <- ggplot(
  cell_stat,
  aes(
    x = sample,
    y = Proportion,
    fill = celltype
  )
) +
  geom_col(width = 0.72) +
  geom_text(
    aes(
      label = ifelse(
        Proportion > 0.05,
        Label,
        ""
      )
    ),
    position = position_stack(vjust = 0.5),
    size = 2.4,
    color = "white"
  ) +
  scale_fill_manual(
    values = cell_colors,
    drop = FALSE
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1)
  ) +
  labs(
    x = "Sample",
    y = "Cell proportion",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(
      color = "black"
    ),
    legend.title = element_blank()
  )

ggsave(
  filename = file.path(output_dir, "Figure_B_celltype_composition.pdf"),
  plot = p_composition,
  width = 9.5,
  height = 5.2
)

# =========================
# 14. LAIR1 FeaturePlot
# =========================
p_lair1_feature <- FeaturePlot(
  seurat_object,
  features = "LAIR1",
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.15
) +
  labs(
    title = "LAIR1",
    x = "UMAP_1",
    y = "UMAP_2"
  ) +
  theme_classic()

ggsave(
  filename = file.path(output_dir, "Figure_C_LAIR1_FeaturePlot.pdf"),
  plot = p_lair1_feature,
  width = 6.2,
  height = 5.2
)

# =========================
# 15. LAIR1 violin plot across all annotated cell types
# =========================
p_lair1_allcells <- VlnPlot(
  seurat_object,
  features = "LAIR1",
  group.by = "celltype",
  pt.size = 0.1
) +
  labs(
    title = "LAIR1",
    x = NULL,
    y = "Expression Level"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    ),
    legend.position = "none"
  )

ggsave(
  filename = file.path(output_dir, "LAIR1_violin_all_celltypes.pdf"),
  plot = p_lair1_allcells,
  width = 8.5,
  height = 5.2
)

# =========================
# 16. LAIR1 expression in monocytes by group
# =========================
mono_subset <- subset(
  seurat_object,
  idents = "Monocytes"
)

mono_subset$group <- factor(
  mono_subset$group,
  levels = c("HC", "S", "NS")
)

# For Seurat v5 differential-expression testing
if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
  mono_subset <- JoinLayers(mono_subset)
}

Idents(mono_subset) <- "group"

# Differential-expression statistics for LAIR1
stats_S_vs_HC <- FindMarkers(
  mono_subset,
  ident.1 = "S",
  ident.2 = "HC",
  features = "LAIR1",
  logfc.threshold = 0,
  min.pct = 0,
  test.use = "wilcox"
)

stats_NS_vs_HC <- FindMarkers(
  mono_subset,
  ident.1 = "NS",
  ident.2 = "HC",
  features = "LAIR1",
  logfc.threshold = 0,
  min.pct = 0,
  test.use = "wilcox"
)

stats_S_vs_NS <- FindMarkers(
  mono_subset,
  ident.1 = "S",
  ident.2 = "NS",
  features = "LAIR1",
  logfc.threshold = 0,
  min.pct = 0,
  test.use = "wilcox"
)

write.csv(
  stats_S_vs_HC,
  file = file.path(output_dir, "LAIR1_monocytes_S_vs_HC.csv")
)

write.csv(
  stats_NS_vs_HC,
  file = file.path(output_dir, "LAIR1_monocytes_NS_vs_HC.csv")
)

write.csv(
  stats_S_vs_NS,
  file = file.path(output_dir, "LAIR1_monocytes_S_vs_NS.csv")
)

# Extract data for publication-style violin plot
plot_data <- FetchData(
  mono_subset,
  vars = c("LAIR1", "group")
)

plot_data$group <- factor(
  plot_data$group,
  levels = c("HC", "S", "NS")
)

my_comparisons <- list(
  c("HC", "S"),
  c("S", "NS"),
  c("HC", "NS")
)

p_lair1_monocytes <- ggplot(
  plot_data,
  aes(
    x = group,
    y = LAIR1,
    fill = group
  )
) +
  geom_violin(
    trim = FALSE,
    scale = "width",
    adjust = 1.2,
    color = "black",
    linewidth = 0.4
  ) +
  geom_boxplot(
    width = 0.10,
    fill = "white",
    color = "black",
    outlier.shape = NA,
    linewidth = 0.4
  ) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.signif",
    symnum.args = list(
      cutpoints = c(
        0,
        0.0001,
        0.001,
        0.01,
        0.05,
        1
      ),
      symbols = c(
        "****",
        "***",
        "**",
        "*",
        "NS"
      )
    )
  ) +
  scale_fill_manual(
    values = c(
      "HC" = "#3182BD",
      "S" = "#E6550D",
      "NS" = "#9E9E9E"
    )
  ) +
  labs(
    x = NULL,
    y = "Log-normalized Expression"
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(
      color = "black"
    ),
    legend.position = "none"
  )

ggsave(
  filename = file.path(output_dir, "LAIR1_monocytes_by_group.pdf"),
  plot = p_lair1_monocytes,
  width = 4.6,
  height = 5.2
)

# =========================
# 17. Save processed Seurat object
# =========================
saveRDS(
  seurat_object,
  file = file.path(output_dir, "processed_seurat_object.rds")
)

# =========================
# 18. Save session information
# =========================
capture.output(
  sessionInfo(),
  file = file.path(output_dir, "sessionInfo.txt")
)

