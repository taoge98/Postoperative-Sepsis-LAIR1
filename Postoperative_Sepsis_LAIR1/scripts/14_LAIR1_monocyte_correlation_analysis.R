############################################################
# LAIR1 Correlation Analysis in Monocytes
#
# Study:
# Integrated Transcriptomic and Functional Analyses Identify
# LAIR1 as a Key Immunoregulatory Molecule in Postoperative Sepsis
#
# This script evaluates the association between LAIR1 and
# immunoregulatory genes in monocytes using:
# 1. Cell-level Spearman correlation
# 2. 50-bin mean-expression Pearson correlation
# 3. Sensitivity analyses using 30, 50, and 80 bins
############################################################


############################################################
# 1. Load required packages
############################################################

library(Seurat)
library(dplyr)
library(purrr)
library(ggplot2)
library(patchwork)


############################################################
# 2. Define input and output paths
############################################################

input_file <- "data/seurat_object.Rda"

output_dir <- "results"

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


############################################################
# 3. Load processed Seurat object
############################################################

load(input_file)

set.seed(1234)


############################################################
# 4. Extract monocytes
############################################################

# If cell-type annotations are stored in metadata,
# use them to define identities.

if ("cell_type" %in% colnames(seurat_object@meta.data)) {
  Idents(seurat_object) <- "cell_type"
}

monocyte_object <- subset(
  seurat_object,
  idents = "Monocytes"
)

# Seurat v5 compatibility
if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
  monocyte_object <- JoinLayers(monocyte_object)
}


############################################################
# 5. Extract normalized expression data
############################################################

genes_to_check <- c(
  "LAIR1",
  "IL10",
  "MRC1",
  "CD163",
  "CD274"
)

# Keep only genes present in the dataset
genes_present <- genes_to_check[
  genes_to_check %in% rownames(monocyte_object)
]

genes_missing <- setdiff(
  genes_to_check,
  genes_present
)

if (length(genes_missing) > 0) {
  message(
    "Genes not detected in the Seurat object: ",
    paste(
      genes_missing,
      collapse = ", "
    )
  )
}

if (!"LAIR1" %in% genes_present) {
  stop(
    "LAIR1 is not present in the Seurat object."
  )
}

# Extract log-normalized expression values
expr_data <- FetchData(
  monocyte_object,
  vars = genes_present,
  layer = "data"
)


############################################################
# 6. Cell-level Spearman correlation analysis
############################################################

# This section reproduces the exploratory cell-level
# correlation analysis.

get_cell_level_cor <- function(
    target_gene,
    data
) {
  
  other_genes <- setdiff(
    colnames(data),
    target_gene
  )
  
  results <- map_dfr(
    other_genes,
    function(g) {
      
      res <- cor.test(
        data[[target_gene]],
        data[[g]],
        method = "spearman",
        exact = FALSE
      )
      
      data.frame(
        Gene = g,
        Correlation = unname(
          res$estimate
        ),
        P_value = res$p.value
      )
    }
  )
  
  return(results)
}


cell_level_cor_results <- get_cell_level_cor(
  target_gene = "LAIR1",
  data = expr_data
)

write.csv(
  cell_level_cor_results,
  file = file.path(
    output_dir,
    "LAIR1_cell_level_Spearman_correlations.csv"
  ),
  row.names = FALSE
)


############################################################
# 7. Function for equal-sized binning
############################################################

get_binned_df <- function(
    data,
    driver,
    target,
    n_bins
) {
  
  data %>%
    select(
      all_of(
        c(
          driver,
          target
        )
      )
    ) %>%
    arrange(
      .data[[driver]]
    ) %>%
    mutate(
      bin = ntile(
        .data[[driver]],
        n_bins
      )
    ) %>%
    group_by(bin) %>%
    summarise(
      LAIR1_mean = mean(
        .data[[driver]],
        na.rm = TRUE
      ),
      Target_mean = mean(
        .data[[target]],
        na.rm = TRUE
      ),
      n_cells = n(),
      .groups = "drop"
    ) %>%
    mutate(
      Target_Gene = target,
      Number_of_bins = n_bins
    )
}


############################################################
# 8. Primary analysis: 50 equal-sized bins
############################################################

target_genes <- intersect(
  c(
    "IL10",
    "MRC1",
    "CD163",
    "CD274"
  ),
  colnames(expr_data)
)

bin_num <- 50

binned_results_50 <- map_dfr(
  target_genes,
  function(g) {
    
    get_binned_df(
      data = expr_data,
      driver = "LAIR1",
      target = g,
      n_bins = bin_num
    )
  }
)

write.csv(
  binned_results_50,
  file = file.path(
    output_dir,
    "LAIR1_binned_expression_50bins.csv"
  ),
  row.names = FALSE
)


############################################################
# 9. Correlation analysis based on 50-bin mean expression
############################################################

# Pearson correlation is applied to bin-level mean
# expression values.

get_binned_cor <- function(
    binned_data,
    method = "pearson"
) {
  
  binned_data %>%
    group_by(Target_Gene) %>%
    group_modify(
      ~ {
        
        test_result <- cor.test(
          .x$LAIR1_mean,
          .x$Target_mean,
          method = method
        )
        
        tibble(
          Correlation = unname(
            test_result$estimate
          ),
          P_value = test_result$p.value,
          Method = method,
          Number_of_bins = unique(
            .x$Number_of_bins
          )
        )
      }
    ) %>%
    ungroup()
}


binned_cor_results_50 <- get_binned_cor(
  binned_results_50,
  method = "pearson"
)

write.csv(
  binned_cor_results_50,
  file = file.path(
    output_dir,
    "LAIR1_binned_Pearson_correlations_50bins.csv"
  ),
  row.names = FALSE
)

print(
  binned_cor_results_50
)


############################################################
# 10. Plot individual binned correlations
############################################################

plot_binned_correlation <- function(
    target_gene,
    binned_data,
    correlation_results
) {
  
  plot_data <- binned_data %>%
    filter(
      Target_Gene == target_gene
    )
  
  stat_data <- correlation_results %>%
    filter(
      Target_Gene == target_gene
    )
  
  r_value <- stat_data$Correlation
  
  p_value <- stat_data$P_value
  
  annotation_text <- paste0(
    "R=",
    formatC(
      r_value,
      format = "f",
      digits = 2
    ),
    ", P=",
    format.pval(
      p_value,
      digits = 2,
      eps = 1e-300
    )
  )
  
  ymax <- max(
    plot_data$Target_mean,
    na.rm = TRUE
  )
  
  ymin <- min(
    plot_data$Target_mean,
    na.rm = TRUE
  )
  
  xpos <- min(
    plot_data$LAIR1_mean,
    na.rm = TRUE
  ) +
    0.05 * diff(
      range(
        plot_data$LAIR1_mean,
        na.rm = TRUE
      )
    )
  
  ypos <- ymax -
    0.05 * (
      ymax - ymin
    )
  
  p <- ggplot(
    plot_data,
    aes(
      x = LAIR1_mean,
      y = Target_mean
    )
  ) +
    geom_point(
      color = "steelblue",
      size = 2.3,
      alpha = 0.85
    ) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      color = "red",
      se = TRUE
    ) +
    annotate(
      "text",
      x = xpos,
      y = ypos,
      label = annotation_text,
      color = "red",
      hjust = 0,
      vjust = 1,
      size = 5
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(
        hjust = 0.5
      )
    ) +
    labs(
      title = paste(
        "LAIR1 vs",
        target_gene
      ),
      x = "Mean LAIR1 Expression (per bin)",
      y = paste0(
        "Mean ",
        target_gene,
        " Expression (per bin)"
      )
    )
  
  return(p)
}


plot_list_50 <- setNames(
  lapply(
    target_genes,
    function(g) {
      
      plot_binned_correlation(
        target_gene = g,
        binned_data = binned_results_50,
        correlation_results = binned_cor_results_50
      )
    }
  ),
  target_genes
)


############################################################
# 11. Save individual correlation plots
############################################################

walk2(
  plot_list_50,
  names(plot_list_50),
  function(p, g) {
    
    ggsave(
      filename = file.path(
        output_dir,
        paste0(
          "LAIR1_vs_",
          g,
          "_50bins.pdf"
        )
      ),
      plot = p,
      width = 5,
      height = 4.5
    )
  }
)


############################################################
# 12. Combined manuscript panel
############################################################

figure_genes <- intersect(
  c(
    "IL10",
    "MRC1",
    "CD163",
    "CD274"
  ),
  names(plot_list_50)
)

combined_plot <- wrap_plots(
  plot_list_50[
    figure_genes
  ],
  ncol = length(
    figure_genes
  )
)

ggsave(
  filename = file.path(
    output_dir,
    "LAIR1_binned_correlations_main_figure.pdf"
  ),
  plot = combined_plot,
  width = 18,
  height = 4.5
)


############################################################
# 13. Sensitivity analyses using 30, 50, and 80 bins
############################################################

bin_numbers <- c(
  30,
  50,
  80
)

binned_results_all <- map_dfr(
  bin_numbers,
  function(n_bins) {
    
    map_dfr(
      target_genes,
      function(g) {
        
        get_binned_df(
          data = expr_data,
          driver = "LAIR1",
          target = g,
          n_bins = n_bins
        )
      }
    )
  }
)

write.csv(
  binned_results_all,
  file = file.path(
    output_dir,
    "LAIR1_binned_expression_30_50_80bins.csv"
  ),
  row.names = FALSE
)


############################################################
# 14. Correlations for sensitivity analyses
############################################################

sensitivity_cor_results <- binned_results_all %>%
  group_by(
    Number_of_bins,
    Target_Gene
  ) %>%
  group_modify(
    ~ {
      
      test_result <- cor.test(
        .x$LAIR1_mean,
        .x$Target_mean,
        method = "pearson"
      )
      
      tibble(
        Correlation = unname(
          test_result$estimate
        ),
        P_value = test_result$p.value,
        Method = "pearson"
      )
    }
  ) %>%
  ungroup()

write.csv(
  sensitivity_cor_results,
  file = file.path(
    output_dir,
    "LAIR1_binned_correlation_sensitivity_30_50_80bins.csv"
  ),
  row.names = FALSE
)

print(
  sensitivity_cor_results
)


############################################################
# 15. Save R session information
############################################################

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "LAIR1_binned_correlation_sessionInfo.txt"
  )
)
