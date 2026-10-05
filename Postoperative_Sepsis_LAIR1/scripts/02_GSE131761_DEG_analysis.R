############################################################	
# Differential Expression Analysis of GSE131761	
#	
# Study:	
# Integrated Transcriptomic and Functional Analyses Identify	
# LAIR1 as a Key Immunoregulatory Molecule in Postoperative Sepsis	
#	
# This script performs differential expression analysis	
# between the sepsis and non-sepsis groups using limma.	
############################################################	


############################################################	
# 1. Load required packages	
############################################################	

library(limma)	
library(ggplot2)	


############################################################	
# 2. Load expression and group data	
############################################################	

data = read.table(	
  "results GSE131761_gene.txt",
  header = T	
)	

group = read.table(	
  "metadata	GSE131761_group.txt",
  header = TRUE	
)	


############################################################	
# 3. Construct design matrix	
############################################################	

design = model.matrix(~0 + factor(group$group))	

colnames(design) = levels(	
  factor(group$group)	
)	

rownames(design) = colnames(data)	


############################################################	
# 4. Define comparison	
############################################################	

contrast.matrix <- makeContrasts(	
  "sepsis-non_sepsis",	
  levels = design	
)	


############################################################	
# 5. Fit linear model	
############################################################	

## step1	
fit <- lmFit(	
  data,	
  design	
)	

## step2	
fit2 <- contrasts.fit(	
  fit,	
  contrast.matrix	
)	

fit2 <- eBayes(	
  fit2	
)	


############################################################	
# 6. Extract differential expression results	
############################################################	

## step3	
DEG = topTable(	
  fit2,	
  coef = 1,	
  n = Inf	
)	

DEG = na.omit(	
  DEG	
)	


############################################################	
# 7. Classify differentially expressed genes	
############################################################	

# Differential expression criteria:	
# P.Value < 0.05 and |logFC| > 0.5	

DEG$regulate = ifelse(	
  DEG$P.Value < 0.05 & DEG$logFC > 0.5,	
  "up",	
  ifelse(	
    DEG$P.Value < 0.05 & DEG$logFC < -0.5,	
    "down",	
    "unchanged"	
  )	
)	

table(	
  DEG$regulate	
)	


############################################################	
# 8. Export differential expression results	
############################################################	

write.table(	
  DEG,	
  file = "results	GSE74224_DEG_results.txt",
  sep = "\t"	
)	


############################################################	
# 9. Volcano plot	
############################################################	

ggplot(	
  DEG,	
  aes(	
    x = logFC,	
    y = -log10(P.Value)	
  )	
) +	
  geom_point(	
    alpha = 0.6,	
    size = 3.5,	
    aes(	
      color = regulate	
    )	
  ) +	
  ylab("-log10(p.Value)") +	
  scale_color_manual(	
    values = c(	
      "blue",	
      "grey",	
      "red"	
    )	
  ) +	
  geom_vline(	
    xintercept = c(-0.5, 0.5),	
    lty = 4,	
    col = "black",	
    lwd = 0.8	
  ) +	
  theme_bw()	
