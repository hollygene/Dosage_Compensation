# make a heatmap 
# read in data 
dosageSensGenes <- read.csv("/Users/hollymcqueary/Dropbox/McQueary/Dosage-Compensation/Paper_Drafts/Dosage_Sensitive_Genes/chrVallDNonE_GC.csv")

library(pheatmap)
library(RColorBrewer)
library(viridis)

chr <- rownames(dosageSensGenes)

col_groups <- substr(chr, 1, 1)
table(col_groups)

# Data frame with column annotations.
mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(dosageSensGenes)

# List with colors for each annotation.
mat_colors <- list(group = brewer.pal(3, "Set1"))
names(mat_colors$group) <- unique(col_groups)

pheatmap(
  mat               = mat,
  color             = inferno(10),
  border_color      = NA,
  show_colnames     = FALSE,
  show_rownames     = FALSE,
  annotation_col    = mat_col,
  annotation_colors = mat_colors,
  drop_levels       = TRUE,
  fontsize          = 14,
  main              = "Default Heatmap"
)