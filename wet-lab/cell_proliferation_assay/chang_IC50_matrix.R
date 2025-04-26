# matt (04.24.2025)

# load packages
library(ggplot2)
library(pheatmap)

# import matrix (.txt) containing IC50 values
import_file_path = "/Users/matthewchang/Library/CloudStorage/GoogleDrive-msc148@case.edu/.shortcut-targets-by-id/1jekwTrBV3zPxqzw9vlL6ZGaqT_qMi-Mk/GryderDrive/Projects/ChemBio/Sub-projects/CBPp300i/Incucyte/Maya_Incucyte/IC50_heatmap/rms_ic50_matrix.txt"
rms_ic50_matrix = read.table(import_file_path,
                             header = T,
                             sep = "\t",
                             row.names = 1,
                             check.names = F)

# drop C2C12 and HSMM
# rms_ic50_matrix = rms_ic50_matrix[, -((ncol(rms_ic50_matrix)-1):ncol(rms_ic50_matrix))]

# add annotations
column_annot = data.frame(Subtype = c(rep("FP-RMS", 5), rep("FN-RMS", 2)))
rownames(column_annot) = colnames(rms_ic50_matrix)
row_annot = data.frame(Class = c("BRDi", "BRDi", "HATi", "HATi", "Degrader", "Degrader", "Degrader"))
rownames(row_annot) = rownames(rms_ic50_matrix)
ann_colors = list(
  Subtype = c("FP-RMS" = "#FF9999", "FN-RMS" = "#9999FF"),
  Class = c("BRDi" = "#FFD700", "HATi" = "#32CD32", "Degrader" = "brown")
)


# custom color scale
colors = colorRampPalette(c("red", "purple2", "black"))(100)
breaks = seq(20, 2000, length.out = 101)

pheatmap(
  rms_ic50_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  angle_col = 315,
  main = "FP-RMS v FN-RMS v myoblast IC50 Heatmap",
  color = colors,
  breaks = breaks,
  annotation_col = column_annot,
  annotation_row = row_annot,
  annotation_colors = ann_colors
)




