##########################
# Updated: 02.24.24 (Matt)
##########################
# This script makes PCA plots, heatmaps, bar plots, box/violin plots, and rank plots
##########################

# import libraries
library(plyr)
library(dplyr)
library(pheatmap)
library(dendsort)
library(tidyverse)
library(ggrepel)
library(ggplot2)
library(reshape2)
library(gridExtra)

# load in Log2(TPM + 1) and Log2FC master matrices
setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC/ExpMatrices/")                         # define the TPM matrix path
EXP.coding.matrix = read.table("IHK_samples.coding.norm.matrix.txt", header = T)                                    # read in TPM matrix
EXP.coding.matrix = EXP.coding.matrix[, 1:(ncol(EXP.coding.matrix) - 4)]                                            # drop IHK-45
EXP.log2FC = read.table("IHK_samples.log2FC.txt", header = T)                                                       # read in log2FC matrix

# load in sample condition info
setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC/SampleList/")
sample.class = read.table("sample-list-IHK-info.txt", header = F, stringsAsFactors = F)
colnames(sample.class) = c("sample.name.list", "drugs", "timepoints", "concentrations", "annotations")

# load in geneset list
setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC/GeneSets/")                            # define the geneset list path
gene.list = read.table("MSC_GRYDER_RH4_CR_TFs_CRISPRTop.genelist.txt", sep = "\t", header = F)                                                   # read in geneset list                                                                                       # save sample name

drug_colors = c(
  "DMSO" = "gray40", "NT" = "gray40", 
  "A485" = "blue", "dCBP" = "darkgreen", "IHK44" = "red", 
  "JQAD" = "darkgoldenrod1", "LS" = "purple", "QL" = "darkorange3"
)

project.folder = "/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC"                            # define the project folder path
sample.set = "IHK_samples"   










########################## PCA ##########################

# subset based on time-point (2h, 6h) and dosage (100nM, 1uM)
EXP.2h.and.6h = EXP.coding.matrix
EXP.2h = EXP.coding.matrix %>% select(gene_id, contains("_2h"))
EXP.6h = EXP.coding.matrix %>% select(gene_id, contains("_6h"))
EXP.100nM = EXP.coding.matrix %>% select(gene_id, matches("(_100nM|_NT_|_DMSO_)"))
EXP.1uM = EXP.coding.matrix %>% select(gene_id, matches("(_1uM|_NT_|_DMSO_)"))

sample.class.2h.and.6h = sample.class %>%
  mutate(condition = as.factor(paste(concentrations, timepoints, sep = "_")))
sample.class.2h = sample.class %>% 
  filter(grepl("_2h$", sample.name.list)) %>%
  mutate(condition = as.factor(paste(concentrations, timepoints, sep = "_")))
sample.class.6h = sample.class %>% 
  filter(grepl("_6h$", sample.name.list)) %>%
  mutate(condition = as.factor(paste(concentrations, timepoints, sep = "_")))
sample.class.100nM = sample.class %>% 
  filter(grepl("(_100nM|_NT_|_DMSO_)", sample.name.list)) %>%
  mutate(condition = as.factor(paste(concentrations, timepoints, sep = "_")))
sample.class.1uM = sample.class %>% 
  filter(grepl("(_1uM|_NT_|_DMSO_)", sample.name.list)) %>%
  mutate(condition = as.factor(paste(concentrations, timepoints, sep = "_")))


# main function for running PCA
run_pca_and_plot = function(EXP.pca_subset, sample.class_subset, time_label) {
  # prepare PCA matrix
  EXP.pca = as.matrix(EXP.pca_subset[, -1])
  rownames(EXP.pca) = EXP.coding.matrix$gene_id
  EXP.pca = log2(EXP.pca + 1)
  
  # perform PCA
  pca = prcomp(t(EXP.pca))
  
  # acquire principle component vectors
  pcv = round((pca$sdev)^2 / sum(pca$sdev^2) * 100, 2)
  
  # merge with metadata
  EXP.pca.df = as.data.frame(pca$x) %>% rownames_to_column("sample.name.list")
  
  EXP.pca.df.meta <- merge(EXP.pca.df, sample.class_subset, by = "sample.name.list") %>%
    mutate(
      sample.name.list = sample.name.list %>%
        gsub("^RH4_", "", .) %>%
        gsub("_100nM_2h|_100nM_6h|_1uM_2h|_1uM_6h|_2h|_6h", "", .)
    )
  
  # plot PCA
  plot.pca = ggplot(EXP.pca.df.meta, aes(PC1, PC2, colour = drugs, shape = condition)) +
    geom_point(size = 3) +
    xlab(paste0("PC1 (", pcv[1], "%)")) +
    ylab(paste0("PC2 (", pcv[2], "%)")) +
    theme_bw() +
    geom_label_repel(aes(label = sample.name.list), size = 5, show.legend = FALSE) +
    scale_colour_manual(values = drug_colors) +
    theme(axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15)) +
    labs(title = paste0("PCA for ", time_label, " samples"),
         subtitle = "Grouped by drug treatment and condition")
  
  return(plot.pca) 
}
  
# run
plot.pca.2h.and.6h = run_pca_and_plot(EXP.2h.and.6h, sample.class.2h.and.6h, "2h and 6h")
print(plot.pca.2h.and.6h)
plot.pca.2h = run_pca_and_plot(EXP.2h, sample.class.2h, "2h")
print(plot.pca.2h)
plot.pca.6h = run_pca_and_plot(EXP.6h, sample.class.6h, "6h")
print(plot.pca.6h)
plot.pca.100nM = run_pca_and_plot(EXP.100nM, sample.class.100nM, "100nM")
print(plot.pca.100nM)
plot.pca.1uM = run_pca_and_plot(EXP.1uM, sample.class.1uM, "1uM")
print(plot.pca.1uM)



























########################## Log2(TPM + 1) heatmap ##########################

sample.order.and.annotations = c(
  #"RH4_DMSO_2h" = "Control/selective degrader",
  "RH4_DMSO_6h" = "Control/selective degrader",
  #"RH4_NT_2h" = "Control/selective degrader",
  "RH4_NT_6h" = "Control/selective degrader",
  #"RH4_JQAD_100nM_2h" = "Control/selective degrader",
  #"RH4_JQAD_1uM_2h" = "Control/selective degrader",
  "RH4_JQAD_100nM_6h" = "Control/selective degrader",
  "RH4_JQAD_1uM_6h" = "Control/selective degrader",
  #"RH4_LS_100nM_2h" = "Control/selective degrader",
  #"RH4_LS_1uM_2h" = "Control/selective degrader",
  "RH4_LS_100nM_6h" = "Control/selective degrader",
  "RH4_LS_1uM_6h" = "Control/selective degrader",
  #"RH4_QL_100nM_2h" = "Control/selective degrader",
  #"RH4_QL_1uM_2h" = "Control/selective degrader",
  "RH4_QL_100nM_6h" = "Control/selective degrader",
  "RH4_QL_1uM_6h" = "Control/selective degrader",
  #"RH4_A485_100nM_2h" = "Dual inhibitor/degrader",
  #"RH4_A485_1uM_2h" = "Dual inhibitor/degrader",
  "RH4_A485_100nM_6h" = "Dual inhibitor/degrader",
  "RH4_A485_1uM_6h" = "Dual inhibitor/degrader",
  #"RH4_dCBP_100nM_2h" = "Dual inhibitor/degrader",
  #"RH4_dCBP_1uM_2h" = "Dual inhibitor/degrader",
  "RH4_dCBP_100nM_6h" = "Dual inhibitor/degrader",
  "RH4_dCBP_1uM_6h" = "Dual inhibitor/degrader",
  #"RH4_IHK44_100nM_2h" = "Dual inhibitor/degrader",
  #"RH4_IHK44_1uM_2h" = "Dual inhibitor/degrader",
  "RH4_IHK44_100nM_6h" = "Dual inhibitor/degrader",
  "RH4_IHK44_1uM_6h" = "Dual inhibitor/degrader"
)

# prepare heatmap
EXP.expressed.matrix = EXP.coding.matrix
annotations = data.frame(
  sample.name.list = names(sample.order.and.annotations),
  annotations = as.vector(sample.order.and.annotations),
  row.names = gsub("RH4_", "", names(sample.order.and.annotations))
)
EXP.expressed.matrix = EXP.expressed.matrix[, c("gene_id", annotations$sample.name.list)]

# filtering
EXP.expressed.matrix$maxTPM = apply(EXP.expressed.matrix[, 2:(ncol(EXP.expressed.matrix) - 1)], 1, FUN = max)   # calculate the maximum TPM value for each gene across all samples
cutoff.expression.min = 10                                                                                      # set the cutoff threshold for minimal expression
EXP.expressed.matrix = subset(EXP.expressed.matrix, EXP.expressed.matrix$maxTPM > cutoff.expression.min)        # subset the matrix to only include genes with a maximum TPM value above the threshold
EXP.expressed.matrix.TFs = subset(EXP.expressed.matrix, EXP.expressed.matrix$gene_id %in% gene.list$V1)         # subset the matrix to only include genes in the gene.list

# clean up
EXP.expressed.matrix.TFs = EXP.expressed.matrix.TFs[match(gene.list$V1, EXP.expressed.matrix.TFs$gene_id), ]    # order genes by gene.list
rownames(EXP.expressed.matrix.TFs) = EXP.expressed.matrix.TFs$gene_id                                           # rownames set to gene_id
EXP.expressed.matrix.TFs = EXP.expressed.matrix.TFs[, 2:(ncol(EXP.expressed.matrix.TFs) - 1)]                   # dropped first (gene_id) and last (maxTPM)
EXP.expressed.matrix.TFs = as.matrix(t(EXP.expressed.matrix.TFs))                                               # transposed matrix
rownames(EXP.expressed.matrix.TFs) = gsub("RH4_", "", rownames(EXP.expressed.matrix.TFs))                       # removed "RH4_"

# plot
color.scale = colorRampPalette(c("white", "pink", "red"))(100)
pheatmap(log2(EXP.expressed.matrix.TFs + 1), 
         cluster_rows = T, 
         cluster_cols = T, 
         scale = 'none', 
         show_rownames = T, 
         show_colnames = T, 
         angle_col = 315, 
         annotation_row = annotations[, "annotations", drop = FALSE], 
         main = paste("Log2(TPM + 1)"),
         color = color.scale)













########################## Log2FC heatmap ##########################

sample.order.and.annotations.FC = c(
  #"RH4_JQAD_100nM_2h_FC" = "Selective degrader",
  #"RH4_JQAD_100nM_6h_FC" = "Selective degrader",
  #"RH4_JQAD_1uM_2h_FC" = "Selective degrader",
  "RH4_JQAD_1uM_6h_FC" = "Selective degrader",
  #"RH4_LS_100nM_2h_FC" = "Selective degrader",
  #"RH4_LS_100nM_6h_FC" = "Selective degrader",
  #"RH4_LS_1uM_2h_FC" = "Selective degrader",
  "RH4_LS_1uM_6h_FC" = "Selective degrader",
  #"RH4_QL_100nM_2h_FC" = "Selective degrader",
  #"RH4_QL_100nM_6h_FC" = "Selective degrader",
  #"RH4_QL_1uM_2h_FC" = "Selective degrader",
  "RH4_QL_1uM_6h_FC" = "Selective degrader",
  #"RH4_dCBP_100nM_2h_FC" = "Dual inhibitor/degrader",
  #"RH4_dCBP_100nM_6h_FC" = "Dual inhibitor/degrader",
  #"RH4_dCBP_1uM_2h_FC" = "Dual inhibitor/degrader",
  "RH4_dCBP_1uM_6h_FC" = "Dual inhibitor/degrader",
  #"RH4_A485_100nM_2h_FC" = "Dual inhibitor/degrader",
  #"RH4_A485_100nM_6h_FC" = "Dual inhibitor/degrader",
  #"RH4_A485_1uM_2h_FC" = "Dual inhibitor/degrader",
  "RH4_A485_1uM_6h_FC" = "Dual inhibitor/degrader",
  #"RH4_IHK44_100nM_2h_FC" = "Dual inhibitor/degrader",
  #"RH4_IHK44_100nM_6h_FC" = "Dual inhibitor/degrader",
  #"RH4_IHK44_1uM_2h_FC" = "Dual inhibitor/degrader",
  "RH4_IHK44_1uM_6h_FC" = "Dual inhibitor/degrader"
)


# Create annotations dataframe
annotations.FC <- data.frame(
  sample.name.list = names(sample.order.and.annotations.FC),
  annotations = as.vector(sample.order.and.annotations.FC),
  row.names = gsub("RH4_", "", names(sample.order.and.annotations.FC))
)

# Prepare heatmap data
EXP.log2FC.heatmap <- EXP.log2FC
EXP.log2FC.heatmap = EXP.log2FC.heatmap[, c("gene_id", annotations.FC$sample.name.list)]

# Filtering and ordering
EXP.log2FC.heatmap <- EXP.log2FC.heatmap[EXP.log2FC.heatmap$gene_id %in% gene.list$V1, ]  # Filter for genes in gene.list
EXP.log2FC.heatmap <- EXP.log2FC.heatmap[match(gene.list$V1, EXP.log2FC.heatmap$gene_id), ] # Order by gene.list
rownames(EXP.log2FC.heatmap) <- EXP.log2FC.heatmap$gene_id  # Set row names to gene_id
EXP.log2FC.heatmap <- EXP.log2FC.heatmap[, 2:ncol(EXP.log2FC.heatmap)]  # Drop first column (gene_id)

# Transpose matrix for consistency with Log2(TPM + 1) structure
EXP.log2FC.heatmap <- as.matrix(t(EXP.log2FC.heatmap))

# Remove "RH4_" prefix from sample names
rownames(EXP.log2FC.heatmap) <- gsub("RH4_", "", rownames(EXP.log2FC.heatmap))

# Plot heatmap
color.scale <- colorRampPalette(c("red", "pink", "white"))(100)
pheatmap(EXP.log2FC.heatmap, 
         cluster_rows = T, 
         cluster_cols = T, 
         scale = "none", 
         show_rownames = T, 
         show_colnames = T, 
         angle_col = 315, 
         annotation_row = annotations.FC[, "annotations", drop = FALSE],
         main = "Log2FC", 
         color = color.scale)









########################## custom Log2(TPM + 1) heatmaps displaying a singular gene across all drug conditions ##########################

# function to extract condition information
extract_info = function(name) {
  # parse through sample name by "_"
  parts = strsplit(name, "_")[[1]]
  drug = parts[2]
  dosage = ""
  timepoint = parts[3]
  
  # check if timepoint is mislabeled
  if (grepl("uM|nM", parts[3])) {
    dosage = parts[3]
    timepoint = parts[4]
  }
  
  return(list(drug = drug, dosage = dosage, timepoint = timepoint))
}

# load in TPM matrix
setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC/ExpMatrices/")
EXP.coding.matrix = read.table("IHK_samples.coding.norm.matrix.txt", header = T)

# choose gene of interest
gene = "MYOD1"
desired.drug.order = c("JQAD", "LS", "QL", 
                       "A485", "dCBP", "IHK44")

# filter
EXP.filtered.matrix = EXP.coding.matrix[EXP.coding.matrix$gene_id == gene, ]
rownames(EXP.filtered.matrix) = EXP.filtered.matrix$gene_id
EXP.filtered.matrix = EXP.filtered.matrix[, -1]

dmso.nt.heatmap = EXP.filtered.matrix[, 1:4]
rownames(dmso.nt.heatmap) = EXP.filtered.matrix$gene_id
EXP.filtered.matrix = EXP.filtered.matrix[, -(1:4)]

dmso.nt.heatmap_log2 = log2(dmso.nt.heatmap + 1)
EXP.filtered.matrix_log2 = log2(EXP.filtered.matrix + 1)

# extract metadeta
drugs = character(length(colnames(EXP.filtered.matrix_log2)))
dosages = character(length(colnames(EXP.filtered.matrix_log2)))
timepoints = character(length(colnames(EXP.filtered.matrix_log2)))
for (j in seq_along(colnames(EXP.filtered.matrix_log2))) {
  info = extract_info(colnames(EXP.filtered.matrix_log2)[j])
  drugs[j] = info$drug
  dosages[j] = info$dosage
  timepoints[j] = info$timepoint
}
unique.drugs = unique(drugs)
unique.timepoints = unique(timepoints)
unique.dosages = unique(dosages)
unique.drugs = intersect(desired.drug.order, unique.drugs)

# create heatmap data
heatmap.data = matrix(NA, nrow = length(unique.dosages) * length(unique.timepoints), ncol = length(unique.drugs))
rownames(heatmap.data) = paste(rep(unique.dosages, each = length(unique.timepoints)), unique.timepoints, sep = "_")
colnames(heatmap.data) = unique.drugs
for (j in seq_along(colnames(EXP.filtered.matrix_log2))) {
  drug = drugs[j]
  dosage = dosages[j]
  timepoint = timepoints[j]
  
  if (paste(dosage, timepoint, sep = "_") %in% rownames(heatmap.data) && drug %in% colnames(heatmap.data)) {
    heatmap.data[paste(dosage, timepoint, sep = "_"), drug] = EXP.filtered.matrix_log2[1, j]
  }
}
heatmap.data = as.data.frame(heatmap.data)
heatmap.data[] = lapply(heatmap.data, as.numeric)
rownames(heatmap.data) = paste(rep(unique.dosages, each = length(unique.timepoints)), unique.timepoints, sep = "_")

colnames(dmso.nt.heatmap_log2) = c("DMSO_2h", "DMSO_6h", "NT_2h", "NT_6h")

combined_data_log2 = cbind(dmso.nt.heatmap_log2, EXP.filtered.matrix_log2)
local.min = min(combined_data_log2, na.rm = T)
local.max = max(combined_data_log2, na.rm = T)



# gradient scales
custom.breaks.EP300 = seq(6,9, length.out =101)
custom.breaks.CREBBP = seq(6,9, length.out =101)
custom.breaks.GAPDH = seq(8,10.5, length.out =101)
custome.breaks.all = seq(local.min, local.max, length.out = 101)

# choose gradient scale from above
gradient.scale = custome.breaks.all

# plot
color.scale = colorRampPalette(c("blue", "lightblue", "white"))(100)
p1 = pheatmap(t(dmso.nt.heatmap_log2), 
              cluster_rows = F, 
              cluster_cols = F, 
              scale = "none", 
              show_rownames = T, 
              show_colnames = T,
              breaks = gradient.scale,
              main = paste(gene),
              legend = F, 
              color = color.scale)
p2 = pheatmap(heatmap.data, 
              cluster_rows = F, 
              cluster_cols = F, 
              scale = "none", 
              show_rownames = T, 
              show_colnames = T,
              angle_col = 315,
              breaks = gradient.scale, 
              main = paste(gene), 
              color = color.scale)
grid.arrange(p1$gtable, p2$gtable, ncol = 2, widths = c(0.75, 2.5))












########################## box/violin plots across multiple genesets ##########################

# load in Log2FC data
setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC/ExpMatrices")
EXP.log2FC = read.table("IHK_samples.log2FC.txt", header = T, row.names = 1)
colnames(EXP.log2FC) = gsub("^RH4_", "", colnames(EXP.log2FC))

# declare samples of interest
samples =        c("A485_100nM_6h_FC", 
                   "IHK44_100nM_6h_FC",
                   "A485_1uM_6h_FC",
                   "IHK44_1uM_6h_FC")
samples.colors = c("A485_100nM_2h_FC" = "blue", 
                   "A485_100nM_6h_FC" = "blue", 
                   "A485_1uM_2h_FC" = "blue", 
                   "A485_1uM_6h_FC" = "blue", 
                   "IHK44_100nM_2h_FC" = "#F6CDCD",
                   "IHK44_100nM_6h_FC" = "#EA4444",
                   "IHK44_1uM_2h_FC" = "#EE7272",
                   "IHK44_1uM_6h_FC" = "#FF0000",
                   "dCBP_100nM_2h_FC" = "green",
                   "dCBP_100nM_6h_FC" = "green",
                   "dCBP_1uM_2h_FC" = "green",
                   "dCBP_1uM_6h_FC" = "green")

# select gene sets
setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC/GeneSets/")
housekeeping_genes = read.table("house_keeping_genes.txt", header = F, stringsAsFactors = F) %>% mutate(gene_set = "Housekeeping genes")
RH4_CR_TFs = read.table("GRYDER_RH4_CR_TFs.genelist.txt", header = F, stringsAsFactors = F) %>% mutate(gene_set = "RH4 CR TFs")
all_genesets = bind_rows(housekeeping_genes, RH4_CR_TFs)

# filter
EXP.log2FC.filter = EXP.log2FC[rownames(EXP.log2FC) %in% all_genesets$V1, ]                                     # subset genes
EXP.log2FC.filter = EXP.log2FC.filter[, colnames(EXP.log2FC.filter) %in% samples]                               # subset samples
EXP.log2FC.filter = EXP.log2FC.filter[, samples]                                                                # set the order of the samples

# clean up
EXP.log2FC.plot = as.data.frame(EXP.log2FC.filter)                                                              # copy over filtered df
EXP.log2FC.plot = rownames_to_column(EXP.log2FC.plot, var = "gene")                                             # transfer rownames (genes) to a column
EXP.log2FC.plot = pivot_longer(EXP.log2FC.plot, cols = -gene, names_to = "sample_name", values_to = "log2FC")   # transform to long format (gene, sample_name, log2FC)
EXP.log2FC.plot = left_join(EXP.log2FC.plot, all_genesets, by = c("gene" = "V1"))
EXP.log2FC.plot$sample_name = factor(EXP.log2FC.plot$sample_name, levels = samples)

# plot
ggplot(EXP.log2FC.plot, aes(x = sample_name, y = log2FC, fill = sample_name)) +
  geom_violin(alpha = 0.6, scale = "width", width = 0.5) + 
  geom_boxplot(width = 0.25, color = "black", alpha = 0.7, outlier.shape = NA) +
  facet_wrap(~ gene_set, scales = "fixed") +
  scale_fill_manual(values = samples.colors) +
  labs(x = "Sample", 
       y = "Log2FC", 
       title = "Boxplot of Log2FC for dual inhibitors/degraders") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))














########################## TPM bar plot for a single gene ##########################

# load in TPM data
setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC/ExpMatrices/")
EXP.coding.matrix = read.table("IHK_samples.coding.norm.matrix.txt", header = T)
cols.to.keep = 1:(ncol(EXP.coding.matrix) - 4)
EXP.coding.matrix = EXP.coding.matrix[, cols.to.keep]

# EDIT HERE
gene = "MYCN"
samples =        c("NT_6h", 
                   "DMSO_6h", 
                   # "QL_100nM_6h", 
                   # "LS_100nM_6h", 
                   # "JQAD_100nM_6h", 
                   "A485_100nM_6h", 
                   "IHK44_100nM_6h", 
                   "dCBP_100nM_6h")
samples.colors = c("NT_6h" = "gray60", 
                   "DMSO_6h" = "gray40", 
                   # "JQAD_100nM_6h" = "#40E0D0", 
                   # "LS_100nM_6h" = "#98FB98", 
                   # "QL_100nM_6h" = "#6B8E23", 
                   "A485_100nM_6h" = "#F4C430", 
                   "IHK44_100nM_6h" = "#FF4500", 
                   "dCBP_100nM_6h" = "#9370DB")

# filter
EXP.single.gene = EXP.coding.matrix[EXP.coding.matrix$gene_id == gene, ]
colnames(EXP.single.gene) = gsub("^RH4_", "", colnames(EXP.single.gene))
EXP.single.gene = EXP.single.gene[, c("gene_id", samples)]
EXP.single.gene = pivot_longer(EXP.single.gene, cols = -gene_id, names_to = "sample_name", values_to = "expression")

# plot
ggplot(EXP.single.gene, aes(x = factor(sample_name, levels = samples), y = expression, fill = sample_name)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) + 
  scale_fill_manual(values = samples.colors) +
  labs(x = "Sample", 
       y = "TPM", 
       title = paste("Expression for", gene)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

################ end of bar plots ################ 












########################## rank plot for a single sample ##########################

# load in geneset
setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC/GeneSets/")
geneset = read.table("GRYDER_RH4_CR_TFs.genelist.txt", header = F, stringsAsFactors = F)
geneset.name = "GRYDER_RH4_CR_TFs.genelist.txt"

# load in Log2FC
setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC/ExpMatrices/")
EXP.log2FC = read.table("IHK_samples.log2FC.txt", header = T, row.names = 1)
# colnames(EXP.log2FC) = gsub("^RH4_", "", colnames(EXP.log2FC))

# define sample
sample = "RH4_IHK44_100nM_6h_FC"

# filter
EXP.log2FC.sample = EXP.log2FC[, sample, drop = F]

# rank and label
EXP.log2FC.df = data.frame(gene = rownames(EXP.log2FC.sample), log2FC = as.numeric(EXP.log2FC.sample[, 1]))
EXP.log2FC.df = EXP.log2FC.df[order(EXP.log2FC.df$log2FC), ]
EXP.log2FC.df$rank = seq_len(nrow(EXP.log2FC.df))
EXP.log2FC.df$target = EXP.log2FC.df$gene %in% geneset$V1
EXP.log2FC.df$label <- ifelse(EXP.log2FC.df$target, EXP.log2FC.df$gene, NA)

# plot
ggplot(EXP.log2FC.df, aes(x = rank, y = log2FC)) +
  geom_point(aes(color = target), alpha = 0.5) +
  geom_text_repel(data = subset(EXP.log2FC.df, target),
                  aes(label = gene),
                  size = 3.5, 
                  max.overlaps = Inf) +
  scale_color_manual(values = c("grey", "red"), name = "Target Gene", labels = c("Non-target", "Target")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))







########################## retrieve lists of genes in groups ########################## 

# define condition to classify
fc_col = "RH4_IHK44_100nM_6h_FC"
tpm_treated_col = "RH4_IHK44_100nM_6h"
tpm_dmso_col = "RH4_DMSO_6h"

# merge TPM and FC columns on gene_id
merged_condition_matrix = merge(EXP.coding.matrix[, c("gene_id", tpm_treated_col, tpm_dmso_col)],
                EXP.log2FC[, c("gene_id", fc_col)],
                by = "gene_id")
colnames(merged_condition_matrix) = c("gene_id", "TPM_treated", "TPM_dmso", "log2FC")

# classification function to assign group status for each gene
classify_gene = function(tpm_treated, tpm_dmso, log2FC) {
  if (tpm_treated < 1 && tpm_dmso < 1) {
    return("Not_Expressed")
  } 
  else if (log2FC >= 1) {
    return("Upregulated")
  } 
  else if (log2FC <= -1) {
    return("Downregulated")
  } 
  else { # -1 <= log2FC <= 1
    return("No_Change")
  }
}

# apply classification
merged_condition_matrix$group = mapply(classify_gene, merged_condition_matrix$TPM_treated, merged_condition_matrix$TPM_dmso, merged_condition_matrix$log2FC)
classify_genes = merged_condition_matrix[, c("gene_id", "group")]

head(classify_genes)
table(classify_genes$group)

write.table(classify_genes,
            file = file.path(project.folder, "genes_group_classification_IHK44_100nM_6hrs.txt"),
            sep = "\t",
            quote = F,
            row.names = F)





