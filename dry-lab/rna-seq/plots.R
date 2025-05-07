###################################
# Matt, 05.07.2025
# build PCA, heatmaps, bar plots, box plots
###################################

library(plyr)
library(dplyr)
library(pheatmap)
library(dendsort)
library(tidyverse)
library(ggrepel)
library(ggplot2)
library(reshape2)
library(gridExtra)

# load TPM and Log2FC matrices
setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC/ExpMatrices/")
EXP.coding.matrix = read.table("IHK_samples.coding.norm.matrix.txt", header=T)
EXP.coding.matrix = EXP.coding.matrix[ , 1:29]
EXP.log2FC = read.table("IHK_samples.log2FC.txt", header=T)

project.folder = "/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC"
sample.set = "IHK_samples"   







#########################################################################################################
################################################## PCA ##################################################
#########################################################################################################

# REQUIRED: load metadata
setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC/SampleList/")
sample.class = read.table("sample-list-IHK-info.txt", header = F, stringsAsFactors = F)
colnames(sample.class) = c("sample.name.list", "drugs", "timepoints", "concentrations", "annotations")

# OPTIONAL: choose color scheme
drug_colors = c(
  "DMSO" = "gray40", "NT" = "gray40", 
  "A485" = "blue", "dCBP" = "darkgreen", "IHK44" = "red", 
  "JQAD" = "darkgoldenrod1", "LS" = "purple", "QL" = "darkorange3"
)

# function to perform PCA
run_pca <- function(EXP.pca_subset, sample.class_subset) {
  # prepare matrix
  EXP.pca <- as.matrix(EXP.pca_subset[, -1])
  rownames(EXP.pca) <- EXP.coding.matrix$gene_id
  EXP.pca <- log2(EXP.pca + 1)
  
  # run PCA
  pca <- prcomp(t(EXP.pca))
  
  # variance explained
  pcv <- round((pca$sdev)^2 / sum(pca$sdev^2) * 100, 2)
  
  # PCA coordinates df
  EXP.pca.df <- as.data.frame(pca$x) %>% rownames_to_column("sample.name.list")
  
  # merge metadata
  EXP.pca.df.meta <- merge(EXP.pca.df, sample.class_subset, by = "sample.name.list") %>% 
    mutate(sample.name.list = sample.name.list %>%
             gsub("^RH4_", "", .) %>%
             gsub("_100nM_2h|_100nM_6h|_1uM_2h|_1uM_6h|_2h|_6h", "", .))
  
  return(list(pca = pca, meta = EXP.pca.df.meta, pcv = pcv))
}

# function to plot PCA
plot_pca <- function(pca_result, time_label) {
  pca <- pca_result$pca
  EXP.pca.df.meta <- pca_result$meta
  pcv <- pca_result$pcv
  
  # plot
  plot.pca <- ggplot(EXP.pca.df.meta, aes(PC1, PC2, colour = drugs, shape = condition)) +
    geom_point(size = 3) +
    xlab(paste0("PC1 (", pcv[1], "%)")) +
    ylab(paste0("PC2 (", pcv[2], "%)")) +
    theme_bw() +
    geom_label_repel(aes(label = sample.name.list), size = 5, show.legend = FALSE) +
    scale_colour_manual(values = drug_colors) +
    theme(axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15)) +
    labs(title = paste0("PCA for ", time_label, " samples"),
         subtitle = "grouped by drug treatment and condition")
  
  return(plot.pca)
}


# run on with different types of groupings
EXP.2h.and.6h <- EXP.coding.matrix
sample.class.2h.and.6h <- sample.class %>% mutate(condition = as.factor(paste(concentrations, timepoints, sep = "_")))
pca_result_2h_and_6h <- run_pca(EXP.2h.and.6h, sample.class.2h.and.6h)
plot.pca.2h.and.6h <- plot_pca(pca_result_2h_and_6h, "2h and 6h")
print(plot.pca.2h.and.6h)

EXP.2h <- EXP.coding.matrix %>% select(gene_id, contains("_2h"))
sample.class.2h <- sample.class %>% filter(grepl("_2h$", sample.name.list)) %>% mutate(condition = as.factor(paste(concentrations, timepoints, sep = "_")))
pca_result_2h <- run_pca(EXP.2h, sample.class.2h)
plot.pca.2h <- plot_pca(pca_result_2h, "2h")
print(plot.pca.2h)

EXP.6h <- EXP.coding.matrix %>% select(gene_id, contains("_6h"))
sample.class.6h <- sample.class %>% filter(grepl("_6h$", sample.name.list)) %>% mutate(condition = as.factor(paste(concentrations, timepoints, sep = "_")))
pca_result_6h <- run_pca(EXP.6h, sample.class.6h)
plot.pca.6h <- plot_pca(pca_result_6h, "6h")
print(plot.pca.6h)

EXP.100nM <- EXP.coding.matrix %>% select(gene_id, matches("(_100nM|_NT_|_DMSO_)"))
sample.class.100nM <- sample.class %>% filter(grepl("(_100nM|_NT_|_DMSO_)", sample.name.list)) %>% mutate(condition = as.factor(paste(concentrations, timepoints, sep = "_")))
pca_result_100nM <- run_pca(EXP.100nM, sample.class.100nM)
plot.pca.100nM <- plot_pca(pca_result_100nM, "100nM")
print(plot.pca.100nM)

EXP.1uM <- EXP.coding.matrix %>% select(gene_id, matches("(_1uM|_NT_|_DMSO_)"))
sample.class.1uM <- sample.class %>% filter(grepl("(_1uM|_NT_|_DMSO_)", sample.name.list)) %>% mutate(condition = as.factor(paste(concentrations, timepoints, sep = "_")))
pca_result_1uM <- run_pca(EXP.1uM, sample.class.1uM)
plot.pca.1uM <- plot_pca(pca_result_1uM, "1uM")
print(plot.pca.1uM)

EXP.PCA_set1_2h <- EXP.coding.matrix %>% select(gene_id, matches("(_NT_2h|_DMSO_2h|_IHK44_.*_2h|_A485_.*_2h|_dCBP_.*_2h)"))
sample.class.set1_2h <- sample.class %>% filter(grepl("(_NT_2h|_DMSO_2h|_IHK44_.*_2h|_A485_.*_2h|_dCBP_.*_2h)", sample.name.list)) %>% mutate(condition = as.factor(paste(concentrations, timepoints, sep = "_")))
pca_result_set1_2h <- run_pca(EXP.PCA_set1_2h, sample.class.set1_2h)
plot.pca.set1_2h <- plot_pca(pca_result_set1_2h, "dual inhibitors/degraders (2h)")
print(plot.pca.set1_2h)

EXP.PCA_set1_6h <- EXP.coding.matrix %>% select(gene_id, matches("(_NT_6h|_DMSO_6h|_IHK44_.*_6h|_A485_.*_6h|_dCBP_.*_6h)"))
sample.class.set1_6h <- sample.class %>% filter(grepl("(_NT_6h|_DMSO_6h|_IHK44_.*_6h|_A485_.*_6h|_dCBP_.*_6h)", sample.name.list)) %>% mutate(condition = as.factor(paste(concentrations, timepoints, sep = "_")))
pca_result_set1_6h <- run_pca(EXP.PCA_set1_6h, sample.class.set1_6h)
plot.pca.set1_6h <- plot_pca(pca_result_set1_6h, "dual inhibitors/degraders (6h)")
print(plot.pca.set1_6h)

EXP.PCA_set2_2h <- EXP.coding.matrix %>% select(gene_id, matches("(_NT_2h|_DMSO_2h|_LS_.*_2h|_QL_.*_2h|_JQAD_.*_2h)"))
sample.class.set2_2h <- sample.class %>% filter(grepl("(_NT_2h|_DMSO_2h|_LS_.*_2h|_QL_.*_2h|_JQAD_.*_2h)", sample.name.list)) %>% mutate(condition = as.factor(paste(concentrations, timepoints, sep = "_")))
pca_result_set2_2h <- run_pca(EXP.PCA_set2_2h, sample.class.set2_2h)
plot.pca.set2_2h <- plot_pca(pca_result_set2_2h, "selective degraders (2h)")
print(plot.pca.set2_2h)

EXP.PCA_set2_6h <- EXP.coding.matrix %>% select(gene_id, matches("(_NT_6h|_DMSO_6h|_LS_.*_6h|_QL_.*_6h|_JQAD_.*_6h)"))
sample.class.set2_6h <- sample.class %>% filter(grepl("(_NT_6h|_DMSO_6h|_LS_.*_6h|_QL_.*_6h|_JQAD_.*_6h)", sample.name.list)) %>% mutate(condition = as.factor(paste(concentrations, timepoints, sep = "_")))
pca_result_set2_6h <- run_pca(EXP.PCA_set2_6h, sample.class.set2_6h)
plot.pca.set2_6h <- plot_pca(pca_result_set2_6h, "selective degraders (6h)")
print(plot.pca.set2_6h)

EXP.PCA_set3 <- EXP.coding.matrix %>% select(gene_id, matches("(_NT_|_DMSO_|_IHK44_.*|_A485_.*|_dCBP_.*)"))
sample.class.set3 <- sample.class %>% filter(grepl("(_NT_|_DMSO_|_IHK44_.*|_A485_.*|_dCBP_.*)", sample.name.list)) %>% mutate(condition = as.factor(paste(concentrations, timepoints, sep = "_")))
pca_result_set3 <- run_pca(EXP.PCA_set3, sample.class.set3)
plot.pca.set3 <- plot_pca(pca_result_set3, "dual inhibitors/degraders")
print(plot.pca.set3)

EXP.PCA_set4 <- EXP.coding.matrix %>% select(gene_id, matches("(_NT_|_DMSO_|_LS_.*|_QL_.*|_JQAD_.*)"))
sample.class.set4 <- sample.class %>% filter(grepl("(_NT_|_DMSO_|_LS_.*|_QL_.*|_JQAD_.*)", sample.name.list)) %>% mutate(condition = as.factor(paste(concentrations, timepoints, sep = "_")))
pca_result_set4 <- run_pca(EXP.PCA_set4, sample.class.set4)
plot.pca.set4 <- plot_pca(pca_result_set4, "selective inhibitors")
print(plot.pca.set4)


















#########################################################################################################
######################################## Log2(TPM + 1) heatmap ##########################################
#########################################################################################################

# REQUIRED: load genelist
setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/GeneSets/List_format/")
gene.list <- read.table("GRYDER_RH4_CR_TFs_CRISPRTop.genelist.txt", sep="\t", header=F)

# REQUIRED: choose samples
samples <- c(
  "RH4_DMSO_6h", 
  "RH4_NT_6h",
  "RH4_JQAD_100nM_6h", 
  "RH4_LS_100nM_6h", 
  "RH4_QL_100nM_6h", 
  "RH4_A485_100nM_6h", 
  "RH4_dCBP_100nM_6h", 
  "RH4_IHK44_100nM_6h"
)

# filtering for samples, minimum TPM, and genes
EXP.coding.matrix.heatmap <- EXP.coding.matrix[ , c("gene_id", samples)]
EXP.coding.matrix.heatmap$maxTPM <- apply(EXP.coding.matrix.heatmap[ , 2:ncol(EXP.coding.matrix.heatmap)], 1, max)
cutoff.expression.min <- 10
EXP.coding.matrix.heatmap <- subset(EXP.coding.matrix.heatmap, maxTPM > cutoff.expression.min)
EXP.coding.matrix.heatmap <- subset(EXP.coding.matrix.heatmap, gene_id %in% gene.list$V1)

# clean up
EXP.coding.matrix.heatmap <- EXP.coding.matrix.heatmap[match(gene.list$V1, EXP.coding.matrix.heatmap$gene_id), ]
rownames(EXP.coding.matrix.heatmap) <- EXP.coding.matrix.heatmap$gene_id
EXP.coding.matrix.heatmap <- EXP.coding.matrix.heatmap[ , 2:(ncol(EXP.coding.matrix.heatmap) - 1)]
EXP.coding.matrix.heatmap <- as.matrix(t(EXP.coding.matrix.heatmap))
rownames(EXP.coding.matrix.heatmap) <- gsub("RH4_", "", rownames(EXP.coding.matrix.heatmap))

# plot
color.scale <- colorRampPalette(c("white", "pink", "red"))(100)
pheatmap(log2(EXP.coding.matrix.heatmap + 1),
         cluster_rows=T,
         cluster_cols=T,
         scale='none',
         show_rownames=T,
         show_colnames=T,
         angle_col=315,
         main="Heatmap of Log2(TPM + 1)",
         color=color.scale)













#########################################################################################################
########################################### Log2FC heatmap ##############################################
#########################################################################################################

# REQUIRED: load genelist
setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/GeneSets/List_format/")
gene.list <- read.table("GRYDER_RH4_CR_TFs_CRISPRTop.genelist.txt", sep="\t", header=F)

# REQUIRED: choose samples
samples.FC <- c(
  "RH4_JQAD_1uM_6h_FC",
  "RH4_LS_1uM_6h_FC",
  "RH4_QL_1uM_6h_FC",
  "RH4_dCBP_1uM_6h_FC",
  "RH4_A485_1uM_6h_FC",
  "RH4_IHK44_1uM_6h_FC"
)

# filtering for samples and genes
EXP.log2FC.heatmap <- EXP.log2FC[ , c("gene_id", samples.FC)]
EXP.log2FC.heatmap <- subset(EXP.log2FC.heatmap, gene_id %in% gene.list$V1)

# clean up
EXP.log2FC.heatmap <- EXP.log2FC.heatmap[match(gene.list$V1, EXP.log2FC.heatmap$gene_id), ]
rownames(EXP.log2FC.heatmap) <- EXP.log2FC.heatmap$gene_id
EXP.log2FC.heatmap <- EXP.log2FC.heatmap[ , -1]
EXP.log2FC.heatmap <- as.matrix(t(EXP.log2FC.heatmap))
rownames(EXP.log2FC.heatmap) <- gsub("RH4_", "", rownames(EXP.log2FC.heatmap))

# plot
color.scale <- colorRampPalette(c("red", "pink", "white"))(100)
pheatmap(EXP.log2FC.heatmap,
         cluster_rows=T,
         cluster_cols=T,
         scale="none",
         show_rownames=T,
         show_colnames=T,
         angle_col=315,
         main="Heatmap of Log2FC",
         color=color.scale)









#########################################################################################################
############################### box/violin plots across multiple genesets ###############################
#########################################################################################################

# REQUIRED: choose samples
samples <- c("IHK44_100nM_2h_FC", 
             "IHK44_1uM_2h_FC", 
             "IHK44_100nM_6h_FC", 
             "IHK44_1uM_6h_FC")

# OPTIONAL: choose color scheme
samples.colors <- c("IHK44_100nM_2h_FC" = "#F6CDCD",
                    "IHK44_100nM_6h_FC" = "#EA4444",
                    "IHK44_1uM_2h_FC" = "#EE7272",
                    "IHK44_1uM_6h_FC" = "#FF0000")

# REQUIRED: load genelists
setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/GeneSets/List_format/")
housekeeping_genes <- read.table("house_keeping_genes.txt", header=F, stringsAsFactors=F) %>% mutate(gene_set="Housekeeping genes")
RH4_CR_TFs <- read.table("GRYDER_RH4_CR_TFs.genelist.txt", header=F, stringsAsFactors=F) %>% mutate(gene_set="RH4 CR TFs")
all_genesets <- bind_rows(housekeeping_genes, RH4_CR_TFs)

# filter for genes and samples
EXP.log2FC.box = EXP.log2FC
colnames(EXP.log2FC.box) <- gsub("^RH4_", "", colnames(EXP.log2FC.box))
rownames(EXP.log2FC.box) <- EXP.log2FC.box$gene_id
EXP.log2FC.box <- EXP.log2FC.box[rownames(EXP.log2FC.box) %in% all_genesets$V1, ]
EXP.log2FC.box <- EXP.log2FC.box[ , colnames(EXP.log2FC.box) %in% samples]
EXP.log2FC.box <- EXP.log2FC.box[, samples]

# clean up
EXP.log2FC.box.plot <- as.data.frame(EXP.log2FC.box)
EXP.log2FC.box.plot <- rownames_to_column(EXP.log2FC.box.plot, var = "gene")
EXP.log2FC.box.plot <- pivot_longer(EXP.log2FC.box.plot, cols = -gene, names_to = "sample_name", values_to = "log2FC")
EXP.log2FC.box.plot <- left_join(EXP.log2FC.box.plot, all_genesets, by = c("gene" = "V1"))
EXP.log2FC.box.plot <- EXP.log2FC.box.plot[, c("gene", "sample_name", "log2FC", "gene_set")]
EXP.log2FC.box.plot$sample_name <- factor(EXP.log2FC.box.plot$sample_name, levels = samples)

# plot
ggplot(EXP.log2FC.box.plot, aes(x=sample_name, y=log2FC, fill=sample_name)) +
  geom_violin(alpha=0.6, scale="width", width=0.5) + 
  geom_boxplot(width=0.25, color="black", alpha=0.7, outlier.shape=NA) +
  facet_wrap(~ gene_set, scales = "fixed") +
  scale_fill_manual(values = samples.colors) +
  labs(x="", y="Log2FC", title="Boxplot of Log2FC for dual inhibitors/degraders") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, hjust=1))








#########################################################################################################
##################################### TPM bar plot for a single gene ####################################
#########################################################################################################

# REQUIRED: choose sample(s)
samples <- c("NT_6h", 
             "DMSO_6h", 
             "IHK44_100nM_2h", 
             "IHK44_1uM_2h", 
             "IHK44_100nM_6h", 
             "IHK44_1uM_6h")

# OPTIONAL: choose color scheme
samples.colors <- c("NT_6h" = "gray60", 
                    "DMSO_6h" = "gray40", 
                    "IHK44_100nM_2h" = "#F6CDCD",
                    "IHK44_100nM_6h" = "#EA4444",
                    "IHK44_1uM_2h" = "#EE7272",
                    "IHK44_1uM_6h" = "#FF0000")

# REQUIRED: choose gene
gene <- "MYCN"

# filter for gene and samples
EXP.coding.matrix.bar <- EXP.coding.matrix[EXP.coding.matrix$gene_id == gene, ]
colnames(EXP.coding.matrix.bar) <- gsub("^RH4_", "", colnames(EXP.coding.matrix.bar))
EXP.coding.matrix.bar <- EXP.coding.matrix.bar[, c("gene_id", samples)]
EXP.coding.matrix.bar <- pivot_longer(EXP.coding.matrix.bar, cols = -gene_id, names_to = "sample_name", values_to = "expression")

# plot
ggplot(EXP.coding.matrix.bar, aes(x=factor(sample_name, levels = samples), y=expression, fill=sample_name)) +
  geom_bar(stat="identity", position="dodge", width=0.8) + 
  scale_fill_manual(values=samples.colors) +
  labs(x="", y="TPM", title=paste("TPM levels for", gene)) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, hjust=1))































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


# REQUIRED: choose gene of interest
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







