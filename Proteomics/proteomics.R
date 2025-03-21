# Last edit: 03.21.2025 (Matt)

# import libraries
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(dplyr)
library(tidyr)

# load data
data_100nM = read.csv("/Users/matthewchang/Documents/R/Proteomics_data/240313_DMSOvs100nM_IHK44_stats_table.csv")
data_150nM = read.csv("/Users/matthewchang/Documents/R/Proteomics_data/240313_DMSOvs150nM_IHK44_stats_table.csv")

# =========================== Volcano Plot ===========================

# define significance cutoffs
log2FC_cutoff_high = 1
log2FC_cutoff_low = 0.5
negLog10Pval_cutoff = 1.3

# proteins of interest
proteins_to_label = c("MYOD1_HUMAN", "MYOG_HUMAN", "SOX8_HUMAN", "PAX3_HUMAN", "FOXO1_HUMAN", "JUN_HUMAN")

# add "Significance" to assign proteins that are overly expressed and "Proteins_of_Interest"
data_100nM_vol = data_100nM %>%
  mutate(
    Significance_vol = case_when(
      negLog10Pval > negLog10Pval_cutoff & log2_fold_change > log2FC_cutoff_high ~ "Upregulated (High)",
      negLog10Pval > negLog10Pval_cutoff & log2_fold_change < -log2FC_cutoff_high ~ "Downregulated (High)",
      negLog10Pval > negLog10Pval_cutoff & abs(log2_fold_change) > log2FC_cutoff_low & abs(log2_fold_change) <= log2FC_cutoff_high ~ "Intermediate",
      TRUE ~ "Low"
    ),
    Proteins_of_Interest = ifelse(Protein.Names %in% proteins_to_label, "Yes", "No")
  )

# filter (for labeling)
data_100nM_vol_labels = data_100nM_vol %>%
  filter(Significance_vol %in% c("Upregulated (High)", "Downregulated (High)"))
         # Proteins_of_Interest == "Yes")

# plot
volcano_plot = ggplot(data_100nM_vol, aes(x = log2_fold_change, y = negLog10Pval, color = Significance_vol)) +
  geom_point(size = 1.3) +
  scale_color_manual(values = c(
    "Low" = "gray",
    "Intermediate" = "darkorange",
    "Upregulated (High)" = "blue",
    "Downregulated (High)" = "red"
  )) +
  geom_hline(yintercept = negLog10Pval_cutoff, linetype = "dashed", color = "gray10", size = 0.5) +
  geom_vline(xintercept = c(-log2FC_cutoff_high, log2FC_cutoff_high), linetype = "dashed", color = "red", size = 0.5) +
  geom_vline(xintercept = c(-log2FC_cutoff_low, log2FC_cutoff_low), linetype = "dashed", color = "darkorange", size = 0.5) +
  theme_minimal() +
  labs(title = "Volcano Plot: DMSO vs IHK-44 100 nM", x = "Log2 Fold Change", y = "-log10(p-value)") +
  theme(legend.title = element_blank()) +
  
  geom_label_repel(
    data = data_100nM_vol_labels,
    aes(label = Protein.Names),
    size = 3,
    color = "black",
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "black",
    segment.size = 0.5,
    arrow = arrow(length = unit(0.01, "npc"), type = "closed"),
    max.overlaps = Inf
  )

print(volcano_plot)

# create lists of overly expressed proteins
upregulated_high = data_100nM_vol %>% filter(Significance_vol == "Upregulated (High)")
downregulated_high = data_100nM_vol %>% filter(Significance_vol == "Downregulated (High)")
intermediate = data_100nM_vol %>% filter(Significance_vol == "Intermediate")

# print
cat("Upregulated (High) Proteins:\n")
print(upregulated_high)

cat("Downregulated (High) Proteins:\n")
print(downregulated_high)

cat("Intermediate Proteins:\n")
print(intermediate)













# =========================== Bar Plot ===========================

# protein of interest
protein_of_interest = "EGR1_HUMAN"

# filter data for the specified protein
protein_data = data_100nM %>%
  filter(Protein.Names == protein_of_interest) %>%
  select(Protein.Names, DMSO_mean, IHK44_100nM_mean) %>%
  pivot_longer(cols = c(DMSO_mean, IHK44_100nM_mean), names_to = "Condition", values_to = "Intensity") %>%
  mutate(Condition = factor(Condition, levels = c("DMSO_mean", "IHK44_100nM_mean"), labels = c("DMSO", "IHK44 100 nM")))

# plot
bar_chart = ggplot(protein_data, aes(x = Condition, y = Intensity, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = paste("Protein Intensity for", protein_of_interest),
       x = "Condition",
       y = "Protein Intensity") +
  scale_fill_manual(values = c("DMSO" = "gray", "IHK44 100 nM" = "red"))

print(bar_chart)










# =========================== MA Plot ===========================

# define significance cutoffs
log2FC_cutoff_high = 1
log2FC_cutoff_low = 0.5
negLog10Pval_cutoff = 1.3 

# proteins of interest
proteins_to_label = c("MYOD1_HUMAN", "MYOG_HUMAN", "SOX8_HUMAN", "PAX3_HUMAN", "FOXO1_HUMAN", "JUN_HUMAN")

# compute "A" the average abundance, and (M) the log2 fold change
data_100nM_ma = data_100nM %>%
  mutate(
    A = (log2(DMSO_mean + 1) + log2(IHK44_100nM_mean + 1)) / 2,
    M = log2_fold_change
  )

# assign significance
data_100nM_ma = data_100nM_ma %>%
  mutate(
    Significance_MA = case_when(
      negLog10Pval > negLog10Pval_cutoff & log2_fold_change > log2FC_cutoff_high ~ "Upregulated (High)",
      negLog10Pval > negLog10Pval_cutoff & log2_fold_change < -log2FC_cutoff_high ~ "Downregulated (High)",
      negLog10Pval > negLog10Pval_cutoff & abs(log2_fold_change) > log2FC_cutoff_low & abs(log2_fold_change) <= log2FC_cutoff_high ~ "Intermediate",
      TRUE ~ "Low"
    ),
    Proteins_of_Interest = ifelse(Protein.Names %in% proteins_to_label, "Yes", "No")
  )

# filter (for labeling)
data_100nM_ma_label = data_100nM_ma %>%
  filter(Significance_MA %in% c("Upregulated (High)", "Downregulated (High)"))
  # filter(Proteins_of_Interest == "Yes")

# plot
ma_plot = ggplot(data_100nM_ma, aes(x = A, y = M, color = Significance_MA)) +
  geom_point(size = 1.3) +
  scale_color_manual(values = c(
    "Low" = "gray",
    "Intermediate" = "darkorange",
    "Upregulated (High)" = "blue",
    "Downregulated (High)" = "red"
  )) +
  geom_hline(yintercept = c(-log2FC_cutoff_high, log2FC_cutoff_high), linetype = "dashed", color = "red", size = 0.5) +  # Threshold lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", size = 0.5) +
  theme_minimal() +
  labs(
    title = "MA Plot: DMSO vs IHK-44 100 nM",
    x = "Average log2 Intensity (A)",
    y = "Log2 Fold Change (M)"
  ) +
  theme(legend.title = element_blank()) +
  
  geom_label_repel(
    data = data_100nM_ma_label,
    aes(label = Protein.Names),
    size = 3,
    color = "black",
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "black",
    segment.size = 0.5,
