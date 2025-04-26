# Last edit: 04.02.2025 (Matt)

# import libraries
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(dplyr)
library(tidyr)
library(VennDiagram)

getwd()
# load summary data. summary data lets NA = 0 --> misleading
#data_100nM = read.csv("/Users/matthewchang/Library/CloudStorage/GoogleDrive-msc148@case.edu/.shortcut-targets-by-id/1jekwTrBV3zPxqzw9vlL6ZGaqT_qMi-Mk/GryderDrive/Projects/ChemBio/Sub-projects/CBPp300i/TALUS_proteomics/Data/240313_DMSOvs100nM_IHK44_stats_table.csv")
#data_150nM = read.csv("/Users/matthewchang/Library/CloudStorage/GoogleDrive-msc148@case.edu/.shortcut-targets-by-id/1jekwTrBV3zPxqzw9vlL6ZGaqT_qMi-Mk/GryderDrive/Projects/ChemBio/Sub-projects/CBPp300i/TALUS_proteomics/Data/240313_DMSOvs150nM_IHK44_stats_table.csv")

# load replicate data
data_replicates = read.csv("/Users/matthewchang/Library/CloudStorage/GoogleDrive-msc148@case.edu/.shortcut-targets-by-id/1jekwTrBV3zPxqzw9vlL6ZGaqT_qMi-Mk/GryderDrive/Projects/ChemBio/Sub-projects/CBPp300i/TALUS_proteomics/TALUS_Data/gryder_protein_matrix_QvalFil.csv")

# export path
export_file_path = "/Users/matthewchang/Library/CloudStorage/GoogleDrive-msc148@case.edu/.shortcut-targets-by-id/1jekwTrBV3zPxqzw9vlL6ZGaqT_qMi-Mk/GryderDrive/Projects/ChemBio/Sub-projects/CBPp300i/TALUS_proteomics"







# =========================== pre-processing of data_replicates ===========================

# remove proteins that do not have at least 3 replicates
DMSO_cols = 2:13
IHK44_100nM_cols = 14:19
IHK44_125nM_cols = 20:25
IHK44_150nM_cols = 26:31
data_replicates_filtered = data_replicates %>%
  rowwise() %>%
  filter(
    sum(!is.na(c_across(all_of(DMSO_cols)))) >= 3,
    sum(!is.na(c_across(all_of(IHK44_100nM_cols)))) >= 3,
    sum(!is.na(c_across(all_of(IHK44_125nM_cols)))) >= 3,
    sum(!is.na(c_across(all_of(IHK44_150nM_cols)))) >= 3
  ) %>%
  ungroup()

# remove dBET columns
dBET_cols = 32:49
data_replicates_filtered = data_replicates_filtered %>% select(-all_of(dBET_cols))

# clean up
colnames(data_replicates_filtered) = sub("^X240301_RH4_", "", colnames(data_replicates_filtered))

# create summary table
data_replicates_filtered_statistics = data_replicates_filtered %>%
  rowwise() %>%
  mutate(
    # extract replicate values without NA
    DMSO_values = list(na.omit(c_across(all_of(DMSO_cols)))),
    IHK44_100nM_values = list(na.omit(c_across(all_of(IHK44_100nM_cols)))),
    IHK44_125nM_values = list(na.omit(c_across(all_of(IHK44_125nM_cols)))),
    IHK44_150nM_values = list(na.omit(c_across(all_of(IHK44_150nM_cols)))),
    
    # mean calculations
    DMSO_mean = mean(DMSO_values),
    IHK44_100nM_mean = mean(IHK44_100nM_values),
    IHK44_125nM_mean = mean(IHK44_125nM_values),
    IHK44_150nM_mean = mean(IHK44_150nM_values),
    
    # SEM calculations
    DMSO_SEM = if(length(DMSO_values) > 1) sd(DMSO_values) / sqrt(length(DMSO_values)) else NA_real_,
    IHK44_100nM_SEM = if(length(IHK44_100nM_values) > 1) sd(IHK44_100nM_values) / sqrt(length(IHK44_100nM_values)) else NA_real_,
    IHK44_125nM_SEM = if(length(IHK44_125nM_values) > 1) sd(IHK44_125nM_values) / sqrt(length(IHK44_125nM_values)) else NA_real_,
    IHK44_150nM_SEM = if(length(IHK44_150nM_values) > 1) sd(IHK44_150nM_values) / sqrt(length(IHK44_150nM_values)) else NA_real_,
    
    # Log2FC calculation
    log2FC_100nM = log2(IHK44_100nM_mean + 1e-6) - log2(DMSO_mean + 1e-6),
    log2FC_125nM = log2(IHK44_125nM_mean + 1e-6) - log2(DMSO_mean + 1e-6),
    log2FC_150nM = log2(IHK44_150nM_mean + 1e-6) - log2(DMSO_mean + 1e-6),
    
    # t-test calculations
    pval_100nM = if(length(DMSO_values) >= 2 && length(IHK44_100nM_values) >= 2){
      t.test(DMSO_values, IHK44_100nM_values)$p.value
    } else NA_real_,
    
    pval_125nM = if(length(DMSO_values) >= 2 && length(IHK44_125nM_values) >= 2){
      t.test(DMSO_values, IHK44_125nM_values)$p.value
    } else NA_real_,
    
    pval_150nM = if(length(DMSO_values) >= 2 && length(IHK44_150nM_values) >= 2){
      t.test(DMSO_values, IHK44_150nM_values)$p.value
    } else NA_real_,
    
    # -log10(p-value) calculation
    negLog10Pval_100nM = -log10(pmax(pval_100nM, 1e-300)),
    negLog10Pval_125nM = -log10(pmax(pval_125nM, 1e-300)),
    negLog10Pval_150nM = -log10(pmax(pval_150nM, 1e-300))
  ) %>%
  ungroup() %>%
  select(
    Protein.Names,
    DMSO_mean, DMSO_SEM,
    IHK44_100nM_mean, IHK44_100nM_SEM, log2FC_100nM, negLog10Pval_100nM,
    IHK44_125nM_mean, IHK44_125nM_SEM, log2FC_125nM, negLog10Pval_125nM,
    IHK44_150nM_mean, IHK44_150nM_SEM, log2FC_150nM, negLog10Pval_150nM
  )







write.table(data_replicates_filtered,
            file = file.path(export_file_path, "data_replicates_filtered.txt"),
            sep = "\t",
            quote = F,
            row.names = F)

write.table(data_replicates_filtered_statistics,
            file = file.path(export_file_path, "data_replicates_filtered_statistics.txt"),
            sep = "\t",
            quote = F,
            row.names = F)



# get removed proteins
data_replicates_removed = data_replicates %>%
  rowwise() %>%
  filter(
    sum(!is.na(c_across(all_of(DMSO_cols)))) < 3 |
      sum(!is.na(c_across(all_of(IHK44_100nM_cols)))) < 3 |
      sum(!is.na(c_across(all_of(IHK44_125nM_cols)))) < 3 |
      sum(!is.na(c_across(all_of(IHK44_150nM_cols)))) < 3
  ) %>%
  ungroup()

data_replicates_removed = data_replicates_removed %>% select(-all_of(dBET_cols))

write.table(data_replicates_removed,
            file = file.path(export_file_path, "data_replicates_removed.txt"),
            sep = "\t",
            quote = F,
            row.names = F)













# =========================== Volcano Plot ===========================

# define conditions to plot for (100nM, 125nM, 150nM)
negLog10Pval_col = "negLog10Pval_100nM"
log2FC_col = "log2FC_100nM"
title_label = "DMSO vs IHK-44 100nM"

# define proteins of interest
proteins_to_label = c("MYOD1_HUMAN", "MYOG_HUMAN", "SOX8_HUMAN", "PAX3_HUMAN", "FOXO1_HUMAN", "JUN_HUMAN")

# significance cutoffs
log2FC_cutoff_high = 1
log2FC_cutoff_low = 0.5
negLog10Pval_cutoff = 1.3

# classification
data_volcano = data_replicates_filtered_statistics %>%
  mutate(
    negLog10Pval = .data[[negLog10Pval_col]],
    log2_fold_change = .data[[log2FC_col]],
    Significance_vol = case_when(
      negLog10Pval > negLog10Pval_cutoff & log2_fold_change > log2FC_cutoff_high ~ "Upregulated (High)",
      negLog10Pval > negLog10Pval_cutoff & log2_fold_change < -log2FC_cutoff_high ~ "Downregulated (High)",
      negLog10Pval > negLog10Pval_cutoff & log2_fold_change > log2FC_cutoff_low & log2_fold_change <= log2FC_cutoff_high ~ "Upregulated (Small)",
      negLog10Pval > negLog10Pval_cutoff & log2_fold_change < -log2FC_cutoff_low & log2_fold_change >= -log2FC_cutoff_high ~ "Downregulated (Small)",
      TRUE ~ "Low"
    ),
    Proteins_of_Interest = ifelse(Protein.Names %in% proteins_to_label, "Yes", "No")
  )


# labeling
data_volcano_labels = data_volcano %>%
  #filter(Significance_vol %in% c("Upregulated (High)", "Downregulated (High)"))
  filter(Proteins_of_Interest == "Yes")

# plot
volcano_plot = ggplot(data_volcano, aes(x = log2_fold_change, y = negLog10Pval, color = Significance_vol)) +
  geom_point(size = 1.3) +
  scale_color_manual(values = c(
    "Low" = "gray",
    "Upregulated (Small)" = "cyan",
    "Downregulated (Small)" = "orange",
    "Upregulated (High)" = "blue",
    "Downregulated (High)" = "red"
  )) +
  geom_hline(yintercept = negLog10Pval_cutoff, linetype = "dashed", color = "gray10", size = 0.5) +
  geom_vline(xintercept = c(-log2FC_cutoff_high, log2FC_cutoff_high), linetype = "dashed", color = "red", size = 0.5) +
  geom_vline(xintercept = c(-log2FC_cutoff_low, log2FC_cutoff_low), linetype = "dashed", color = "darkorange", size = 0.5) +
  theme_minimal() +
  labs(title = title_label, x = "Log2 Fold Change", y = "-log10(p-value)") +
  theme(legend.title = element_blank()) +
  geom_label_repel(
    data = data_volcano_labels,
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

# create lists of proteins based on "Significance_vol"
upregulated_high = data_volcano %>% filter(Significance_vol == "Upregulated (High)")
downregulated_high = data_volcano %>% filter(Significance_vol == "Downregulated (High)")
upregulated_intermediate = data_volcano %>% filter(Significance_vol == "Upregulated (Small)")
downregulated_intermediate = data_volcano %>% filter(Significance_vol == "Downregulated (Small)")

cat("Upregulated (High) Proteins:\n")
print(upregulated_high)
cat("\nDownregulated (High) Proteins:\n")
print(downregulated_high)
cat("\nUpregulated (Small) Proteins:\n")
print(upregulated_intermediate)
cat("\nDownregulated (Small) Proteins:\n")
print(downregulated_intermediate)















# =========================== Bar Plot ===========================

# define protein of interest
proteins_to_label = "CBP_HUMAN"

# filtering
data_bar = data_replicates_filtered_statistics %>%
  filter(Protein.Names == proteins_to_label) %>%
  select(Protein.Names, 
         DMSO_mean, DMSO_SEM, 
         IHK44_100nM_mean, IHK44_100nM_SEM,
         IHK44_125nM_mean, IHK44_125nM_SEM,
         IHK44_150nM_mean, IHK44_150nM_SEM) %>%
  pivot_longer(cols = c(DMSO_mean, IHK44_100nM_mean, IHK44_125nM_mean, IHK44_150nM_mean),
               names_to = "Condition",
               values_to = "Intensity") %>%
  mutate(
    SEM = case_when(
      Condition == "DMSO_mean" ~ DMSO_SEM,
      Condition == "IHK44_100nM_mean" ~ IHK44_100nM_SEM,
      Condition == "IHK44_125nM_mean" ~ IHK44_125nM_SEM,
      Condition == "IHK44_150nM_mean" ~ IHK44_150nM_SEM
    ),
    Condition = factor(Condition,
                       levels = c("DMSO_mean", "IHK44_100nM_mean", "IHK44_125nM_mean", "IHK44_150nM_mean"),
                       labels = c("DMSO", "IHK44 100nM", "IHK44 125nM", "IHK44 150nM"))
  )

# plot
bar_chart = ggplot(data_bar, aes(x = Condition, y = Intensity, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = Intensity - SEM, ymax = Intensity + SEM), width = 0.2, position = position_dodge(0.8)) +
  theme_minimal() +
  labs(title = paste("Protein Intensity for", proteins_to_label),
       x = "Condition",
       y = "Protein Intensity") +
  scale_fill_manual(values = c("DMSO" = "gray", 
                               "IHK44 100nM" = "red1", 
                               "IHK44 125nM" = "red3", 
                               "IHK44 150nM" = "red4"))
print(bar_chart)







# =========================== Box Plot ===========================

# define protein of interest
proteins_to_label = "CBP_HUMAN"

# filtering
DMSO_cols_names = colnames(data_replicates_filtered)[2:13]
IHK44_100nM_cols_names = colnames(data_replicates_filtered)[14:19]
IHK44_125nM_cols_names = colnames(data_replicates_filtered)[20:25]
IHK44_150nM_cols_names = colnames(data_replicates_filtered)[26:31]
data_box = data_replicates_filtered %>%
  filter(Protein.Names == proteins_to_label) %>%
  select(Protein.Names, all_of(c(DMSO_cols_names, IHK44_100nM_cols_names, IHK44_125nM_cols_names, IHK44_150nM_cols_names))) %>%
  pivot_longer(
    cols = -Protein.Names,
    names_to = "Sample",
    values_to = "Intensity"
  ) %>%
  mutate(
    Condition = case_when(
      Sample %in% DMSO_cols_names ~ "DMSO",
      Sample %in% IHK44_100nM_cols_names ~ "IHK44 100nM",
      Sample %in% IHK44_125nM_cols_names ~ "IHK44 125nM",
      Sample %in% IHK44_150nM_cols_names ~ "IHK44 150nM",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Intensity), !is.na(Condition))

# plot
box_violin_plot = ggplot(data_box, aes(x = Condition, y = Intensity, fill = Condition)) +
  #geom_violin(alpha = 0.4, trim = FALSE, scale = "width", width = 0.5) +
  geom_boxplot(width = 0.5, alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.05, size = 2, shape = 21, color = "black", alpha = 0.9) +
  theme_minimal() +
  labs(title = paste("Protein Intensity for", proteins_to_label),
       x = "Condition",
       y = "Protein Intensity") +
  scale_fill_manual(values = c("DMSO" = "gray", 
                               "IHK44 100nM" = "red1", 
                               "IHK44 125nM" = "red3", 
                               "IHK44 150nM" = "red4"))
print(box_violin_plot)










# =========================== MA Plot (not updated) ===========================

# define significance cutoffs
log2FC_cutoff_high = 1
log2FC_cutoff_low = 0.5
negLog10Pval_cutoff = 1.3 

# proteins of interest
proteins_to_label = c("MYOD1_HUMAN", "MYOG_HUMAN", "SOX8_HUMAN", "PAX3_HUMAN", "FOXO1_HUMAN", "JUN_HUMAN")

# compute "A" abundance, and (M) Log2FC
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
  # geom_hline(yintercept = c(-log2FC_cutoff_high, log2FC_cutoff_high), linetype = "dashed", color = "red", size = 0.5) +  # Threshold lines
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
    arrow = arrow(length = unit(0.01, "npc"), type = "closed"),
    max.overlaps = Inf
  )

print(ma_plot)






# =========================== Scatter Plot (not updated) ===========================

# define significance cutoffs
log2FC_cutoff_high = 1
log2FC_cutoff_low = 0.5
negLog10Pval_cutoff = 1.3

# proteins of interest
proteins_to_label = c("MYOD1_HUMAN", "MYOG_HUMAN", "SOX8_HUMAN", "PAX3_HUMAN", "FOXO1_HUMAN", "JUN_HUMAN")

# compute Log10 of DMSO and IHK44 expression
data_100nM_cor = data_100nM %>%
  mutate(
    log10_DMSO = log10(DMSO_mean + 1),
    log10_IHK44 = log10(IHK44_100nM_mean + 1),
    Proteins_of_Interest = ifelse(Protein.Names %in% proteins_to_label, "Yes", "No")
  )

# apply significance cutoffs
data_100nM_cor = data_100nM_cor %>%
  mutate(
    Significance_MA = case_when(
      negLog10Pval > negLog10Pval_cutoff & log2_fold_change > log2FC_cutoff_high ~ "Upregulated (High)",
      negLog10Pval > negLog10Pval_cutoff & log2_fold_change < -log2FC_cutoff_high ~ "Downregulated (High)",
      negLog10Pval > negLog10Pval_cutoff & abs(log2_fold_change) > log2FC_cutoff_low & abs(log2_fold_change) <= log2FC_cutoff_high ~ "Intermediate",
      TRUE ~ "Low"
    )
  )

# fit linear regression model and residuals
lm_fit = lm(log10_IHK44 ~ log10_DMSO, data = data_100nM_cor)

data_100nM_cor = data_100nM_cor %>%
  mutate(
    predicted_IHK44 = predict(lm_fit, newdata = data_100nM_cor),
    residual = log10_IHK44 - predicted_IHK44
  )

residual_sd = sd(data_100nM_cor$residual, na.rm = TRUE)
outlier_threshold = 10 * residual_sd
data_filtered_for_R2 = data_100nM_cor %>%
  filter(abs(residual) <= outlier_threshold)
correlation_filtered = cor(data_filtered_for_R2$log10_DMSO, data_filtered_for_R2$log10_IHK44)
R2_filtered = correlation_filtered^2

# filter (for labeling)
data_100nM_scatter_labels = data_100nM_cor %>%
  filter(Significance_MA %in% c("Upregulated (High)", "Downregulated (High)"))
  # filter(Proteins_of_Interest == "Yes")

# plot
scatter_plot = ggplot(data_100nM_cor, aes(x = log10_DMSO, y = log10_IHK44, color = Significance_MA)) +
  geom_point(size = 1.3, alpha = 0.8) +
  scale_color_manual(values = c(
    "Low" = "gray",
    "Intermediate" = "darkorange",
    "Upregulated (High)" = "blue",
    "Downregulated (High)" = "red"
  )) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  theme_minimal() +
  xlim(1, max(data_100nM_cor$log10_DMSO)) +
  ylim(1, max(data_100nM_cor$log10_IHK44)) +
  labs(
    title = "Scatter Plot: DMSO vs IHK-44 100nM",
    x = "log10(DMSO Intensity)",
    y = "log10(IHK-44 Intensity)"
  ) +
  theme(legend.title = element_blank()) +
  
  geom_label_repel(
    data = data_100nM_scatter_labels,
    aes(label = Protein.Names),
    size = 3,
    color = "black",
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "black",
    segment.size = 0.5,
    arrow = arrow(length = unit(0.01, "npc"), type = "closed"),
    max.overlaps = Inf
  ) +
  
  annotate("text", x = 1.5,  
           y = max(data_100nM_cor$log10_IHK44) - 0.5, 
           label = paste0("R² (10×SD) = ", round(R2_filtered, 4)), 
           size = 5, hjust = 0)
print(scatter_plot)

# extreme outliers (not detected)
extreme_outliers = data_100nM_cor %>%
  filter(abs(residual) > outlier_threshold) %>%
  select(Protein.Names, log10_DMSO, log10_IHK44, predicted_IHK44, residual)
cat("Extreme Outliers (Residual > 10 * SD):\n")
print(extreme_outliers)









# =========================== Venn diagram (not updated) ===========================
# Replace these with your actual file paths
file1 = "/Users/matthewchang/Library/CloudStorage/GoogleDrive-msc148@case.edu/.shortcut-targets-by-id/1jekwTrBV3zPxqzw9vlL6ZGaqT_qMi-Mk/GryderDrive/Projects/ChemBio/Sub-projects/CBPp300i/TALUS_proteomics/downregulated_high_IHK44_100nM.txt"
file2 = "/Users/matthewchang/Library/CloudStorage/GoogleDrive-msc148@case.edu/.shortcut-targets-by-id/1jekwTrBV3zPxqzw9vlL6ZGaqT_qMi-Mk/GryderDrive/Projects/ChemBio/Sub-projects/CBPp300i/TALUS_proteomics/downregulated_high_IHK44_150nM.txt"

# Read the data
data1 = read.table(file1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
data2 = read.table(file2, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Extract unique protein names
set1_proteins = unique(data1$Protein.Names)
set2_proteins = unique(data2$Protein.Names)

venn.plot = venn.diagram(
  x = list(IHK44_100nM = set1_proteins, IHK44_150nM = set2_proteins),
  filename = NULL,  # prevents writing to file, draws in RStudio viewer
  fill = c("skyblue", "salmon"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  #cat.pos = c(0, 180),        # adjust angle of labels
  cat.dist = c(0.07, 0.07),   # push labels away from circles
  cat.just = list(c(1, 0), c(0, 0)),  # left/right align
  main = "Significant and downregulated proteins"
)

# Plot it
grid::grid.newpage()  # clears previous plot before drawing a new one
grid::grid.draw(venn.plot)



