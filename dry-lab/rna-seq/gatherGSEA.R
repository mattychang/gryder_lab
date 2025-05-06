##########################
# Updated: 05.06.2025 (Matt)
##########################

library(plyr)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(grid)
library(gridExtra)
library(tibble)

GSEA.folder <- "/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC/GSEA_results/GSEA_results_against_DMSO"
GSEA.folder.name <- "GSEA_results_against_DMSO"
project.folder <- "/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC/GSEA_results"

# parse sample information
# note: do NOT use "." in when exporting GSEA results --> will cause incorrect parsing
setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC/GSEA_results")
samples <- list.dirs(path = GSEA.folder, full.names = F, recursive = F)
sample.df <- read.table(text = samples, sep=".")
sample.df$V3 <- as.character(sample.df$V3)
colnames(sample.df) <- c("Drug","GSEA_Style","Digits")
sample.df$neg.xls.path <- paste(GSEA.folder, "/", paste(sample.df$Drug, sample.df$GSEA_Style, sample.df$Digits, sep="."), "/gsea_report_for_na_neg_", sample.df$Digits,".tsv", sep="")
sample.df$pos.xls.path <- paste(GSEA.folder, "/", paste(sample.df$Drug, sample.df$GSEA_Style, sample.df$Digits, sep="."), "/gsea_report_for_na_pos_", sample.df$Digits, ".tsv", sep="")

# import raw GSEA results
column.names <- c("NAME", "ES", "NES", "NOM.p.val", "FDR.q.val", "FWER.p.val", "Drug")
allsamples <- data.frame(matrix(ncol = length(column.names), nrow = 0))
colnames(allsamples) <- column.names
lapply(sample.df$Drug, function(x) {
  
    temp.df = subset(sample.df,sample.df$Drug %in% x)
    
    if(file.exists(temp.df$neg.xls.path) == 'TRUE') {
      df.neg = read.table(temp.df$neg.xls.path,sep="\t",header=T)
    }
    else {
      df.neg = df.empty
    }
    if(file.exists(temp.df$pos.xls.path) == 'TRUE') {
      df.pos = read.table(temp.df$pos.xls.path,sep="\t",header=T)
    }
    else {
      df.pos = df.empty
    }
    
    df = rbind(df.pos,df.neg)                                                                                     # combine df.neg and df.pos
    df = df[,c("NAME","ES","NES","NOM.p.val","FDR.q.val","FWER.p.val")]                                           # subselect
    df$Drug = x                                                                                                   # add Drug column for identification
    
    allsamples <<- rbind(allsamples, df)
})

# clean up and transform
row.names(allsamples) <- 1:nrow(allsamples)
allsamples$NES <- as.numeric(allsamples$NES)
allsamples$NOM.p.val <- as.numeric(allsamples$NOM.p.val) + 0.00001
allsamples <- na.omit(allsamples)
allsamples$log10.NOM.p.val <- log10(allsamples$NOM.p.val)
    
# order genesets by mean NES
NAME.ranks <- aggregate(NES ~ NAME, allsamples, mean)
colnames(NAME.ranks) <- c("NAME", "rankmetric")
allsamples <- join(allsamples, NAME.ranks, by = "NAME")
allsamples$NAME <- factor(allsamples$NAME, levels = allsamples[order(unique(allsamples$rankmetric, decreasing=F)),]$NAME)

ES.min <- min(abs(allsamples$ES))
NES.min <- min(abs(allsamples$NES))

# OPTIONAL: simplify naming convention
allsamples$NAME <- gsub(allsamples$NAME, pattern = "PAX3FOXO1", replacement = "P3F")
allsamples$NAME <- gsub(allsamples$NAME, pattern = "ENHANCERS", replacement = "ENH")
allsamples$NAME <- gsub(allsamples$NAME, pattern = "DIFFERENTIATION", replacement = "DIFF")

# OPTIONAL: export compiled/curated GSEA results
write.table(allsamples, paste(GSEA.folder, ".GSEA.allsamples.summary.txt", sep = ""),
                              col.names = T,
                              row.names = F,
                              quote = F,
                              sep="\t")











########################## bubble charts and heatmaps ##########################

# REQUIRED: choose samples
included_samples <- c("IHK44_100nM_6h", "A485_100nM_6h")
GENE.plot.subset <- subset(allsamples, allsamples$Drug %in% included_samples) 

# bubble plot 1
ggplot(GENE.plot.subset, aes(x = Drug, y = NAME, size = abs(ES) - ES.min, colour = ES)) + 
  geom_point() +
  scale_colour_gradient2(low = "tomato",mid = "white", high = "lightskyblue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# bubble plot 2
ggplot(GENE.plot.subset, aes(x = Drug, y = NAME, size = -log10.NOM.p.val, colour = ES)) + 
  geom_point() +
  scale_colour_gradient2(low = "blue",mid = "white", high = "orange") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# heatmap
heatmap.matrix <- GENE.plot.subset %>%
  select(NAME, NES, Drug) %>%
  pivot_wider(names_from = NAME, values_from = NES) %>%
  column_to_rownames("Drug") %>%
  as.matrix()
pheatmap(heatmap.matrix,
         scale = "row", 
         cluster_rows = T, 
         cluster_cols = T,
         main = paste(GSEA.folder.name," heatmap (metric: NES)", sep=""),
         legend = T)








########################## enrichment score and density plots ##########################

# REQUIRED: choose geneset
Geneset = "GRYDER_RH4_TOP_CRTFS"

# get rank and enrichment score for all genes in all samples for Geneset
Geneset.allsamples = data.frame(RANK.IN.GENE.LIST = numeric(), 
                                RUNNING.ES = numeric(), 
                                SYMBOL = character(), 
                                Drug = character())
lapply(sample.df$Drug, function(x) {
  ESplot.sample = subset(sample.df, sample.df$Drug %in% x)
  ESplot.sample$geneset.path = paste(GSEA.folder, 
                                     "/", 
                                     paste(ESplot.sample$Drug, ESplot.sample$GSEA_Style, ESplot.sample$Digits, sep="."),
                                     "/", 
                                     Geneset, 
                                     ".tsv", 
                                     sep="")
  Geneset.df = read.table(ESplot.sample$geneset.path, sep = "\t", header = T)
  Geneset.df = Geneset.df[,c("RANK.IN.GENE.LIST", "RUNNING.ES", "SYMBOL")]
  Geneset.df$Drug = x
  Geneset.allsamples <<- rbind(Geneset.allsamples, Geneset.df)
})


# REQUIRED: choose sample(s)
included_samples <- c("IHK44_100nM_6h")
Geneset.allsamples <- subset(Geneset.allsamples, Geneset.allsamples$Drug %in% included_samples)
drugcolors <- c(
  "A485_100nM_2h" = "blue", "A485_100nM_6h" = "blue", "A485_1uM_2h" = "blue", "A485_1uM_6h" = "blue", 
  "dCBP_100nM_2h" = "green", "dCBP_100nM_6h" = "green", "dCBP_1uM_2h" = "green", "dCBP_1uM_6h" = "green", 
  "IHK44_100nM_2h" = "red", "IHK44_100nM_6h" = "red", "IHK44_1uM_2h" = "red", "IHK44_1uM_6h" = "red", 
  "JQAD_100nM_2h" = "yellow", "JQAD_100nM_6h" = "yellow", "JQAD_1uM_2h" = "yellow", "JQAD_1uM_6h" = "yellow", 
  "LS_100nM_2h" = "purple", "LS_100nM_6h" = "purple", "LS_1uM_2h" = "purple", "LS_1uM_6h" = "purple", 
  "QL_100nM_2h" = "orange", "QL_100nM_6h" = "orange", "QL_1uM_2h" = "orange", "QL_1uM_6h" = "orange"
)

# OPTIONAL: choose gene(s) of interest
#Geneset.allsamples.pick = subset(Geneset.allsamples, Geneset.allsamples$SYMBOL %in% c("MYCN", "MYOG", "MYOD1"))

# plot for enrichment score
p1 <- ggplot(Geneset.allsamples, aes(x = RANK.IN.GENE.LIST, y = RUNNING.ES, group = Drug, color = Drug)) + 
  geom_line(size = 1.2) + 
  geom_hline(yintercept = 0, lty = 'dashed') +
  scale_color_manual(values = drugcolors) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.25, 0.25)) +
  ylab("Enrichment Score") + 
  xlab(paste("Genes Ranked by log2 fold change")) +
  ggtitle(paste("Geneset:", Geneset, sep=" ")) #+
  #geom_point(data = Geneset.allsamples.pick, color = "black") +
  #geom_text(data = Geneset.allsamples.pick, label = Geneset.allsamples.pick$SYMBOL)

# plot for density
p2 <- ggplot(Geneset.allsamples,group = Drug) + 
  stat_density(aes(x = RANK.IN.GENE.LIST, y = 0.5, fill = ..density..), geom = "tile", position="identity") + 
  facet_wrap(~Drug, ncol=1) + 
  scale_fill_gradient(low = "white", high = "darkgrey") +
  theme(axis.line = element_blank(), 
        axis.title.x = element_blank(), 
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background=element_blank()
        )

# combine enrichment score and density plots
p3 <- p2 + geom_linerange(aes(x = RANK.IN.GENE.LIST, ymin = 0, ymax = 1, color = Drug)) +
  facet_wrap(~Drug, ncol=1) +
  theme_bw() + 
  scale_color_manual(values = drugcolors) +
  theme(axis.line = element_blank(),
        axis.title.x = element_blank(), 
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background=element_blank()
        )

# display plots
grid.arrange(arrangeGrob(p1,ncol=1))
grid.arrange(arrangeGrob(p1,p3,ncol=1))



# OPTIONAL: export
export_path <- file.path(project.folder, paste0("RH4_CRTFs.GSEA_matrix.txt"))
write.table(Geneset.allsamples, file = export_path, sep = "\t", row.names = F, col.names = T, quote = F, append = F)
