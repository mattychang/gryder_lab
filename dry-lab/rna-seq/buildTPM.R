###################################
# Matt, 05.07.2025
# build TPM matrix, log2FC matrix, GSEA ranklists
###################################

project.folder = "/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC"
sample.set = "IHK_samples"






#########################################################################################################
############################## create TPM matrix (coding genes, normalized) #############################
#########################################################################################################

# load sample list (TXT)
setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC/SampleList")
sample.list.all <- read.table("sample-list-IHK-name.txt", header=T, sep="\t")
sample.list <- sample.list.all$Sample
#file.exists = file.exists(paste0("DATA/", sample.list,"/", sample.list, ".gene.TPM.txt", sep=""))
#print(file.exists)

# load RSEM files for all samples and calculate TPM
setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/")
lapply(sample.list, function(x) {
  
  # load more stringent protein coding list (from Diana)
  coding <- read.table("ref/RSEM_hg38/HGNC_protein-coding_gene_19229_2022.txt", sep="\t", header=T)
  
  # load RSEM (genes.results) output with expected count, TPM, FPKM
  EXP <- read.table(paste("DATA/",x,"/",x,".genes.results",sep=""), sep="\t", header=T)
  
  # remove non-coding RNA entries
  EXP$coding <- EXP$gene_id %in% coding$symbol
  EXP.coding <<- subset(EXP, EXP$coding %in% c("TRUE"))
  
  # add up all coding expected counts
  expected_sum <- sum(EXP.coding$expected_count)
  
  # make new column for expected/total
  EXP.coding$count_norm <- EXP.coding$expected_count / expected_sum * 1000000
  
  # export
  write.table(EXP.coding[,1:9], file=paste("DATA/",x,"/",x,".gene.coding.norm.txt",sep=""), sep="\t", row.names=F, col.names=T, quote=F)
})

# build TPM matrix
EXP.coding.matrix <- EXP.coding["gene_id"]
lapply(sample.list, function(x) {
  
  # load TPM values
  EXP.sample <- read.table(paste("DATA/",x,"/",x,".gene.coding.norm.txt",sep=""), sep="\t", header=T)
  EXP.sample <- as.data.frame(EXP.sample[,9])
  
  # clean up (change as needed)
  sample.name <- gsub("Sample_","",x)
  sample.name <- gsub("-", "", sample.name)
  sample.name <- gsub("_RNA_022924_CWRU", "", sample.name)
  colnames(EXP.sample) <- c(sample.name)
  
  # append into TPM matrix
  EXP.coding.matrix <<- cbind(EXP.coding.matrix, EXP.sample)
})

# export
#dir.create(file.path(project.folder, "ExpMatrices"))
file_path <- file.path(project.folder, "ExpMatrices", paste0(sample.set, ".coding.norm.matrix.txt"))
write.table(EXP.coding.matrix, file=file_path, row.names=F, col.names=T, quote=F)









#########################################################################################################
########################################## create Log2FC matrix #########################################
#########################################################################################################

EXP.log2FC = EXP.coding.matrix

# Log2FC calculation
EXP.log2FC$RH4_A485_100nM_2h_FC = log2(EXP.log2FC$RH4_A485_100nM_2h + 1) - log2(EXP.log2FC$RH4_DMSO_2h + 1)
EXP.log2FC$RH4_A485_1uM_2h_FC = log2(EXP.log2FC$RH4_A485_1uM_2h + 1) - log2(EXP.log2FC$RH4_DMSO_2h + 1)
EXP.log2FC$RH4_A485_100nM_6h_FC = log2(EXP.log2FC$RH4_A485_100nM_6h + 1) - log2(EXP.log2FC$RH4_DMSO_6h + 1)
EXP.log2FC$RH4_A485_1uM_6h_FC = log2(EXP.log2FC$RH4_A485_1uM_6h + 1) - log2(EXP.log2FC$RH4_DMSO_6h + 1)
EXP.log2FC$RH4_JQAD_100nM_2h_FC = log2(EXP.log2FC$RH4_JQAD_100nM_2h + 1) - log2(EXP.log2FC$RH4_DMSO_2h + 1)
EXP.log2FC$RH4_JQAD_1uM_2h_FC = log2(EXP.log2FC$RH4_JQAD_1uM_2h + 1) - log2(EXP.log2FC$RH4_DMSO_2h + 1)
EXP.log2FC$RH4_JQAD_100nM_6h_FC = log2(EXP.log2FC$RH4_JQAD_100nM_6h + 1) - log2(EXP.log2FC$RH4_DMSO_6h + 1)
EXP.log2FC$RH4_JQAD_1uM_6h_FC = log2(EXP.log2FC$RH4_JQAD_1uM_6h + 1) - log2(EXP.log2FC$RH4_DMSO_6h + 1)
EXP.log2FC$RH4_dCBP_100nM_2h_FC = log2(EXP.log2FC$RH4_dCBP_100nM_2h + 1) - log2(EXP.log2FC$RH4_DMSO_2h + 1)
EXP.log2FC$RH4_dCBP_1uM_2h_FC = log2(EXP.log2FC$RH4_dCBP_1uM_2h + 1) - log2(EXP.log2FC$RH4_DMSO_2h + 1)
EXP.log2FC$RH4_dCBP_100nM_6h_FC = log2(EXP.log2FC$RH4_dCBP_100nM_6h + 1) - log2(EXP.log2FC$RH4_DMSO_6h + 1)
EXP.log2FC$RH4_dCBP_1uM_6h_FC = log2(EXP.log2FC$RH4_dCBP_1uM_6h + 1) - log2(EXP.log2FC$RH4_DMSO_6h + 1)
EXP.log2FC$RH4_QL_100nM_2h_FC = log2(EXP.log2FC$RH4_QL_100nM_2h + 1) - log2(EXP.log2FC$RH4_DMSO_2h + 1)
EXP.log2FC$RH4_QL_1uM_2h_FC = log2(EXP.log2FC$RH4_QL_1uM_2h + 1) - log2(EXP.log2FC$RH4_DMSO_2h + 1)
EXP.log2FC$RH4_QL_100nM_6h_FC = log2(EXP.log2FC$RH4_QL_100nM_6h + 1) - log2(EXP.log2FC$RH4_DMSO_6h + 1)
EXP.log2FC$RH4_QL_1uM_6h_FC = log2(EXP.log2FC$RH4_QL_1uM_6h + 1) - log2(EXP.log2FC$RH4_DMSO_6h + 1)
EXP.log2FC$RH4_LS_100nM_2h_FC = log2(EXP.log2FC$RH4_LS_100nM_2h + 1) - log2(EXP.log2FC$RH4_DMSO_2h + 1)
EXP.log2FC$RH4_LS_1uM_2h_FC = log2(EXP.log2FC$RH4_LS_1uM_2h + 1) - log2(EXP.log2FC$RH4_DMSO_2h + 1)
EXP.log2FC$RH4_LS_100nM_6h_FC = log2(EXP.log2FC$RH4_LS_100nM_6h + 1) - log2(EXP.log2FC$RH4_DMSO_6h + 1)
EXP.log2FC$RH4_LS_1uM_6h_FC = log2(EXP.log2FC$RH4_LS_1uM_6h + 1) - log2(EXP.log2FC$RH4_DMSO_6h + 1)
EXP.log2FC$RH4_IHK44_100nM_2h_FC = log2(EXP.log2FC$RH4_IHK44_100nM_2h + 1) - log2(EXP.log2FC$RH4_DMSO_2h + 1)
EXP.log2FC$RH4_IHK44_1uM_2h_FC = log2(EXP.log2FC$RH4_IHK44_1uM_2h + 1) - log2(EXP.log2FC$RH4_DMSO_2h + 1)
EXP.log2FC$RH4_IHK44_100nM_6h_FC = log2(EXP.log2FC$RH4_IHK44_100nM_6h + 1) - log2(EXP.log2FC$RH4_DMSO_6h + 1)
EXP.log2FC$RH4_IHK44_1uM_6h_FC = log2(EXP.log2FC$RH4_IHK44_1uM_6h + 1) - log2(EXP.log2FC$RH4_DMSO_6h + 1)
#EXP.log2FC$RH4_IHK45_100nM_2h_FC = log2(EXP.log2FC$RH4_IHK45_100nM_2h + 1) - log2(EXP.log2FC$RH4_DMSO_2h + 1)
#EXP.log2FC$RH4_IHK45_1uM_2h_FC = log2(EXP.log2FC$RH4_IHK45_1uM_2h + 1) - log2(EXP.log2FC$RH4_DMSO_2h + 1)
#EXP.log2FC$RH4_IHK45_100nM_6h_FC = log2(EXP.log2FC$RH4_IHK45_100nM_6h + 1) - log2(EXP.log2FC$RH4_DMSO_6h + 1)
#EXP.log2FC$RH4_IHK45_1uM_6h_FC = log2(EXP.log2FC$RH4_IHK45_1uM_6h + 1) - log2(EXP.log2FC$RH4_DMSO_6h + 1)



# clean up
EXP.log2FC <- EXP.log2FC[ , -c(2:33)]
EXP.log2FC[, 2:25] <- round(EXP.log2FC[, 2:25], digits=5)

# export
file_path <- file.path(project.folder, "ExpMatrices", paste0(sample.set, ".log2FC.txt"))
write.table(EXP.log2FC, file=file_path, row.names=F, col.names=T, quote=F)


#########################################################################################################
######################################### create GSEA ranklists #########################################
#########################################################################################################

EXP.GSEA <- EXP.log2FC

#dir.create(file.path(project.folder, "GSEA_ranklist", "vs_DMSO"), recursive=T)

# write out RNK files for each sample to run GSEA
for (i in 2:ncol(EXP.GSEA)) {
  
  # create df with two columns (gene_id, DeltaTPM)
  Ranklist <- data.frame(EXP.GSEA[, 1])
  Ranklist$DeltaTPM <- EXP.GSEA[, i]
  
  # sort in descending order
  Ranklist <- Ranklist[rev(order(Ranklist$DeltaTPM)), ]
  
  # export
  SampleName <- colnames(EXP.GSEA)[i]
  mytime <- format(Sys.time(), "%b_%d_%Y")
  myfile <- file.path(project.folder, "GSEA_ranklist", "vs_DMSO", paste0(SampleName, "_", mytime, ".rnk"))
  write.table(Ranklist, file=myfile, sep="\t", row.names=F, col.names=F, quote=F, append=F)
}

###  GSEA instructions: use the following datasets and settings
###  Under "Tools", click "Run GSEApreranked"
###  Gene sets database: GRYDERLAB_RMS_gene_sets.gmt
###  Collapse/Remap to gene symbols: No_collapse
###  Chip platoform: Human_Gene_Symbol_with_Remapping_MSigDBv2023.2.Hs.chip
###  Max size: 5000
###  Min size: 10





