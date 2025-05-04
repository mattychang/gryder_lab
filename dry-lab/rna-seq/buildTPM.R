##########################
# Updated: 02.24.24 (Matt)
##########################
# This script makes the TPM matrix, Log2FC matrix, and GSEA ranklists
##########################

getwd()
setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC/SampleList")                       # define the sample list path
sample.list.all = read.table(file.choose(), header = T, sep = "\t")                                             # acquire the sample list (.txt file, "Sample" as first entry)
sample.list = sample.list.all$Sample

project.folder = "/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC"                        # define the project folder path
sample.set = "IHK_samples"                                                                                      # save sample name

##########################
# create TPM matrix
##########################

setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/")                                                               # define the data path to retrieve samples
file.exists = file.exists(paste0("DATA/", sample.list,"/", sample.list, ".gene.TPM.txt", sep=""))               # check to see if TPM.txt exists for each sample (TRUE for all existing files)
sample.list = sample.list[file.exists]                                                                          # remove samples whose data file does not exist

lapply(sample.list, function(x) {                                                                               # make coding TPM file, normalize based on expected counts
    coding <- read.table("ref/RSEM_hg38/HGNC_protein-coding_gene_19229_2022.txt", sep="\t", header=T)           # load in more stringent protein coding list (from Diana)
    EXP <- read.table(paste("DATA/",x,"/",x,".genes.results",sep=""), sep="\t", header=T)                       # load RSEM (genes.results) output file for the current sample with expected count, TPM, FPKM
    EXP$coding = EXP$gene_id %in% coding$symbol                                                                 # remove non-coding RNA entries
    EXP.coding <<- subset(EXP, EXP$coding %in% c("TRUE"))
    expected_sum = sum(EXP.coding$expected_count)                                                               # sum all coding expected counts
    EXP.coding$count_norm = EXP.coding$expected_count / expected_sum * 1000000                                  # normalize expected counts to counts per million
    write.table(EXP.coding[,1:9],
                file=paste("DATA/",x,"/",x,".gene.coding.norm.txt",sep=""), 
                sep="\t", 
                row.names=FALSE, 
                col.names=TRUE, 
                quote=FALSE)
})                                                                                                              # note: outputs should all be NULL

EXP.coding.matrix = EXP.coding["gene_id"]                                                                       # create coding matrix
lapply(sample.list, function(x) {                                                                               
    EXP.sample = read.table(paste("DATA/",x,"/",x,".gene.coding.norm.txt",sep=""), sep="\t", header=T)          # read the normalized counts
    EXP.sample = as.data.frame(EXP.sample[,9])                                                                  # extract the normalized counts
    removable.string = "Sample_"
    sample.name = gsub(removable.string,"",x)
    colnames(EXP.sample) = c(sample.name)
    EXP.coding.matrix <<- cbind(EXP.coding.matrix, EXP.sample)                                                  # combine the current sample data into the matrix
})

# clean up
colnames(EXP.coding.matrix) = gsub("-", "", colnames(EXP.coding.matrix))
colnames(EXP.coding.matrix) = gsub("_RNA_022924_CWRU", "", colnames(EXP.coding.matrix))

# export
dir.create(file.path(project.folder, "ExpMatrices"))
file_path = file.path(project.folder, "ExpMatrices", paste0(sample.set, ".coding.norm.matrix.txt"))
write.table(EXP.coding.matrix, 
            file = file_path, 
            row.names = F, 
            col.names = T, 
            quote = F)







##########################
# create Log2FC matrix
##########################

EXP.log2FC = EXP.coding.matrix

# Log2FC calculation step
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

# clean up
EXP.log2FC = EXP.log2FC[, -c(2:33)]
EXP.log2FC[,(ncol(EXP.coding.matrix) + 1) : ncol(EXP.log2FC)] = round(EXP.log2FC[,(ncol(EXP.coding.matrix) + 1) : ncol(EXP.log2FC)], digits = 5)

# export
file_path = file.path(project.folder, "ExpMatrices", paste0(sample.set, ".log2FC.txt"))
write.table(EXP.log2FC, 
            file = file_path, 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE)


##########################
# create GSEA ranklists
##########################

EXP.GSEA = EXP.log2FC

dir.create(file.path(project.folder, "GSEA_ranklist", "vs_DMSO"), recursive = T)

# write out rank (.rnk) files for GSEA
for (i in 2:ncol(EXP.GSEA)) {
  Ranklist = data.frame(EXP.GSEA[, 1])                                                                            # create Ranklist df, 1st column is gene_id
  Ranklist$DeltaTPM = EXP.GSEA[, i]                                                                               # 2nd column is DeltaTPM
  
  Ranklist = Ranklist[rev(order(Ranklist$DeltaTPM)), ]                                                            # sort in descending order by DeltaTPM
  
  SampleName = colnames(EXP.GSEA)[i]                                                                              # acquire name
  mytime = format(Sys.time(), "%b_%d_%Y")                                                                         # acquire date
  
  myfile = file.path(project.folder, 
                     "GSEA_ranklist", 
                     "vs_DMSO", 
                     paste0(SampleName, "_", mytime, ".rnk"))
  write.table(Ranklist, 
              file = myfile, 
              sep = "\t", 
              row.names = F, 
              col.names = F, 
              quote = F, 
              append = F)
}

###  GSEA instructions: use the following datasets and settings
###  Under "Tools", click "Run GSEApreranked"
###  Gene sets database: GRYDERLAB_RMS_gene_sets.gmt
###  Collapse/Remap to gene symbols: No_collapse
###  Chip platoform: Human_Gene_Symbol_with_Remapping_MSigDBv2023.2.Hs.chip
###  Max size: 5000
###  Min size: 10





