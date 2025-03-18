#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.18")
#install.packages("hoardeR")

library(AnnotationDbi)
library(hoardeR)

circ_counts <- read.csv("results/filtered_expression/PPMI_circ_logCPM.csv", header = TRUE)
miRNA_counts <- read.csv("results/filtered_expression/PPMI_miRNA_logCPM.csv", header = TRUE)

gene_counts <- read.csv("results/filtered_expression/PPMI_gene_logCPM.csv", header = TRUE)
gene_counts <- gene_counts %>%
  separate(col = X, into = c("ID", "Version"), sep = "_")%>%
  separate(col = Version, into = c("Version", "Gene"), sep = "\\|")

correlation_pairs <- read.table("results/filtered_expression/filtered_circRNA_miRNA_correlation.tsv", sep = "\t", header = TRUE)

#list of miRNA targets from the correlation pairs
cor_final <- correlation_pairs %>%
  dplyr::filter(adj_pval < 0.05 & spearman_R < 0)

miRNA_list <- cor_final$miRNA
miRNA_list <- miRNA_list[!duplicated(miRNA_list)]
mis_anno <- c("hsa-miR-103a-2-5p", "hsa-miR-103a-1-5p","hsa-miR-6881-3p", "hsa-miR-103a-3p", "hsa-miR-181a-5p", "hsa-miR-16-5p", 
              "hsa-miR-130a-5p", "hsa-miR-454-3p", "hsa-miR-3529-3p", "hsa-miR-5585-3p", "hsa-miR-548av-3p", "hsa-miR-19a-3p", 
              "hsa-miR-30a-3p", "hsa-miR-16-2-3p", "hsa-miR-145-5p", "hsa-miR-30e-3p", "hsa-miR-15a-5p", "hsa-miR-142-3p",
              "hsa-miR-320a-3p", "hsa-miR-122b-3p", "hsa-miR-133a-3p", "hsa-miR-183-5p", "hsa-miR-320a-5p")

miRNA_list_final <- miRNA_list[!miRNA_list%in%mis_anno]

#Use targetscan to get the list of resulting genes 
targets <- list()

for(i in 1:length(miRNA_list_final)){
  targets[[i]] <- targetScan(mirna = miRNA_list_final[i], species = "Human", release = "7.2")
  targets[[i]]['miRNA'] <- miRNA_list_final[i]
}

target_df_final <- do.call(rbind, targets)
target_df_final$`Target Gene` <- target_df_final$Ortholog
target_df_final$interaction <- paste0(target_df_final$miRNA,".", target_df_final$`Target Gene`)
target_df_final <- target_df_final[!duplicated(target_df_final$interaction), ]
#write.csv(target_df_final, file = "results/target_genes_targetscan.csv", row.names = FALSE)

#Using miRTarBase for target miRNAs
miRTarBase_human <- readxl::read_xlsx("hsa_MTI.xlsx")
miRTarBase_human$interaction <- paste0(miRTarBase_human$miRNA, ".", miRTarBase_human$`Target Gene`)
miRTarBase_final <- miRTarBase_human[miRTarBase_human$miRNA%in%miRNA_list, ]

#Get list of miRNA-gene targets common in both miRTarBase and targetscan
#miRNA_target_df <- left_join(target_df_final, miRTarBase_human, by = "interaction")
miRNA_target_df <- full_join(target_df_final, miRTarBase_final, by = "interaction")

#Removing the predicted target by targetscan with missing information from miRTarBase
target_final <- miRNA_target_df[!is.na(miRNA_target_df$`miRTarBase ID`)&!duplicated(miRNA_target_df$interaction), ]
#write.csv(target_final, "results/target_genes_targetscan_miRTarBase.csv", row.names = FALSE)

#Creating a target scan symbol table for SPONGE software
targets <- read.csv("results/target_genes_targetscan_miRTarBase.csv", header = TRUE)

#Creating a frequency table 
targets_mtx <- as.matrix(table(targets$Target.Gene.y, targets$miRNA.y))
print(targets_mtx)
head(targets_mtx)
dim(targets_mtx)

rownames(targets_mtx)
colnames(targets_mtx)
#write.csv(targets_mtx, file = "results/target_genes_mtx.csv")

#Target matrix 
#BiocManager::install("org.Hs.eg.db")
target_mtx <- read.csv("results/target_genes_mtx.csv")
dim(target_mtx)

library(AnnotationDbi)
library(org.Hs.eg.db)

target_mtx$ENSEMBL <- mapIds(org.Hs.eg.db, keys = target_mtx$X,  column = "ENSEMBL", keytype = "SYMBOL") 
target_mtx <- target_mtx[!is.na(target_mtx$ENSEMBL), ]
rownames(target_mtx) <- target_mtx$ENSEMBL
dim(target_mtx)

#write.csv(target_mtx, file = "results/target_genes_mtx_ENSEMBL.csv")
