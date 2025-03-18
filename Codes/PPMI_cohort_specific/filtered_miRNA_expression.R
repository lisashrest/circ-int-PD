library(edgeR)
setwd("miRNA_smRNAseq/edger_qc/")

#import data from ciriquant
details <- read.csv("meta/sample_details_1.csv", stringsAsFactors = FALSE, header = TRUE)
details$filename <- paste0("P1_", details$PATNO, "_", details$EVENT_ID)
rownames(details) <- details$filename

#Converting variables of interests into factors
details$Group <- factor(details$APPRDX, 
                        levels = c(2,1),
                        labels = c("Control", "Parkinson"))

#Creating Object List 
details_up <- subset(details, details$genetic == "LRRK2-/SNCA-/GBA-")
details_up$filename <- paste0("P1_", details_up$PATNO, "_", details_up$EVENT_ID)

miRNA_mtx <- read.csv("mature_counts.csv", header =  TRUE, row.names = 1)
miRNA_mtx <- t(miRNA_mtx)

#Import seq metadata
meta <- read.csv("meta/metaDataIR3.csv", header = TRUE)
head(meta)

meta <- meta %>%
  mutate(filename = paste0("P1_", PATNO, "_", CLINICAL_EVENT))

#Extract information from meta file
meta_up <- inner_join(details_up, meta, by = "filename")
meta_up$Plate <- factor(meta_up$Plate)
meta_up$GENDER <- factor(meta_up$GENDER)
meta_up$Condition <- factor(paste0(meta_up$Group, ".", meta_up$EVENT_ID))

#Columns of interest
column <- c("PATNO.x", "Group", "Condition", "EVENT_ID", "GENDER", "Plate", "age", "filename")
meta_final <- meta_up[, column]
rownames(meta_final) <- meta_final$filename

#Subset the original matrix column
miRNA_mtx_up <- miRNA_mtx[, meta_final$filename]
dim(miRNA_mtx_up)

#CreatingDGEList
miRNA_DGE <- DGEList(counts = miRNA_mtx_up, group = meta_final$Condition)

#Checking library sizes
head(miRNA_DGE$samples)

#Filtering low reads
#circ Filtering
#AverageCPM 
AveLogCPM <- aveLogCPM(miRNA_DGE)
hist(AveLogCPM)

#Library Size Normalization
miRNA_DGE_norm <- calcNormFactors(miRNA_DGE, method = "TMM")

#Cluster the samples together 
logCPM <- cpm(miRNA_DGE_norm, log=TRUE)
#write.csv(logCPM, "~/Desktop/CircRNA-Project/Post_review/circRNA-sponging_pipeline/filtered_expression/PPMI_miRNA_logCPM.csv", row.names = TRUE)

CPM <- cpm(miRNA_DGE_norm, log = FALSE)
#write.csv(CPM, "~/Desktop/CircRNA-Project/Post_review/circRNA-sponging_pipeline/filtered_expression/PPMI_miRNA_CPM.csv", row.names = TRUE)
