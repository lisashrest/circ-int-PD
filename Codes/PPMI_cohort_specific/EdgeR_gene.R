##########For gene expression analysis
#BiocManager::install("edgeR")
#BiocManager::install("org.Hs.eg.db")
#install.packages("dplyr")
#install.packages("plyr")
#install.packages("tidyverse")
#install.packages("splitstackshape")
#devtools::install_github('kevinblighe/EnhancedVolcano')

library(splitstackshape)
library(plyr)
library(tidyverse)
library(config)
library(edgeR)
library(org.Hs.eg.db)
library(dplyr)
library(plyr)
library(tidyverse)
library(ggrepel)


setwd("~/Desktop/CircRNA-Project/PPMI/data/")

#import data from ciriquant
lib_mtx <- read.csv("library_info.csv", row.names = 1)
gene_mtx <- read.csv("gene_count_matrix.csv", row.names = 1, check.names=FALSE)
gene_mtx <- gene_mtx[ , rownames(lib_mtx)]
bsj_mtx <- read.csv("circRNA_bsj.csv", row.names = 1, check.names=FALSE)
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

#Import seq metadata
meta <- read.csv("metaDataIR3.csv", header = TRUE)
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
gene_mtx_up <- gene_mtx[, meta_final$filename]
bsj_mtx_up <- bsj_mtx[, meta_final$filename]

gene_DGE <- DGEList(counts = gene_mtx_up, group = meta_final$Condition)
#Checking library sizes
head(gene_DGE$samples)

#Filtering low reads
#circ Filtering
keep <- filterByExpr(gene_DGE)
gene_DGE <- gene_DGE[keep, ]
dim(gene_DGE)

#AverageCPM 
AveLogCPM_gene <- aveLogCPM(gene_DGE)
hist(AveLogCPM_gene)

#Library Size Normalization
gene_DGE_norm <- calcNormFactors(gene_DGE, method = "TMM")

#Cluster the samples together 
CPM <- cpm(gene_DGE_norm, log = FALSE)
logCPM <- cpm(gene_DGE_norm, log=TRUE)
plotMDS(logCPM, col = as.numeric(meta_up$GENDER))

#Creating Design Matrix
design_factor <- as.data.frame(meta_final)

PATNO <- factor(design_factor$PATNO.x)

EVENT_ID <- factor(design_factor$EVENT_ID)

Group <- factor(design_factor$Group)

Sex <- factor(design_factor$GENDER, 
              levels = c("Male","Female"))

design <- model.matrix(~ PATNO + Sex)

HC.V06 <- Group=="Control" & EVENT_ID=="V06" 
PD.V06 <- Group=="Parkinson" & EVENT_ID=="V06"
PD.V08 <- Group=="Parkinson" & EVENT_ID=="V08"
HC.V08 <- Group=="Control" & EVENT_ID=="V08"

design<-cbind(design,HC.V06, HC.V08, PD.V06, PD.V08)
colnames(design)
is.fullrank(design)

design <- design[,!grepl("SexFemale", colnames(design))] # getting to full rank.
is.fullrank(design)

#Estimating Dispersion
gene_DGE <- estimateDisp(gene_DGE_norm, design, robust= TRUE)
plotBCV(gene_DGE)

#Estimating the QL dispersion 
fit_gene <- glmQLFit(gene_DGE, design, robust=TRUE)
head(fit_gene$coefficients)

#Differential Expression Analysis
V06vsBL.P_res_gene <- glmQLFTest(fit_gene, coef = "PD.V06") #Year 2 vs Baseline
V08vsBL.P_res_gene <- glmQLFTest(fit_gene, coef = "PD.V08") #Year 3 vs Baseline
V08vsV06.P_res_gene <- glmQLFTest(fit_gene, contrast =  c(rep(0,142), -1, 1)) #Year 3 vs Year 2

V06vsBL.C_res_gene <- glmQLFTest(fit_gene, coef = "HC.V06") #Year 2 vs Baseline
V08vsBL.C_res_gene <- glmQLFTest(fit_gene, coef = "HC.V08") #Year 3 vs Baseline
V08vsV06.C_res_gene <- glmQLFTest(fit_gene, contrast =  c(rep(0,140), -1, 1, 0, 0)) #Year 3 vs Year 2

V06vsBL.P_is.gene <- decideTestsDGE(V06vsBL.P_res_gene, adjust.method = "BH", p.value = 0.05)
V08vsBL.P_is.gene <- decideTestsDGE(V08vsBL.P_res_gene, adjust.method = "BH", p.value = 0.05)
V08vsV06.P_is.gene <- decideTestsDGE(V08vsV06.P_res_gene, adjust.method = "BH", p.value = 0.05)

V06vsBL.C_is.gene <- decideTestsDGE(V06vsBL.C_res_gene, adjust.method = "BH", p.value = 0.05)
V08vsBL.C_is.gene <- decideTestsDGE(V08vsBL.C_res_gene, adjust.method = "BH", p.value = 0.05)
V08vsV06.C_is.gene <- decideTestsDGE(V08vsV06.C_res_gene, adjust.method = "BH", p.value = 0.05)

summary(V06vsBL.P_is.gene)
summary(V08vsBL.P_is.gene)
summary(V08vsV06.P_is.gene)

summary(V06vsBL.C_is.gene)
summary(V08vsBL.C_is.gene)
summary(V08vsV06.C_is.gene)

#MA plots PD groups
plotMD(V06vsBL.P_res_gene, status=V06vsBL.P_is.gene, main = "DE Genes: Parkinson's Participants Samples at Year 2 vs Baseline (Paired)")+ 
  abline(h=c(-0.2,0.2), col="red")

plotMD(V08vsBL.P_res_gene, status=V08vsBL.P_is.gene, main = "DE Genes: Parkinson's Participants Samples at Year 3 vs Baseline (Paired)")+ 
  abline(h=c(-0.2,0.2), col="red")

plotMD(V08vsV06.P_res_gene, status=V08vsV06.P_is.gene, main = "DE Genes: Parkinson's Participants Samples at Year 3 vs Year 2 (Paired)")+ 
  abline(h=c(-0.2,0.2), col="red")

#MA plots Control group
plotMD(V06vsBL.C_res_gene, status=V06vsBL.C_is.gene, main = "DE Genes: Control Participants Samples at Year 2 vs Baseline (Paired)")+ 
  abline(h=c(-0.2,0.2), col="red")

plotMD(V08vsBL.C_res_gene, status=V08vsBL.C_is.gene, main = "DE Genes: Control Participants Samples at Year 3 vs Baseline (Paired)")+ 
  abline(h=c(-0.2,0.2), col="red")

plotMD(V08vsV06.C_res_gene, status=V08vsV06.C_is.gene, main = "DE Genes: Control Participants Samples at Year 3 vs Year 2 (Paired)")+ 
  abline(h=c(-0.2,0.2), col="red")

#Check the top DEGs in PD participants
#circ
topg_V06 <- as.data.frame(topTags(V06vsBL.P_res_gene, n = Inf))
#write.csv(topg_V06, file = "~/Desktop/CircRNA-Project/Post_review/EdgeR/results/PPMI_gene/DE_Analysis_P_PD_V06vsBL_gene_final.csv")

topg_V08 <- as.data.frame(topTags(V08vsBL.P_res_gene, n= Inf))
#write.csv(topg_V08, file = "~/Desktop/CircRNA-Project/Post_review/EdgeR/results/PPMI_gene/DE_Analysis_P_PD_V08vsBL_gene_final.csv")

topg_V06_HC <- as.data.frame(topTags(V06vsBL.C_res_gene, n = Inf))
#write.csv(topg_V06_HC, file = "~/Desktop/CircRNA-Project/Post_review/EdgeR/results/DE_Analysis_C_V06vsBL_gene_final.csv")

topg_V08_HC <- as.data.frame(topTags(V08vsBL.C_res_gene, n= Inf))
#write.csv(topg_V08_HC, file = "~/Desktop/CircRNA-Project/Post_review/EdgeR/results/DE_Analysis_C_V08vsBL_gene_final.csv")

topg_V06_V08_PD <- as.data.frame(topTags(V08vsV06.P_res_gene, n=Inf))
#write.csv(topg_V06_V08_PD, file = "~/Desktop/CircRNA-Project/Post_review/EdgeR/results/PPMI_gene/DE_Analysis_P_V08vsV06_gene_final.csv")

topg_V06_V08_HC <- as.data.frame(topTags(V08vsV06.C_res_gene, n=Inf))
#write.csv(topg_V06_V08_PD, file = "~/Desktop/CircRNA-Project/Post_review/EdgeR/results/PPMI_gene/DE_Analysis_C_V08vsV06_gene_final.csv")
