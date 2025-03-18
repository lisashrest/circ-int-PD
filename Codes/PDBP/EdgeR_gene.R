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

setwd("data/")

#import data from ciriquant
lib_mtx <- read.csv("library_info_PDBP_final.csv", row.names = 1)
gene_mtx <- read.csv("gene_count_matrix.csv", row.names = 1, check.names=FALSE)
gene_mtx <- gene_mtx[ , rownames(lib_mtx)]
bsj_mtx <- read.csv("circRNA_bsj_PDBP_final.csv", row.names = 1, check.names=FALSE)
details <- read.csv("meta/sample_details_outliers_removed.csv", stringsAsFactors = TRUE, header = TRUE)
details$filename <- paste0(details$Sample.Participant, "_", details$Visit)
rownames(details) <- details$filename

details$Group <- factor(details$Enroll.Case.Control,
                        levels = c("Control", "Case"),
                        labels = c("Control", "Parkinson"))

details$HY.Stage <- factor(details$HY.Stage)

details$Condition <- factor(paste0(details$Group, ".", details$EVENT_ID))

#Subset the original matrix column
gene_mtx_up <- gene_mtx[, details$filename]

#Create EdgeR Object
gene_DGE <- DGEList(counts = gene_mtx_up, group = details$Condition)

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
logCPM <- cpm(gene_DGE_norm, log=TRUE)

plotMDS(logCPM, col = as.numeric(details$Lab.Plate))
plotMDS(logCPM, col = as.numeric(details$Enroll.Sex))

#Creating Design Matrix
design_factor <- as.data.frame(details)

PATNO <- factor(design_factor$Sample.Participant)

EVENT_ID <- factor(design_factor$EVENT_ID,
                   levels = c("BL", "V06"))

Group <- factor(design_factor$Enroll.Case.Control,
                levels = c("Control", "Case"))

Sex <- factor(design_factor$Enroll.Sex, 
              levels = c("Male","Female"))

design <- model.matrix(~ PATNO + Sex)

HC.V06 <- Group=="Control" & EVENT_ID=="V06" 
PD.V06 <- Group=="Case" & EVENT_ID=="V06"

design<-cbind(design,HC.V06, PD.V06)
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
V06vsBL.C_res_gene <- glmQLFTest(fit_gene, coef = "HC.V06") #Year 2 vs Baseline

V06vsBL.P_is.gene <- decideTestsDGE(V06vsBL.P_res_gene, adjust.method = "BH", p.value = 0.05)

V06vsBL.C_is.gene <- decideTestsDGE(V06vsBL.C_res_gene, adjust.method = "BH", p.value = 0.05)

summary(V06vsBL.P_is.gene)
summary(V06vsBL.C_is.gene)

#MA plots PD groups
plotMD(V06vsBL.P_res_gene, status=V06vsBL.P_is.gene, main = "DE Genes: Parkinson's Participants Samples at Year 2 vs Baseline (Paired)")+ 
  abline(h=c(-0.2,0.2), col="red")

#MA plots Control group
plotMD(V06vsBL.C_res_gene, status=V06vsBL.C_is.gene, main = "DE Genes: Control Participants Samples at Year 2 vs Baseline (Paired)")+ 
  abline(h=c(-0.2,0.2), col="red")

#Check the top DEGs in PD participants
#circ
setwd("~/Desktop/CircRNA-Project/Post_review/EdgeR/results/PDBP_gene/")

topc_V06 <- as.data.frame(topTags(V06vsBL.P_res_gene, n = Inf))
#write.csv(topc_V06, file = "DE_Analysis_P_PD_V06vsBL_gene_final.csv")

#Check the top DEGs in HC participants
#circ
topc_V06_HC <- as.data.frame(topTags(V06vsBL.C_res_gene, n = Inf))
#write.csv(topc_V06_HC, file = "DE_Analysis_C_V06vsBL_gene_final.csv")

