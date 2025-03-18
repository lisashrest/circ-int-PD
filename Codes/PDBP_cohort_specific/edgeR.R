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
bsj_mtx <- read.csv("sample_details_outliers_removed.csv", stringsAsFactors = TRUE, header = TRUE)
details$filename <- paste0(details$Sample.Participant, "_", details$Visit)
rownames(details) <- details$filename

details$Group <- factor(details$Enroll.Case.Control,
                        levels = c("Control", "Case"),
                        labels = c("Control", "Parkinson"))

details$HY.Stage <- factor(details$HY.Stage)

details$Condition <- factor(paste0(details$Group, ".", details$EVENT_ID))

#Subset the original matrix column
gene_mtx_up <- gene_mtx[, details$filename]
bsj_mtx_up <- bsj_mtx[, details$filename]

gene_DGE <- DGEList(counts = gene_mtx_up, group = details$Condition)
circ_DGE <- DGEList(counts = bsj_mtx_up,
                    group = details$Condition,
                    lib.size = gene_DGE$samples[, "lib.size"],
                    norm.factors = gene_DGE$samples[, "norm.factors"])

#Checking library sizes
head(circ_DGE$samples)

#Filtering low reads
#circ Filtering
keep <- filterByExpr(circ_DGE)
circ_DGE <- circ_DGE[keep, ]
dim(circ_DGE)

#AverageCPM 
AveLogCPM_circ <- aveLogCPM(circ_DGE)
hist(AveLogCPM_circ)

#Library Size Normalization
circ_DGE_norm <- calcNormFactors(circ_DGE, method = "TMM")

#Cluster the samples together 
logCPM <- cpm(circ_DGE_norm, log=TRUE)
#write.csv(logCPM, file = "~/Desktop/CircRNA-Project/Post_review/EdgeR/results/normalized_counts/PDBP_circ_normalized.csv", row.names = TRUE)

plotMDS(logCPM, col = as.numeric(details$Lab.Plate))

#Creating Design Matrix
design_factor <- as.data.frame(details)

PATNO <- factor(design_factor$Sample.Participant)

EVENT_ID <- factor(design_factor$EVENT_ID)

Group <- factor(design_factor$Group)

design <- model.matrix(~PATNO)

HC.V06 <- Group=="Control" & EVENT_ID=="V06" 
PD.V06 <- Group=="Parkinson" & EVENT_ID=="V06"

design<-cbind(design,HC.V06,PD.V06)
colnames(design)
is.fullrank(design)

#Estimating Dispersion
circ_DGE <- estimateDisp(circ_DGE_norm, design, robust= TRUE)
plotBCV(circ_DGE)
circ_DGE$common.dispersion

#Estimating the QL dispersion 
fit_circ <- glmQLFit(circ_DGE, design, robust=TRUE)
head(fit_circ$coefficients)

#Differential Expression Analysis
V06vsBL.P_res_circ <- glmQLFTest(fit_circ, coef = "PD.V06") #Year 2 vs Baseline
V06vsBL.C_res_circ <- glmQLFTest(fit_circ, coef = "HC.V06") #Year 2 vs Baseline


V06vsBL.P_is.circ <- decideTestsDGE(V06vsBL.P_res_circ, adjust.method = "BH", p.value = 0.05)
V06vsBL.C_is.circ <- decideTestsDGE(V06vsBL.C_res_circ, adjust.method = "BH", p.value = 0.05)

summary(V06vsBL.P_is.circ)
summary(V06vsBL.C_is.circ)

#MA plots PD groups
plotMD(V06vsBL.P_res_circ, status=V06vsBL.P_is.circ, main = "DE CircRNAs: Parkinson's Participants Samples at Year 2 vs Baseline (Paired)")+ 
  abline(h=c(-0.2,0.2), col="red")

#MA plots Control group
plotMD(V06vsBL.C_res_circ, status=V06vsBL.C_is.circ, main = "DE CircRNAs: Control Participants Samples at Year 2 vs Baseline (Paired)")+ 
  abline(h=c(-0.2,0.2), col="red")


#Check the top DEGs in PD participants
#circ
setwd("results/PDBP_circ/")

topc_V06 <- as.data.frame(topTags(V06vsBL.P_res_circ, n = Inf))
#write.csv(topc_V06, file = "DE_Analysis_P_PD_V06vsBL_circ_final.csv")

#Check the top DEGs in PD participants
#circ
topc_V06_HC <- as.data.frame(topTags(V06vsBL.C_res_circ, n = Inf))
#write.csv(topc_V06_HC, file = "DE_Analysis_C_V06vsBL_circ_final.csv")
