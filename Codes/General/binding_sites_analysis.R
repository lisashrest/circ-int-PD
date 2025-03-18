library(tidyverse)
setwd("results/miRanda/")

raw_bindSites <- read.table("MirandaOutput.tsv", header = T, sep = '\t', stringsAsFactors = F)
bindSites <- raw_bindSites[,c(1,2)]

allBindSites <- dplyr::count(bindSites, Seq1, name="freq")
distinct <- distinct(bindSites)
distinctBindSites <- dplyr::count(distinct, Seq1, name="freq")

#write.table(allBindSites, file = paste0("bindsites_per_circRNA.tsv"), sep = "\t", quote = F, row.names = F)

# miRanda score distribution
scores <- raw_bindSites[,c(1,3)]

# fiter 25% worst scores
filtered_scores <- raw_bindSites[raw_bindSites$Tot.Score > quantile(raw_bindSites$Tot.Score, 0.25),]
#write.table(filtered_scores, file = paste0("bindsites_25%_filtered.tsv"), sep = "\t", quote = F, row.names = F)

#create an id for miRNA-circRNA pairs
filtered_scores$id <- paste0(filtered_scores$Seq2,".", filtered_scores$Seq1)

#import TarPmiR 
res_TarPmiR <- read.csv("results/TarPmiR.csv")
res_TarPmiR$id <- paste0(res_TarPmiR$Seq2, ".", res_TarPmiR$miRNA)

#Only select the observations in filtered scores also selected by TarPmiR
filtered_scores <- filtered_scores[filtered_scores$id%in%res_TarPmiR$id, ]

#add circRNA information to the filtered score table 
circ_info <- read.csv("results/miRanda/circRNAs.csv", header = TRUE)
head(circ_info)

#Create an id for circRNA-miRNA pairs
filtered_scores <- left_join(filtered_scores, circ_info, by = "Seq2")
#write.csv(filtered_scores, "miRanda/filtered_circRNA_miRNAs.csv")

pairs_freq <- table(filtered_scores$ID, filtered_scores$Seq1)
pairs_mtx <- as.matrix(pairs_freq)
#write.table(pairs_mtx, file = paste0("circRNA_miRNA_pairs.tsv"), sep = "\t", quote = F, row.names = T)
