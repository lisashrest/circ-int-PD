library(tidyverse)
library(clusterProfiler)
library(ggVennDiagram)
library(readxl)
library(AnnotationDbi)
library(org.Hs.eg.db)

#Importing the files
#Loading DEGs from PPMI
output_genes_PPMI <- read.csv("results/PPMI/DE_Analysis_P_PD_V06vsBL_gene_final.csv")

#Significant Genes
#output_PPMI <- separate(data = output_genes_PPMI, col = Ensembl, into = c("Ensembl", "version"), sep = "\\.")

#Down-regulated genes
#Selecting genes with FDR < 0.05 
de_genes_PPMI <- output_genes_PPMI %>%
  filter(FDR < 0.05 & logFC < -0.25)

#Target genes
#Creating a df with miRNA targets
#Import the files
#setwd("circRNA-sponging_pipeline/")

#Collapse all the tables into one
miRNA_data <- read.csv("results/sponge_output/sig_circ_ceRNA_interactions.csv", header = TRUE)

#Visualize the intersection in a venn diagram
list_PPMI <- list(de_genes_PPMI$Gene, miRNA_data$`Target gene`)

#Import the WGCNA module file 
module_genes <- read.csv("results/WGCNA/geneInfo_Y2.csv", header = TRUE)
head(module_genes)

#Subsset gene modules with significant weak negative correlation to Year 2
turquoise <- module_genes %>%
  dplyr::filter(moduleColor == "turquoise")%>%
  dplyr::filter(GS.Year_2 < -0.15 & p.GS.Year_2 < 0.05 & MM.turquoise > 0.8)

#Subset gene modules with significant weak negative correlation to Year 2
pink <- module_genes %>%
  dplyr::filter(moduleColor == "pink")%>%
  dplyr::filter(GS.Year_2 < -0.15 & p.GS.Year_2 < 0.05 & MM.pink > 0.8)

#Plotting the genes of interest in the turquoise module
list_PPMI <- list(de_genes_PPMI$Gene, miRNA_data$geneB, turquoise$geneSymbol)

#Creative Venn Diagram 
ggVennDiagram(list_PPMI, set_size = 5, label = "count", 
              label_alpha = 0, label_size = 10, label_txtWidth = 10,
              edge_size = 2, 
              category.names = c("Down-reg", "Target genes (miRNAs)", 
                                 "Turquoise"))+
  scale_fill_distiller(palette = "Pastel2", direction = 1)+
  theme(plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),
        legend.position = "none",
        )+
  scale_x_continuous(expand = expansion(mult = .2))

#Pink
list_PPMI_pink <- list(de_genes_PPMI$Gene, miRNA_data$geneB, pink$geneSymbol)

#Creative Venn Diagram 
ggVennDiagram(list_PPMI_pink, set_size = 5, label = "count", label_alpha = 0, label_size = 10, label_txtWidth = 10, edge_size = 2,
              category.names = c("Down-reg", "Target genes (miRNAs)", "Pink"))+
  scale_fill_distiller(palette = "Pastel2", direction = 1)+
  theme(legend.position = "none")+
  theme(plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5))+
  scale_x_continuous(expand = expansion(mult = .2))

