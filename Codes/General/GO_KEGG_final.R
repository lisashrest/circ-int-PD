library(tidyverse)
library(clusterProfiler)
library(ggVennDiagram)
library(readxl)
library(AnnotationDbi)
library(org.Hs.eg.db)

#Importing the files
#Loading the list of genes in the ceRNA network
ceRNA <- read.csv("sponge_output/sig_circ_ceRNA_interactions.csv", header = TRUE, row.names = 1)
ceRNA$Gene <- ceRNA$geneB

#loading differential expression data for the year 2 and baseline comparison in PPMI data
norm_cts <- read.csv("normalized_counts/PPMI_gene_logCPM.csv", header = TRUE, row.names =1)
de_genes <- read.csv("results/DE_Analysis_P_PD_V06vsBL_gene_final.csv", header = TRUE)

de_genes$ENTREZID <- mapIds(org.Hs.eg.db, keys = de_genes$Gene, keytype = "SYMBOL",
                            column = "ENTREZID")

de_genes_dr <- de_genes %>%
  filter(FDR < 0.05)%>%
  filter(logFC < -0.25)

#add the de_genes information to ceRNA Network data
full_ceRNA_dr <- inner_join(ceRNA, de_genes_dr, by = "Gene")

#Import the WGCNA 
WGCNA_res <- read.csv("results/geneInfo_Y2.csv", header = TRUE)

#Selecting only the genes in turquoise and pink modules 
WGCNA_res_sig_pink <- WGCNA_res %>%
  filter(moduleColor == "pink")%>%
  filter(GS.Year_2 < -0.15 & p.GS.Year_2 < 0.05 & MM.pink > 0.8)

WGCNA_res_sig_turquoise <- WGCNA_res %>%
  filter(moduleColor == "turquoise")%>%
  dplyr::filter(GS.Year_2 < -0.15 & p.GS.Year_2 < 0.05 & MM.turquoise > 0.8)

WGCNA_sig_pink <- WGCNA_res_sig_pink %>%
  mutate(MM = MM.pink)%>%
  select(c(1:7,44))

WGCNA_sig_turquoise <- WGCNA_res_sig_turquoise %>%
  mutate(MM = MM.turquoise)%>%
  select(c(1:7,44))

WGCNA_sig <- rbind(WGCNA_res_sig_pink, WGCNA_res_sig_turquoise)

#Only select ceRNA genes present in WGCNA results
ceRNA_dr <- full_ceRNA_dr[full_ceRNA_dr$Gene%in%WGCNA_sig$geneSymbol, ]
ceRNA_dr <- ceRNA_dr[order(ceRNA_dr$logFC, decreasing = TRUE),]

#list of de genes in ceRNA network
genes_of_interest <- ceRNA_dr$ID
genes_of_interest_mtx <- norm_cts[rownames(norm_cts)%in%genes_of_interest,]
#write.csv(genes_of_interest_mtx, file = "results/genes_of_interest.csv")

#ClusterProfiler
geneList_PPMI <- ceRNA_dr$PValue
names(geneList_PPMI) <- ceRNA_dr$ENTREZID

#import annotation folder
gene_universe <- de_genes$ENTREZID
head(gene_universe)

#Creating input for over-representation analysis
gene_PPMI <- names(geneList_PPMI)
head(gene_PPMI)

#Creating input for enrichment analysis
lfc_list_PPMI <- ceRNA_dr$logFC
names(lfc_list_PPMI) <- ceRNA_dr$ENTREZ_ID

lfc_list_final_PPMI <- sort(lfc_list_PPMI, decreasing =  TRUE)

#GO ORA using clusterProfiler
ego_PPMI <- enrichGO(gene = gene_PPMI,
                     universe = gene_universe,
                     OrgDb = org.Hs.eg.db,
                     ont = "All",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
head(ego_PPMI)
write.csv(ego_PPMI, file = "results/GO_clusterProfiler_downreg.csv", row.names = FALSE)


#Gene set enrichment analysis using clusterProfiler
gsea_all_PPMI <- gseGO(geneList = lfc_list_final_PPMI,
                       OrgDb = org.Hs.eg.db,
                       ont = "All",
                       minGSSize = 10,
                       maxGSSize = 500,
                       pvalueCutoff = 0.05,
                       verbose = FALSE, 
                       eps = 0, 
                       nPermSimple = 10000,
                       scoreType = "pos")

head(gsea_all_PPMI)

#Visualization of the GO Analysis
p1 <- clusterProfiler::dotplot(ego_PPMI, split = "ONTOLOGY", showCategory = 10, font.size = 12) +
  facet_grid(ONTOLOGY~., scale="free") + 
  ggtitle("Downregulated genes")

p1

#KEGG ORA using clusterProfiler
kk_PPMI <- enrichKEGG(gene         = gene_PPMI,
                      organism     = 'hsa',
                      pvalueCutoff = 0.05)
head(kk_PPMI)
write.csv(kk_PPMI, file = "results/KEGG_ORA_downreg.csv")

#KEGG Enrichment Analysis
kk2_PPMI <- gseKEGG(geneList = lfc_list_final_PPMI,
                    organism = 'hsa',
                    minGSSize = 100,
                    maxGSSize = 500,  
                    pvalueCutoff = 1,
                    verbose = FALSE)
head(kk2_PPMI)

#Visualization of the clusterProfiler analysis
p5 <- clusterProfiler::dotplot(kk_PPMI, showCategory = 10, font.size = 12) + 
  ggtitle("Downregulated genes")

p5

cowplot::plot_grid(p1, p5, 
                   labels = c("a", "b"),
                   label_size = 15)
