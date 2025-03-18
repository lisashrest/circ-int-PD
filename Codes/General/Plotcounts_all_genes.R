library(tidyverse)
library(dplyr)
library(ggpubr)
library(RColorBrewer)
library(reshape2)

load("env/PPMI_gene.RData")

#TMM Normalized count matrix from EdgeR
counts <- gene_DGE_norm$counts
head(counts)

host_genes <- c("SLC8A1", "EMB", "MAN1A2", "VRK1", "PICALM", "ZNF91")

#Names of the host genes in the host_genes vector
gene_counts <- counts[grepl(paste(host_genes, collapse = "|"), rownames(counts)), ]
gene_counts <- gene_counts[-c(4,5,7), ]
gene_names <- rownames(gene_counts)

#create a dataframe with ids and repective names of host genes
gene_data <- data.frame(host_genes, gene_names)
colnames(gene_data) <- c("Gene_Name", "ID")

# Extract out the logCPM expression of down-reg circRNAs from the `dge` object:
gene_exp <- gene_DGE_norm$counts %>% 
  cpm(log = TRUE) %>%
  as.data.frame %>%
  rownames_to_column("ID")

gene_exp <- gene_DGE_norm$counts %>% 
  cpm(log = TRUE) %>%
  as.data.frame %>%
  rownames_to_column("ID") %>%
  dplyr::filter(ID %in% gene_data$ID)%>%
  separate(col = "ID", into = c("ID", "Version"), sep = "[_]")%>%
  separate(col = "Version", into = c("Version", "Gene_Name"), sep = "\\|")%>%
  melt()

colnames(gene_exp) <- c("ID", "Version", "Gene_Name", "sample", "Normalized reads")

samples <- gene_DGE_norm$samples %>%
  as.data.frame() %>%
  rownames_to_column("sample")

meta <- data.frame(meta_final$Group, meta_final$filename, meta_final$EVENT_ID)
colnames(meta) <- c("Condition", "sample", "EVENT_ID")

# Add sample metadata information so that we can 
# use the `group` column in dge$samples for plotting. 
gene_exp <- left_join(gene_exp, samples, by = "sample")
gene_exp <- left_join(gene_exp, meta, by = "sample")
gene_exp <- left_join(gene_exp, gene_data, by = "ID")

gene_exp_PD <- subset(gene_exp, Group == "Parkinson")

# Plot boxplots of the logCPM expression:
my_comparisons <- list(c("BL", "V06"), c("BL", "V08"), c("V06", "V08"))
time_point <- c("Baseline", "Year 2", "Year 3")

plotcounts_PDBP <- gene_exp_PD %>%
  ggplot(aes(x = EVENT_ID, y = `Normalized Reads`, colour = EVENT_ID)) +
  geom_boxplot() +
  geom_jitter() +
  theme(aspect.ratio = 1) +
  ggtitle("PPMI Cohort (Samples from PD participants only)")+
  theme(axis.text.x = element_text(size = 16, angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 18),
        legend.text = element_text(size =16), 
        legend.title = element_text(size = 18))+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif")+
  facet_wrap(~ Gene_Name.x)+
  xlab("Time Point in the Study")+
  scale_color_brewer(palette="Set2",name = "Visit", labels = c("Baseline", "Year 2", "Year 3"))+
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  )+
  scale_x_discrete(labels = time_point)

#PDBP Cohort
load("~/Desktop/CircRNA-Project/Post_review/EdgeR/env/PDBP_gene.RData")

#TMM Normalized count matrix from EdgeR
counts <- gene_DGE_norm$counts
head(counts)

host_genes <- c("SLC8A1", "EMB", "MAN1A2", "VRK1", "PICALM", "ZNF91")

#Names of the host genes in the host_genes vector
gene_counts <- counts[grepl(paste(host_genes, collapse = "|"), rownames(counts)), ]
gene_counts <- gene_counts[-c(4,5,7), ]

#create a dataframe with ids and repective names of host genes
gene_data <- data.frame(host_genes, rownames(gene_counts))
colnames(gene_data) <- c("Gene_Name", "ID")

# Extract out the logCPM expression of down-reg circRNAs from the `dge` object:
gene_exp <- gene_DGE_norm$counts %>% 
  cpm(log = TRUE) %>%
  as.data.frame%>%
  rownames_to_column("ID")

gene_exp <- gene_DGE_norm$counts %>% 
  cpm(log = TRUE) %>%
  as.data.frame %>%
  rownames_to_column("ID") %>%
  dplyr::filter(ID %in% gene_data$ID)%>%
  separate(col = "ID", into = c("ID", "Version"), sep = "[_]")%>%
  separate(col = "Version", into = c("Version", "Gene_Name"), sep = "\\|")%>%
  melt()

colnames(gene_exp) <- c("ID", "Version", "Gene_Name", "sample", "Normalized Reads")

samples <- gene_DGE_norm$samples %>%
  as.data.frame() %>%
  rownames_to_column("sample")

meta <- data.frame(details$Enroll.Case.Control, details$filename, details$EVENT_ID)
colnames(meta) <- c("Condition", "sample", "EVENT_ID")

# Add sample metadata information so that we can 
# use the `group` column in dge$samples for plotting. 
gene_exp <- left_join(gene_exp, samples, by = "sample")
gene_exp <- left_join(gene_exp, meta, by = "sample")
gene_exp <- left_join(gene_exp, gene_data, by = "ID")

gene_exp_PD <- subset(gene_exp, Condition == "Case")

# Plot boxplots of the logCPM expression:
my_comparisons <- list(c("BL", "V06"))
time_point <- c("Baseline", "Year 2")

plotcounts_PDBP <- gene_exp_PD %>%
  ggplot(aes(x = EVENT_ID, y = `Normalized Reads`, colour = EVENT_ID)) +
  geom_boxplot() +
  geom_jitter() +
  theme(aspect.ratio = 1) +
  ggtitle("PDBP Cohort (Samples from PD participants only)")+
  theme(axis.text.x = element_text(size = 16, angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 18),
        legend.text = element_text(size =16), 
        legend.title = element_text(size = 18))+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif")+
  facet_wrap(~ Gene_Name.x)+
  xlab("Time Point in the Study")+
  scale_color_brewer(palette="Set2",name = "Visit", labels = c("Baseline", "Year 2", "Year 3"))+
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  )+
  scale_x_discrete(labels = time_point)

cowplot::plot_grid(plotcounts_PPMI, plotcounts_PDBP,
                   labels = c("d", "e"))
