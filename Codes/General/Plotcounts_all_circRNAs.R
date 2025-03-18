library(tidyverse)
library(dplyr)
library(ggpubr)
library(RColorBrewer)
library(reshape2)

#PPMI Cohort
load("env/PPMI_circRNA.RData")

#Create a vector with position coordinates of the down-regulated circRNAs of interest
circ_oi_down <- c("chr2:40655613|40657444", "chr1:117944808|117963271", "chr19:23541232|23545527",  "chr14:97299804|97327072", "chr5:49694941|49707217", "chr11:85718585|85742653")

#Create a vector with names of the down-regulated circRNAs of interest
circ_name_down <- c("circSLC8A1", "circMAN1A2", "circZNF91", "circVRK1", "circEMB", "circPICALM")

#create a dataframe with ids and repective names of down-regulated circRNAs of interest
circ_data_down <- data.frame(circ_oi_down, circ_name_down)
colnames(circ_data_down) <- c("circRNA", "CircName")

# Extract out the logCPM expression of down-reg circRNAs from the `dge` object:
circ_exp <- circ_DGE %>% 
  cpm(log = TRUE) %>%
  as.data.frame %>%
  rownames_to_column("CircID") %>%
  dplyr::filter(CircID %in% circ_oi_down) %>%
  melt()
colnames(circ_exp) <- c("circRNA", "sample", "logCPM")

samples <- circ_DGE$samples %>%
  as.data.frame() %>%
  rownames_to_column("sample")

meta <- data.frame(details_up$APPRDX, details_up$PD_MED_USE, details_up$filename)
colnames(meta) <- c("Condition", "PD_medication", "sample")
meta

# Add sample metadata information so that we can 
# use the `group` column in dge$samples for plotting. 
circ_exp <- left_join(circ_exp, samples, by = "sample")
circ_exp <- left_join(circ_exp, meta, by = "sample")
circ_exp <- left_join(circ_exp, circ_data_down, by = "circRNA")

circ_exp_PD <- subset(circ_exp, Condition == "1")

# Plot boxplots of the logCPM expression:
my_comparisons <- list(c("Parkinson.BL", "Parkinson.V06"), c("Parkinson.BL", "Parkinson.V08"), c("Parkinson.V06", "Parkinson.V08"))

plotcounts_PPMI_circ <- circ_exp_PD %>%
  ggplot(aes(x = group, y = logCPM, colour = group)) +
  geom_boxplot() +
  geom_jitter() +
  theme(aspect.ratio = 1) +
  ggtitle("PPMI Cohort")+
  theme(axis.text.x = element_text(size = 16, angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 18),
        axis.text.y = element_text(size = 16),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 18),
        legend.text = element_text(size =16), 
        legend.title = element_text(size = 18))+
  stat_compare_means(comparisons = my_comparisons, paired = TRUE, method = "wilcox.test", label = "p.signif")+
  facet_wrap(~ CircName)+
  xlab("Time Point in the Study")+
  scale_x_discrete(labels = c("Baseline", "Year 2", "Year 3"))+
  scale_color_brewer(palette="Set2",name = "Visit", labels = c("Baseline", "Year 2", "Year 3"))+
  ylim(-6,4)

plotcounts_PPMI_circ

#PDBP Cohort
load("env/PDBP_circRNA.RData")

#Create a vector with position coordinates of the down-regulated circRNAs of interest
circ_oi_down <- c("chr2:40655613|40657444", "chr1:117944808|117963271", "chr19:23541232|23545527",  "chr14:97299804|97327072", "chr5:49694941|49707217", "chr11:85718585|85742653")

#Create a vector with names of the down-regulated circRNAs of interest
circ_name_down <- c("circSLC8A1", "circMAN1A2", "circZNF91", "circVRK1", "circEMB", "circPICALM")

#create a dataframe with ids and repective names of down-regulated circRNAs of interest
circ_data_down <- data.frame(circ_oi_down, circ_name_down)
colnames(circ_data_down) <- c("circRNA", "CircName")

# Extract out the logCPM expression of down-reg circRNAs of interest from the `dge` object:
circ_exp <- circ_DGE %>% 
  cpm(log = TRUE) %>%
  as.data.frame %>%
  rownames_to_column("CircID") %>%
  dplyr::filter(CircID %in% circ_oi_down) %>%
  melt()

colnames(circ_exp) <- c("circRNA", "sample", "logCPM")

samples <- circ_DGE$samples %>%
  as.data.frame() %>%
  rownames_to_column("sample")

meta <- data.frame(details$Enroll.Case.Control, details$filename)
colnames(meta) <- c("Condition", "sample")

# Add sample metadata information so that we can 
# use the `group` column in dge$samples for plotting. 
circ_exp <- left_join(circ_exp, samples, by = "sample")
circ_exp <- left_join(circ_exp, meta, by = "sample")
circ_exp <- left_join(circ_exp, circ_data_down, by = "circRNA")

circ_exp_PD <- subset(circ_exp, Condition == "Case")

# Plot boxplots of the logCPM expression:
my_comparisons <- list(c("Parkinson.BL", "Parkinson.V06"))

plotcounts_PDBP_circ <- circ_exp_PD %>%
  ggplot(aes(x = group, y = logCPM, colour = group)) +
  geom_boxplot() +
  geom_jitter()+
  theme(aspect.ratio = 1) +
  ggtitle("PDBP Cohort")+
  theme(axis.text.x = element_text(size = 16, angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 18),
        axis.text.y = element_text(size = 16),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 18),
        legend.text = element_text(size =16), 
        legend.title = element_text(size = 18))+
  stat_compare_means(comparisons = my_comparisons, paired = TRUE, method = "wilcox.test", label = "p.signif")+
  facet_wrap(~ CircName)+
  xlab("Time Point in the Study")+
  scale_x_discrete(labels = c("Baseline", "Year 2"))+
  scale_color_brewer(palette="Set2",name = "Visit", labels = c("Baseline", "Year 2"))+
  ylim(-6,4)

plotcounts_PDBP_circ

#Figure2_Part2
cowplot::plot_grid(plotcounts_PPMI_circ, plotcounts_PDBP_circ, 
                   labels = c("c","d"),
                   ncol = 2, nrow = 1, label_size = 20)
