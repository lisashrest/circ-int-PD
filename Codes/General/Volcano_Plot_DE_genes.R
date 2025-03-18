library(ggrepel)
library(tidyverse)

setwd("results/PPMI_gene/")
circ1 <- read.csv("DE_Analysis_P_PD_V06vsBL_gene_final.csv", header = TRUE)

circ1 <- circ1 %>% 
  mutate(
    Expression = case_when(logFC >= log(1.285) & FDR <= 0.05 ~ "Up-regulated",
                           logFC <= -log(1.285) & FDR <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )

circ1$genelabels <- ""
circ1$genelabels <- ifelse(circ1$Gene == "MAN1A2" 
                           | circ1$Gene == "SLC8A1"
                           | circ1$Gene == "ZNF91"
                           | circ1$Gene == "VRK1"
                           | circ1$Gene == "EMB"
                           | circ1$Gene == "PICALM", TRUE, FALSE)

PPMI_V06 <- ggplot(data=circ1, aes(x=logFC, y=-log10(FDR), 
                                   col=Expression, label=ifelse(circ1$genelabels,
                                                                circ1$Gene, ""))) +
  geom_point() + 
  theme_minimal() +
  geom_label_repel(alpha = 1, size = 4, 
                   fontface = 'bold',
                   box.padding = 0.5, point.padding = 0.5,
                   na.rm=TRUE, max.overlaps = Inf)+
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.25, 0.25), col="red", linetype = "dotdash") +
  geom_hline(yintercept=-log10(0.05), col="red", linetype = "dotdash")+
  ylim(0,5)+
  xlim(-1,1)+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text = element_text(size =12))+
  ggtitle("PPMI Cohort (Year 2 vs Baseline)")+
  labs(caption = "Total variables = 51,701")+
  theme(plot.caption = element_text(size = 15),
        axis.text = element_text(size = 16),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 18),
        legend.text = element_text(size =16), 
        legend.title = element_text(size = 18))

PPMI_V06

circ2 <- read.csv("DE_Analysis_P_PD_V08vsBL_gene_final.csv", header = TRUE)
circ2 <- circ2 %>% 
  mutate(
    Expression = case_when(logFC >= log(1.285) & FDR <= 0.05 ~ "Up-regulated",
                           logFC <= -log(1.285) & FDR <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )

circ2$genelabels <- ""
circ2$genelabels <- ifelse(circ2$Gene == "MAN1A2" 
                           | circ2$Gene == "SLC8A1"
                           | circ2$Gene == "ZNF91"
                           | circ2$Gene == "VRK1"
                           | circ2$Gene == "EMB"
                           | circ2$Gene == "PICALM", TRUE, FALSE)

PPMI_V08 <- ggplot(data=circ2, aes(x=logFC, y=-log10(FDR), 
                                   col=Expression, label=ifelse(circ2$genelabels,
                                                                circ2$Gene, ""))) +
  geom_point() + 
  theme_minimal() +
  geom_label_repel(alpha = 1, size = 4, 
                   fontface = 'bold',
                   box.padding = 0.5, point.padding = 0.5,
                   na.rm=TRUE, max.overlaps = Inf)+
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.25, 0.25), col="red", linetype = "dotdash") +
  geom_hline(yintercept=-log10(0.05), col="red", linetype = "dotdash")+
  ylim(0,5)+
  xlim(-1,1)+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text = element_text(size =12))+
  ggtitle("PPMI Cohort (Year 3 vs Baseline)")+
  labs(caption = "Total variables = 51,701")+
  theme(plot.caption = element_text(size = 15),
        axis.text = element_text(size = 16),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 18),
        legend.text = element_text(size =16), 
        legend.title = element_text(size = 18))

PPMI_V08


circ3 <- read.csv("DE_Analysis_P_V08vsV06_gene_final.csv", header = TRUE)
circ3 <- circ3 %>% 
  mutate(
    Expression = case_when(logFC >= log(1.285) & FDR <= 0.05 ~ "Up-regulated",
                           logFC <= -log(1.285) & FDR <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )

circ3$genelabels <- ""
circ3$genelabels <- ifelse(circ3$Gene == "MAN1A2" 
                           | circ3$Gene == "SLC8A1"
                           | circ3$Gene == "ZNF91"
                           | circ3$Gene == "VRK1"
                           | circ3$Gene == "EMB"
                           | circ3$Gene == "PICALM", TRUE, FALSE)

PPMI_V08_V06 <- ggplot(data=circ3, aes(x=logFC, y=-log10(FDR), 
                                   col=Expression, label=ifelse(circ3$genelabels,
                                                                circ3$Gene, ""))) +
  geom_point() + 
  theme_minimal() +
  geom_label_repel(alpha = 1, size = 4, 
                   fontface = 'bold',
                   box.padding = 0.5, point.padding = 0.5,
                   na.rm=TRUE, max.overlaps = Inf)+
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.25, 0.25), col="red", linetype = "dotdash") +
  geom_hline(yintercept=-log10(0.05), col="red", linetype = "dotdash")+
  ylim(0,5)+
  xlim(-1,1)+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text = element_text(size =12))+
  ggtitle("PPMI Cohort (Year 3 vs Year 2)")+
  labs(caption = "Total variables = 51,701")+
  theme(plot.caption = element_text(size = 15),
        axis.text = element_text(size = 16),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 18),
        legend.text = element_text(size =16), 
        legend.title = element_text(size = 18))

PPMI_V08_V06


#Plotting Volcano plots for control samples
circ4 <- read.csv("DE_Analysis_C_V06vsBL_gene_final.csv", header = TRUE)

circ4 <- circ4 %>% 
  mutate(
    Expression = case_when(logFC >= log(1.285) & FDR <= 0.05 ~ "Up-regulated",
                           logFC <= -log(1.285) & FDR <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )

PPMI_V06_HC <- ggplot(data=circ4, aes(x=logFC, y=-log10(FDR), col=Expression)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("black", "red")) +
  geom_vline(xintercept=c(-0.2, 0.2), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")+
  ylim(0,5)+
  xlim(-2,2)+
  ggtitle("PPMI Cohort (Year 2 vs Baseline)")+
  labs(caption = "Total variables = 51,701")+
  theme(plot.caption = element_text(size = 15),
        axis.text = element_text(size = 16),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 18),
        legend.text = element_text(size =16), 
        legend.title = element_text(size = 18))

PPMI_V06_HC


circ5 <- read.csv("DE_Analysis_C_V08vsBL_gene_final.csv", header = TRUE)
circ5 <- circ5 %>% 
  mutate(
    Expression = case_when(logFC >= log(1.285) & FDR <= 0.05 ~ "Up-regulated",
                           logFC <= -log(1.285) & FDR <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )

#PDBP_V06
setwd("results/PDBP_gene/")
circ7 <- read.csv("DE_Analysis_P_PD_V06vsBL_gene_final.csv", header = TRUE)

circ7 <- circ7 %>% 
  mutate(
    Expression = case_when(logFC >= log(1.285) & FDR <= 0.05 ~ "Up-regulated",
                           logFC <= -log(1.285) & FDR <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )

circ7$genelabels <- ""
circ7$genelabels <- ifelse(circ7$Gene == "MAN1A2" 
                           | circ7$Gene == "SLC8A1"
                           | circ7$Gene == "ZNF91"   
                           | circ7$Gene == "VRK1"
                           | circ7$Gene == "EMB"
                           | circ7$Gene == "PICALM", TRUE, FALSE)



PDBP_V06 <- ggplot(data=circ7, aes(x=logFC, y=-log10(FDR), 
                                   col=Expression,
                                   label= ifelse(circ7$genelabels == TRUE, circ7$Gene, ""))) +
  geom_point() + 
  theme_minimal() +
  geom_label_repel(alpha = 1, size = 4, 
                   fontface = 'bold',
                   box.padding = 0.5, point.padding = 0.5,
                   na.rm=TRUE, max.overlaps = Inf)+
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.25, 0.25), col="red", linetype = "dotdash") +
  geom_hline(yintercept=-log10(0.05), col="red", linetype = "dotdash")+
  ylim(0,5)+
  xlim(-1,1)+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text = element_text(size =12))+
  ggtitle("PDBP Cohort (Year 2 vs Baseline)")+
  labs(caption = "Total variables = 42,043")+
  theme(plot.caption = element_text(size = 15),
        axis.text = element_text(size = 16),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 18),
        legend.text = element_text(size =16), 
        legend.title = element_text(size = 18))

PDBP_V06

#Healthy Control
circ8 <- read.csv("DE_Analysis_C_V06vsBL_circ_final.csv", header = TRUE)

circ8 <- circ8 %>% 
  mutate(
    Expression = case_when(logFC >= log(1.285) & FDR <= 0.05 ~ "Up-regulated",
                           logFC <= -log(1.285) & FDR <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )

PDBP_V06_HC <- ggplot(data=circ8, aes(x=logFC, y=-log10(FDR), col=Expression)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("black", "red")) +
  geom_vline(xintercept=c(-0.25, 0.25), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")+
  ylim(0,5)+
  xlim(-2,2)+
  ggtitle("PDBP Cohort (Year 2 vs Baseline)")+
  labs(caption = "Total variables = 42043")+
  theme(plot.caption = element_text(size = 15),
        axis.text = element_text(size = 16),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 18),
        legend.text = element_text(size =16), 
        legend.title = element_text(size = 18))

cowplot::plot_grid(PPMI_V06, PDBP_V06, 
                   labels = c("a","b"),
                   ncol = 2, nrow = 1, label_size = 20)
