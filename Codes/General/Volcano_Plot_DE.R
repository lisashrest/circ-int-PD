library(ggrepel)
setwd("results/PPMI_circ/")
circ1 <- read.csv("DE_Analysis_P_PD_V06vsBL_circ_final.csv", header = TRUE)

circ1 <- circ1 %>% 
  mutate(
    Expression = case_when(logFC >= log(1.285) & FDR <= 0.05 ~ "Up-regulated",
                           logFC <= -log(1.285) & FDR <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )

circ1$genelabels <- ""
circ1$genelabels <- ifelse(circ1$Symbol == "MAN1A2" 
                          | circ1$Symbol == "SLC8A1"
                          | circ1$Symbol == "ZNF91"
                          | circ1$Symbol == "VRK1"
                          | circ1$Symbol == "EMB"
                          | circ1$Symbol == "PICALM", TRUE, FALSE)

PPMI_V06 <- ggplot(data=circ1, aes(x=logFC, y=-log10(FDR), 
                                  col=Expression, label=ifelse(circ1$genelabels,
                                                               paste0("circ",circ1$Symbol), ""))) +
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
  labs(caption = "Total variables = 1105")+
  theme(plot.caption = element_text(size = 15),
        axis.text = element_text(size = 16),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 18),
        legend.text = element_text(size =16), 
        legend.title = element_text(size = 18))

PPMI_V06

#V08 vs BL
circ2 <- read.csv("DE_Analysis_P_PD_V08vsBL_circ_final.csv", header = TRUE)


circ2 <- circ2 %>% 
  mutate(
    Expression = case_when(logFC >= log(1.285) & FDR <= 0.05 ~ "Up-regulated",
                           logFC <= -log(1.285) & FDR <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )

circ2$genelabels <- ""
circ2$genelabels <- ifelse(circ2$Symbol == "MAN1A2" 
                          | circ2$Symbol == "SLC8A1"
                          | circ2$Symbol == "ZNF91"
                          | circ2$Symbol == "VRK1"
                          | circ2$Symbol == "EMB"
                          | circ2$Symbol == "PICALM", TRUE, FALSE)


# plot adding up all layers we have seen so far
PPMI_V08 <- ggplot(data=circ2, aes(x=logFC, y=-log10(FDR), 
                                   col=Expression, 
                                   label=ifelse(circ2$genelabels, paste0("circ",circ2$Symbol), ""))) +
  geom_point() + 
  theme_minimal() +
  geom_label_repel(alpha = 1, size = 4, 
                   fontface = 'bold',
                   box.padding = 0.5, point.padding = 0.3,
                   na.rm=TRUE, max.overlaps = 30, seed=45)+
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.25, 0.25), col="red", linetype = "dotdash") +
  geom_hline(yintercept=-log10(0.05), col="red", linetype = "dotdash")+
  ylim(0,5)+
  xlim(-1,1)+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text = element_text(size =12))+
  ggtitle("PPMI Cohort (Year 3 vs Baseline)")+
  labs(caption = "Total variables = 1105")+
  theme(plot.caption = element_text(size = 15),
        axis.text = element_text(size = 16),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 18),
        legend.text = element_text(size =16), 
        legend.title = element_text(size = 18))


PPMI_V08

#Plotting Volcano plots for control samples
circ3 <- read.csv("DE_Analysis_C_V06vsBL_circ_final.csv", header = TRUE)

circ3 <- circ3 %>% 
  mutate(
    Expression = case_when(logFC >= log(1.285) & FDR <= 0.05 ~ "Up-regulated",
                           logFC <= -log(1.285) & FDR <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )

PPMI_V06_HC <- ggplot(data=circ3, aes(x=logFC, y=-log10(FDR), col=Expression)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("black", "red")) +
  geom_vline(xintercept=c(-0.2, 0.2), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")+
  ylim(0,5)+
  xlim(-2,2)+
  ggtitle("PPMI Cohort (Year 2 vs Baseline)")+
  labs(caption = "Total variables = 1105")+
  theme(plot.caption = element_text(size = 15),
        axis.text = element_text(size = 16),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 18),
        legend.text = element_text(size =16), 
        legend.title = element_text(size = 18))

PPMI_V06_HC

circ4 <- read.csv("DE_Analysis_C_V08vsBL_circ_final.csv", header = TRUE)
circ4 <- circ4 %>% 
  mutate(
    Expression = case_when(logFC >= log(1.285) & FDR <= 0.05 ~ "Up-regulated",
                           logFC <= -log(1.285) & FDR <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )

# plot adding up all layers we have seen so far
PPMI_V08_HC <- ggplot(data=circ4, aes(x=logFC, y=-log10(FDR), col=Expression)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("black", "red")) +
  geom_vline(xintercept=c(-0.2, 0.2), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")+
  ylim(0,5)+
  xlim(-2,2)+
  ggtitle("PPMI Cohort (Year 3 vs Baseline)")+
  labs(caption = "Total variables = 1105")+
  theme(plot.caption = element_text(size = 15),
        axis.text = element_text(size = 16),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 18),
        legend.text = element_text(size =16), 
        legend.title = element_text(size = 18))

PPMI_V08_HC

#PDBP_V06
setwd("results/PDBP_circ/")
circ7 <- read.csv("DE_Analysis_P_PD_V06vsBL_circ_final.csv", header = TRUE)

circ7 <- circ7 %>% 
  mutate(
    Expression = case_when(logFC >= log(1.285) & FDR <= 0.05 ~ "Up-regulated",
                           logFC <= -log(1.285) & FDR <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )

circ7$genelabels <- ""
circ7$genelabels <- ifelse(circ7$Symbol == "MAN1A2" 
                           | circ7$Symbol == "ZNF124"
                           | circ7$Symbol == "SLC8A1"
                           | circ7$Symbol == "ZNF91"   
                           | circ7$Symbol == "VRK1"
                           | circ7$Symbol == "EMB"
                           | circ7$Symbol == "PICALM", TRUE, FALSE)



PDBP_V06 <- ggplot(data=circ7, aes(x=logFC, y=-log10(FDR), 
                                  col=Expression, label=ifelse(circ7$genelabels,
                                                               paste0("circ",circ7$Symbol), ""))) +
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
  labs(caption = "Total variables = 774")+
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
  labs(caption = "Total variables = 774")+
  theme(plot.caption = element_text(size = 15),
        axis.text = element_text(size = 16),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 18),
        legend.text = element_text(size =16), 
        legend.title = element_text(size = 18))

cowplot::plot_grid(PPMI_V06, PDBP_V06, 
                   labels = c("a","b"),
                   ncol = 2, nrow = 1, label_size = 20)

