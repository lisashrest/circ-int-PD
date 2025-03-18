#Spearman Rank Correlarion on the expression of circRNAs in two cohorts: PPMI and PDBP
#Aim: To validate the observations of circRNAs (logFC) on both cohorts is correlated. Validity of the observation in the independent cohort
library(tidyverse)
library(smplot2)

PPMI_circ <- read.csv("results/PPMI_circ/DE_Analysis_P_PD_V06vsBL_circ_final.csv", header = TRUE)

PDBP_circ <- read.csv("results/PDBP_circ/DE_Analysis_P_PD_V06vsBL_circ_final.csv", header = TRUE)

#Get a final list of circRNAs identified in both cohorts
final_circ <- inner_join(PPMI_circ, PDBP_circ, by = "X")  

#Spearmans Rank Correlation 
p <- cor.test(final_circ$logFC.x, final_circ$logFC.y, method = "spearman")


logFC_concordance <- ggplot(data = final_circ, aes(x = logFC.x, y = logFC.y)) +
  geom_point(shape = 21, fill = "black", size = 3) +
  geom_smooth(method = loess, formula = y~x, fill = "grey") +
  ylab("PDBP cohort - logFC (Paired Analysis in PD samples: Year 2 vs Baseline)")+
  xlab("PPMI cohort - logFC (Paired Analysis in PD samples: Year 2 vs Baseline)")+
  annotate("text", x=-0.25, y=1.65, label=paste0("r = ", round(p$estimate,2)), hjust=0, size = 6) +
  annotate("text", x=-0.25, y=1.75, label=paste0("p = < 2.2e-16"), hjust=0, size = 6) +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size =18),
        axis.title.y = element_text(size = 18))

cowplot::plot_grid(logFC_concordance, PPMI_V08, 
                   labels = c("d","e"),
                   ncol = 2, nrow = 1, label_size = 20)
