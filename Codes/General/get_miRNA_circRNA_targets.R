miRNA_circRNA <- read.table("results/filtered_circRNA_miRNA_correlation_v2.tsv", header = TRUE, sep = "\t")
miRNA_circRNA <- miRNA_circRNA[,1:4]

miRNA_circRNA <- miRNA_circRNA %>%
  filter(miRNA_binding_sites > 0)

miRNA_target_list <- list()
circRNAs <- miRNA_circRNA$circRNA[!duplicated(miRNA_circRNA$circRNA)]

for(i in 1:length(circRNAs)){
  circRNA <- circRNAs[i]
  miRNA_targets <- miRNA_circRNA %>%
    filter(circRNA == circRNA)
  miRNA_target_list[[i]] <- miRNA_targets
  names(miRNA_target_list)[i] <- circRNA
}

saveRDS(miRNA_target_list, file = "results/miRNA_target.rds")
