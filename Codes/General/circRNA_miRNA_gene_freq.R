#Getting the miRNA-gene_target matrix
target_genes <- read.csv("results/target_genes_mtx.csv", header =TRUE, row.names = 1) #import miRNA-genes frequency table
dim(target_genes)

target_miRNAs <- read.table("results/circRNA_miRNA_pairs.tsv", header =TRUE, sep ="\t") #import circRNA-miRNA frequency table
dim(target_miRNAs)

#Getting the circRNA annotations
circ_info <- read.csv("metadata/circRNAs.csv", header = TRUE)
circ_info <- circ_info[order(circ_info$ID), ]

#Get the rownames from CIRCBaseID 
target_miRNAs <- target_miRNAs %>%
  rownames_to_column("ID_1")

target_miRNAs$ID <- circ_info$CircBaseID
target_miRNAs <- target_miRNAs[, -1]
colnames(target_miRNAs)
rownames(target_miRNAs) <- target_miRNAs$ID

colnames(target_genes)
target_genes <- target_genes %>%
  rownames_to_column("ID")
colnames(target_genes)

#from target_miRNA pairs only take the 101 miRNAs of interest 
target_miRNAs_final <- target_miRNAs[, colnames(target_miRNAs)%in%colnames(target_genes)]

#List both frequency 
combine <- bind_rows(target_miRNAs_final, target_genes) %>%
  group_by(ID) %>% 
  summarise(across(everything(), sum, na.rm = TRUE))
dim(combine)

#write.csv(combine, file = "results/miRNA_targets_circRNA_genes.csv", row.names = FALSE)
