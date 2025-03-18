library(org.Hs.eg.db)
library(DESeq2)
library(tidyverse)

#Set your working directory
setwd("/data/")

#import data from ciriquant
lib_mtx <- read.csv("library_info.csv", row.names = 1)
gene_mtx <- read.csv("gene_count_matrix.csv", row.names = 1, check.names=FALSE)
gene_mtx <- gene_mtx[ , rownames(lib_mtx)]
bsj_mtx <- read.csv("circRNA_bsj.csv", row.names = 1, check.names=FALSE)
details <- read.csv("sample_details_1.csv", stringsAsFactors = FALSE, header = TRUE)
details$filename <- paste0("P1_", details$PATNO,"_", details$EVENT_ID)

details$Sex <- factor(details$gen, 
                      levels = c(1,2),
                      labels = c("Male", "Female"))

details$Visit <- factor(details$EVENT_ID, 
                        levels = c("BL", "V06", "V08"),
                        labels = c("Baseline", "Year 2", "Year 3"))

details$Group <- factor(details$APPRDX, 
                        levels = c(1,2),
                        labels = c("Parkinson", "Control"))

details$Family_history <- factor(details$fampd_new,
                                 levels = c(1,2,3),
                                 labels = c("1st_degree", "2nd_degree", "No History"))


details$race <- factor(details$race,
                       levels = c(1,2,3,4),
                       labels = c("White", "Black", "Asian", "Other"))


details$HY_stage <- factor(details$NHY, 
                           levels = c(0,1,2,3,4,5))

details$SITE <- factor(details$SITE)


#Creating DESeq2 object
bsj_mtx_up <- bsj_mtx[, details$filename]

#Removing chrX and chrY circRNAs from the circRNA count matrix
filtered_counts <-  bsj_mtx_up %>%
  as.data.frame() %>%
  rownames_to_column("ID")

selected_circRNAs <- filtered_counts[!grepl("chrX:", filtered_counts$ID), ]
selected_circRNAs <- selected_circRNAs[!grepl("chrY:", selected_circRNAs$ID), ]
ID <- selected_circRNAs$ID

#final count  after chromosome X/Y derived circRNAs
bsj_mtx_final <- bsj_mtx_up[rownames(bsj_mtx_up)%in%ID,]
dim(bsj_mtx_final)

#Creating a DESeq2 object 
dds <- DESeqDataSetFromMatrix(countData = bsj_mtx_final,
                              colData = details,
                              design = ~ PATNO + EVENT_ID)
dds

#Pre-filtering
keep <- rowSums(counts(dds) >= 2) >= 140
dds <- dds[keep,]
dds

#Normalizing counts using VST
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
dim(vsd)

S1_a <- plotPCA(vsd, intgroup = "age", ntop = 500)+
  labs(colour = "Age")+
  theme(legend.title = element_text(size = 18),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18), 
        legend.text = element_text(size = 15))


S1_b <- plotPCA(vsd, intgroup = "Sex", ntop = 500)+
  labs(colour = "Sex")+
  theme(legend.title = element_text(size = 18),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 15))

S1_c <- plotPCA(vsd, intgroup = "SITE", ntop = 500)+
  labs(colour = "SITE")+
  theme(legend.title = element_text(size = 18),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 15))

S1_d <- plotPCA(vsd, intgroup = "Family_history", ntop = 500)+
  labs(colour = "Family History")+
  theme(legend.title = element_text(size = 18),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 15))

S1_e <- plotPCA(vsd, intgroup = "race", ntop = 500)+
  labs(colour = "Race")+
  theme(legend.title = element_text(size = 18),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 15))

S1_f <- plotPCA(vsd, intgroup = "HY_stage", ntop = 500)+
  labs(colour = "HY Stage")+
  theme(legend.title = element_text(size = 18),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 15))

#Plot Grid plots
cowplot::plot_grid(S1_a, S1_b,
                   labels = c("a","b"),
                   ncol = 2, nrow = 1, label_size = 20)

cowplot::plot_grid(S1_c, S1_d, 
                   labels = c("c","d"),
                   ncol = 2, nrow = 1, label_size = 20)

cowplot::plot_grid(S1_e, S1_f, 
                   labels = c("e","f"),
                   ncol = 2, nrow = 1, label_size = 20)



#For gene count matrix
gene_mtx_up <- gene_mtx[, details$filename]

#Removing chrX and chrY circRNAs from the filtered deseq2 object
filtered_counts <- gene_mtx_up %>%
  as.data.frame() %>%
  rownames_to_column("ID")%>%
  separate(col = "ID", into = c("ENSEMBL", "Gene"), sep = "\\|")%>%
  separate(col = "ENSEMBL", into = c("ensembl_gene_id_version", "version"), sep = "_")

#Get the chromosome info
library(biomaRt)
grch37 <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org",
                  path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=grch37)

t2g <-getBM(attributes=c('ensembl_gene_id',"ensembl_gene_id_version",
                         "external_gene_name",'chromosome_name','start_position',
                         'end_position'), mart = ensembl)

gene1 <- filtered_counts$ensembl_gene_id_version #vector of your gene names
gene1 <- as.data.frame(gene1)
colnames(gene1) <- "ensembl_gene_id_version"

my_ids.version <- merge(gene1, t2g, by= 'ensembl_gene_id_version')
dim(my_ids.version)
head(my_ids.version)

#Get the columns of interest
id_info <- my_ids.version[, c("ensembl_gene_id_version", "chromosome_name")]
id_info$Chromosome <- paste0("chr", id_info$chromosome_name)
id_info <- id_info %>%
  dplyr::select(c("ensembl_gene_id_version", "Chromosome"))

#merge chromosome info 
filtered_counts_up <- left_join(filtered_counts, id_info, by = "ensembl_gene_id_version")

selected_genes <- filtered_counts_up[!grepl("chrX", filtered_counts_up$Chromosome), ]
selected_genes <- selected_genes[!grepl("chrY", selected_genes$Chromosome), ]

selected_genes$ID <- paste0(selected_genes$ensembl_gene_id_version, "_",
                               selected_genes$version, "|", selected_genes$Gene)

ID <- selected_genes$ID

#final count  after chromosome X/Y derived circRNAs
gene_mtx_final <- gene_mtx_up[rownames(gene_mtx_up)%in%ID,]
dim(gene_mtx_final)

#Creating DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = gene_mtx_final,
                              colData = details,
                              design = ~ PATNO + EVENT_ID)
dds

#Pre-filtering genes with with 10 read counts in at least 30 % of the samples
keep <- rowSums(counts(dds) >= 10) >= 140
dds <- dds[keep,]
dds

#Transformation of gene counts using VST
vsd <- varianceStabilizingTransformation(dds)
dim(vsd)
set.seed(1)

S2_g <- plotPCA(vsd, intgroup = "age", ntop = 500)+
  labs(colour = "Age")+
  theme(legend.title = element_text(size = 18),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18), 
        legend.text = element_text(size = 15))

S2_h <- plotPCA(vsd, intgroup = "Sex", ntop = 500)+
  labs(colour = "Sex")+
  theme(legend.title = element_text(size = 18),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 15))

S2_i <- plotPCA(vsd, intgroup = "SITE", ntop = 500)+
  labs(colour = "SITE")+
  theme(legend.title = element_text(size = 18),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 15))

S2_j <- plotPCA(vsd, intgroup = "Family_history", ntop = 500)+
  labs(colour = "Family History")+
  theme(legend.title = element_text(size = 18),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 15))

S2_k <- plotPCA(vsd, intgroup = "race", ntop = 500)+
  labs(colour = "Race")+
  theme(legend.title = element_text(size = 18),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 15))

S2_l <- plotPCA(vsd, intgroup = "HY_stage", ntop = 500)+
  labs(colour = "HY Stage")+
  theme(legend.title = element_text(size = 18),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 15))

#Plot Grid plots
cowplot::plot_grid(S2_g, S2_h, S2_i, S2_j, S2_k, S2_l,
                   labels = c("a","b","c","d","e","f"),
                   ncol = 3, nrow = 2, label_size = 20)


cowplot::plot_grid(S2_g, S2_h, 
                   labels = c("a","b"),
                   ncol = 2, nrow = 1, label_size = 20)

cowplot::plot_grid(S2_i, S2_j, 
                   labels = c("c","d"),
                   ncol = 2, nrow = 1, label_size = 20)

cowplot::plot_grid(S2_k, S2_l, 
                   labels = c("e","f"),
                   ncol = 2, nrow = 1, label_size = 20)

#Checking lower PCs
S2_g <- plotPCA.san(vsd, intgroup = "age", ntop = 500)+
  labs(colour = "Age")+
  theme(legend.title = element_text(size = 18),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18), 
        legend.text = element_text(size = 15))

S2_h <- plotPCA.san(vsd, intgroup = "Sex", ntop = 500)+
  labs(colour = "Sex")+
  theme(legend.title = element_text(size = 18),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 15))

S2_i <- plotPCA.san(vsd, intgroup = "SITE", ntop = 500)+
  labs(colour = "SITE")+
  theme(legend.title = element_text(size = 18),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 15))

S2_j <- plotPCA.san(vsd, intgroup = "Family_history", ntop = 500)+
  labs(colour = "Family History")+
  theme(legend.title = element_text(size = 18),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 15))

S2_k <- plotPCA.san(vsd, intgroup = "race", ntop = 500)+
  labs(colour = "Race")+
  theme(legend.title = element_text(size = 18),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 15))

S2_l <- plotPCA.san(vsd, intgroup = "HY_stage", ntop = 500)+
  labs(colour = "HY Stage")+
  theme(legend.title = element_text(size = 18),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 15))

#Plot Grid plots
cowplot::plot_grid(S2_g, S2_h, S2_i, S2_j, S2_k, S2_l,
                   labels = c("a","b","c","d","e","f"),
                   ncol = 3, nrow = 2, label_size = 20)
