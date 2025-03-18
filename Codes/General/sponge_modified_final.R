#setting the miRNA-gene_target matrix
target_scan_symbols_counts <- read.csv("results/miRNA_targets_circRNA_genes.csv", header =TRUE, row.names = 1)
dim(target_scan_symbols_counts)

#Change the miRNA names to miRBase ID
target_scan_symbol_t <- target_scan_symbols_counts %>% 
  t() %>%
  as.data.frame()%>%
  rownames_to_column("ID")

target_scan_symbol_t$ID <- gsub("\\.", "-", target_scan_symbol_t$ID)
rownames(target_scan_symbol_t) <- target_scan_symbol_t$ID

target_scan_symbols_counts <- as.matrix(t(target_scan_symbol_t[,-1]))
#class(target_scan_symbols_counts) <- "numeric"
head(target_scan_symbols_counts)
str(target_scan_symbols_counts)

# SET MIRNA EXPRESSION
print("reading miRNA expression...")
mi_rna_expr <- read.csv(file = "normalized_counts/PPMI_miRNA_logCPM.csv", header = T, check.names = F, row.names = 1)
head(mi_rna_expr)

mi_rna_expr <- mi_rna_expr %>%
  rownames_to_column("ID")

mi_rna_expr$ID <- gsub("\\.", "-", mi_rna_expr$ID)
mi_rna_expr <- mi_rna_expr %>%
  column_to_rownames("ID")

head(mi_rna_expr)

# SET GENE EXPRESSION
print("reading gene expression...")
gene_expr <- read.csv(file = "normalized_counts/PPMI_gene_logCPM.csv", header = T, check.names = F, row.names = 1)

#Get the row Ids as a column
gene_expr <- gene_expr %>%
  rownames_to_column("ID")%>%
  separate(col = "ID", into = c("ENSEMBL", "GENE"), sep = "\\|")%>%
  separate(col = "ENSEMBL", into = c("ENSEMBL", "ENSEMBL_ID", "Version"), sep = "\\.|_")

#For duplicated gene names, add number at the second one
gene_expr <- gene_expr %>%
  mutate(dupl = if_else(duplicated(GENE), 1, 0)) %>%
  group_by(GENE) %>%
  mutate(dupl = cumsum(dupl),
         GENE = if_else(dupl > 0, paste(GENE, dupl, sep = "_"), GENE)) %>%
  dplyr::select(-dupl)

rownames(gene_expr) <- gene_expr$GENE
head(gene_expr)

# if (argv$normalize & !argv$tpm) {
#   # normalize expressions if not already done
#   print("normalizing gene expression")
#   gene_expr <- normalize.data(gene_expr)
# }

# READ CIRC_RNA EXPRESSION AND COMBINE THEM
print("adding circRNA expression...")
circ_filtered_raw <- read.csv(file = "normalized_counts/PPMI_circ_logCPM.csv", header = T, check.names = F, row.names =1)
circ_filtered_raw <- circ_filtered_raw %>%
  rownames_to_column("ID")

#import circRNA annotation 
circ_anno <- read.csv("normalized_counts/PPMI_circ_annotations.csv", header =TRUE)
circ_anno <- circ_anno[, c(1,6)]
circ_filtered_raw <- inner_join(circ_filtered_raw, circ_anno, by = "ID")

circ_filtered_raw$CirBaseID[circ_filtered_raw$CirBaseID == ""] <- NA
na_count <- sum(is.na(circ_filtered_raw$CirBaseID))
replacements <- paste0("unannotated", seq_len(na_count))

#Replace rownames with circbaseID
circ_filtered_raw$CirBaseID[is.na(circ_filtered_raw$CirBaseID)] <- replacements
rownames(circ_filtered_raw) <- circ_filtered_raw$CirBaseID

print("reading samplesheet...")
meta <- read.csv("metadata/metadata.csv", check.names = F, row.names = 1)

# filter for expressions only
circ_filtered <- circ_filtered_raw[,meta$filename]

# remove NAs
circ_filtered <- circ_filtered[complete.cases(circ_filtered),]
head(circ_filtered)

#select gene expression only 
gene_filtered <- as.matrix(gene_expr[, meta$filename])
rownames(gene_filtered) <- rownames(gene_expr)
head(gene_filtered)

# combine linear and circular expressions
gene_expr_final <- rbind(circ_filtered, gene_filtered)
mi_rna_expr <- as.matrix(mi_rna_expr[, meta$filename])
dim(mi_rna_expr)

# filter for matching samples
print("gene_expr samples:")
dim(gene_expr_final)
print("miRNA expr samples:")
dim(mi_rna_expr)
print("target scan symbols samples:")
dim(target_scan_symbols_counts)

# shared samples
# shared.samples <- intersect(colnames(gene_expr), colnames(mi_rna_expr))
# all.samples <- unique(c(colnames(gene_expr), colnames(mi_rna_expr)))
# if (length(shared.samples) < length(all.samples)) {
#   cat("Discarding samples", all.samples[!all.samples %in% shared.samples], "because they are missing in either gene or miRNA expression\n")
# }
# gene_expr <- gene_expr[,shared.samples]
# mi_rna_expr <- mi_rna_expr[,shared.samples]
# dim(gene_expr)

# transform for sponge
gene_expr_final[is.na(gene_expr_final)] <- 0
mi_rna_expr[is.na(mi_rna_expr)] <- 0

# # use tpms instead of counts
# if (argv$tpm) {
#   # convert circRNA and linear expression
#   print("using TPMs instead of counts for gene and miRNA expression")
#   TPM.map <- read.table(argv$tpm_map, header = T, sep = "\t", check.names = F)
#   gene_expr <- TPM.map[rownames(gene_expr), colnames(gene_expr)]
#   # convert miRNA expression
#   mir_fasta <- readDNAStringSet(argv$mir_fasta)
#   mi_rna_expr <- mi_rna_expr[rownames(mi_rna_expr) %in% names(mir_fasta),]
#   lengths <- mir_fasta[names(mir_fasta)%in%rownames(mi_rna_expr)]@ranges@width
#   mi_rna_expr <- mi_rna_expr/lengths
#   mi_rna_expr <- log2(t(t(mi_rna_expr)*1e6/colSums(mi_rna_expr)) + argv$pseudocount)
  
#   # transform for sponge
#   gene_expr[is.na(gene_expr)] <- 0
#   mi_rna_expr[is.na(mi_rna_expr)] <- 0
# }

# annotate circRNAs if possible
if("CirBaseID" %in% colnames(circ_filtered_raw)) {
  # get all circBase IDs for row names in the circRNA expression file
  IDs <- rownames(gene_expr_final)
  try_an <- circ_filtered_raw[IDs, "circBaseID"]
  new <- which(!is.na(try_an) & try_an != "None")
  # annotate
  IDs[new] <- try_an[new]
  rownames(gene_expr_final) <- IDs
}

# transpose for SPONGE
mi_rna_expr <- as.matrix(t(mi_rna_expr))
gene_expr_final <- as.matrix(t(gene_expr_final))

# cast to matrix
target_scan_symbols_counts <- as.matrix(target_scan_symbols_counts)

print("Gene expression:")
print(gene_expr_final[1:5, 1:5])
print("miRNA expression:")
print(mi_rna_expr[1:5, 1:5])
print("target scan symbols:")
print(target_scan_symbols_counts[1:5, 1:5])

print("calculating gene-miRNA interactions...")
# (A) gene-miRNA interactions
genes_miRNA_candidates <- SPONGE::sponge_gene_miRNA_interaction_filter(
  gene_expr = gene_expr_final,
  mir_expr = mi_rna_expr,
  mir_predicted_targets = target_scan_symbols_counts,
  log.level = "INFO",
  elastic.net = TRUE,
  F.test.p.adj.threshold = 0.05)

save.image(file = file.path(out, "sponge_final.RData"))
print("calculating ceRNA interactions...")
# (B) ceRNA interactions
ceRNA_interactions <- SPONGE::sponge(gene_expr = gene_expr_final,
                                     mir_expr = mi_rna_expr,
                                     mir_interactions = genes_miRNA_candidates)
save.image(file = file.path(out, "sponge_final.RData"))

print("building null model...")
# (C) Null-model-based p-value computation
mscor_null_model <- sponge_build_null_model(number_of_datasets = 100, number_of_samples = nrow(gene_expr_final))

# simulation plot
sim_plot <- sponge_plot_simulation_results(mscor_null_model)
png(file.path(out, "plots/simulation.png"))
plot(sim_plot)

# ceRNA interaction signs
ceRNA_interactions_sign <- sponge_compute_p_values(sponge_result = ceRNA_interactions, 
                                                   null_model = mscor_null_model)
print("building ceRNA network...")

# (D) ceRNA interaction network
fdr <- as.double(0.05)
min.interactions <- 5000
ceRNA_interactions_fdr <- ceRNA_interactions_sign[which(ceRNA_interactions_sign$p.adj < fdr),]
if (nrow(ceRNA_interactions_fdr)<min.interactions && nrow(ceRNA_interactions_sign)>min.interactions) {
  print("Warning: fdr setting too strict, no significant interactions detected; min of padj is:")
  print(min(ceRNA_interactions_sign$p.adj))
  print("adjusting...")
  fdr <- min(ceRNA_interactions_sign$p.adj) * 1.01
  ceRNA_interactions_fdr <- ceRNA_interactions_sign[which(ceRNA_interactions_sign$p.adj < fdr),]
  while (nrow(ceRNA_interactions_fdr)<min.interactions) {
    fdr <- fdr * 1.01
    ceRNA_interactions_fdr <- ceRNA_interactions_sign[which(ceRNA_interactions_sign$p.adj < fdr),]
  }
  cat("adjusted fdr to :", fdr, "to allow for a minimum", min.interactions, "interactions", "\n")
  cat("current fdr relevant interactions:", nrow(ceRNA_interactions_fdr))
  ceRNA_interactions_fdr <- ceRNA_interactions_fdr[order(ceRNA_interactions_fdr$p.adj),]
}

dir.create("circRNA/plots", recursive = T)
dir.create("total/plots", recursive = T)

# save R objects
save.image(file = file.path(out, "sponge_final.RData"))

# MOST SIGNIFICANT SPONGES
network_centralities <- sponge_node_centralities(ceRNA_interactions_fdr)
ceRNA_interactions_weight <- ceRNA_interactions_fdr
ceRNA_interactions_weight$weight <- -log10(ceRNA_interactions_fdr$p.adj)
weighted_network_centralities <- sponge_node_centralities(ceRNA_interactions_weight)
weighted_network_plot <- sponge_plot_network_centralities(weighted_network_centralities, top = 3)
png("total/plots/centralities.png")
plot(weighted_network_plot)
dev.off()

# CIRC RNA MRNA ONLY
circ.mRNA.only <- circ.mRNA.subnetwork(ceRNA_interactions_fdr, "hsa_circ_")
ceRNA_interactions_all_circ <- circ.mRNA.only
write.table(ceRNA_interactions_all_circ, "circRNA/circRNAs_as_ceRNAs.tsv", sep = "\t", row.names = F)

# NETWORK ANALYSIS CIRC
ceRNA_interactions_circ_weight <- ceRNA_interactions_all_circ
ceRNA_interactions_circ_weight$weight <- -log10(ceRNA_interactions_all_circ$p.val)
weighted_network_centralities_circ <- sponge_node_centralities(ceRNA_interactions_circ_weight)
# plot top n samples
n = 10
# betweeness
top_network_plot_btw <- sponge_plot_network_centralities(weighted_network_centralities_circ, top = n, measure = "btw")
png(file = "circRNA/plots/circ_btw.png")
plot(top_network_plot_btw)
dev.off()

# eigenvector
top_network_plot_ev <- sponge_plot_network_centralities(weighted_network_centralities_circ, top = n, measure = "ev")
png(file = "circRNA/plots/circ_ev.png")
plot(top_network_plot_ev)
dev.off()

# counts
top_network_plot_c <- sponge_plot_network_centralities(weighted_network_centralities_circ, top = n, measure = "count")
png(file = "circRNA/plots/circ_counts.png")
plot(top_network_plot_c)
dev.off()
stopCluster(cl) # stop cluster
# save R objects
save.image(file = file.path(out, "sponge_final.RData"))

#Create the network plots
ceRNA_interactions_fdr <- ceRNA_interactions_fdr[ceRNA_interactions_fdr$p.adj < 0.05,]
ceRNA_interactions_fdr_circ_1 <- ceRNA_interactions_fdr[grep("hsa_circ", ceRNA_interactions_fdr$geneB),]
ceRNA_interactions_fdr_circ_2 <- ceRNA_interactions_fdr[grep("hsa_circ", ceRNA_interactions_fdr$geneA),]
ceRNA_interactions_fdr_circ <- rbind(ceRNA_interactions_fdr_circ_1, ceRNA_interactions_fdr_circ_2)

sponge_plot_network(ceRNA_interactions_fdr_circ, genes_miRNA_candidates)
