#Compute correlations
details <- read.csv("meta/sample_details_1.csv", stringsAsFactors = FALSE, header = TRUE)

#Creating Object List 
details_up <- subset(details, details$genetic == "LRRK2-/SNCA-/GBA-")
details_up$filename <- paste0("P1_", details_up$PATNO, "_", details_up$EVENT_ID)
samples <- details_up$filename

#binding sites
pairBindSites <- read.table("results/miRanda/circRNA_miRNA_pairs.tsv", header = T, sep = "\t", stringsAsFactors = F, check.names = F)

#Import filtered miRNA normalized counts file
miRNA_expression <- read.csv("results/filtered_expression/PPMI_miRNA_logCPM.csv", header = T, stringsAsFactors = F, check.names = F, row.names = 1)

miRNA_expression <- miRNA_expression %>%
  rownames_to_column("miRNA")

miRNA_expression$miRNA <- gsub("\\.", "-", miRNA_expression$miRNA)

#import filteres circRNA normalized counts file
circRNA_expression <- read.csv("results/filtered_expression/PPMI_circ_CPM.csv", header = T, stringsAsFactors = F, check.names = F, row.names = 1)

#add annotation to normalized counts
circ_anno <- read.csv("annotation/PPMI_circ_annotations.csv", header = TRUE)
circ_anno <- circ_anno %>%
  column_to_rownames("ID")

colnames(circ_anno)
circ_anno <- circ_anno[, c(1:9, 12)]

#merge annotation folder with counts folder
circRNA_expression <- merge(circ_anno, circRNA_expression, by = 0, all = TRUE)
circRNA_expression <- circRNA_expression %>%
  column_to_rownames("Row.names")

#add unannotated labels to missing circBaseID 
circRNA_expression$CirBaseID[circRNA_expression$CirBaseID == ""] <- NA
empty_cells <- sum(is.na(circRNA_expression$CirBaseID))
replacements <- paste0("unannotated", "_", seq_len(empty_cells))

circRNA_expression$CirBaseID[is.na(circRNA_expression$CirBaseID)] <- replacements

if(length(samples) < 5){
  stop("Cannot perform correlation on less than 5 samples")
}

# filter expression for miRNA binding pairs
valid_circRNAs <- intersect(rownames(circRNA_expression), rownames(pairBindSites))
valid_miRNAs <- intersect(miRNA_expression$miRNA, colnames(pairBindSites))

circRNA_expression <- circRNA_expression[valid_circRNAs,]
write.csv(circRNA_expression, file = "~/Desktop/CircRNA-Project/Post_review/circRNA-sponging_pipeline/filtered_expression/circRNA_correlation_input.csv")

miRNA_expression <- miRNA_expression[miRNA_expression$miRNA %in% valid_miRNAs,]
write.csv(miRNA_expression, file = "~/Desktop/CircRNA-Project/Post_review/circRNA-sponging_pipeline/filtered_expression/miRNA_correlation_input.csv")

# check for annotation
annotation <- "CirBaseID" %in% colnames(circRNA_expression)

#header <- "circRNA\tmiRNA\tcircRNA_miRNA_ratio\tmiRNA_binding_sites\tpearson_R\tcorr_pval\tRSS_norm\tintercept\tintercept_pval\tslope\tslope_pval\tadj_r_squared"
#write(header, file=paste0("filtered_circRNA_miRNA_correlation_libSizeEstNorm_directwritten.tsv"), append = F)

#Change rownames of pairbind sites by CircBaseID
rows <- c("hsa_circ_0000118", "hsa_circ_0023936", "hsa_circ_0000566", 
         "hsa_circ_0109315", "hsa_circ_0000994","hsa_circ_0001481")

rownames(pairBindSites) <- rows


miRNA_for_row <- function(miRNA_expr_line, circRNA, circRNA_counts){
  mirna <- as.character(miRNA_expr_line[1])
  message("processing:", mirna)
  # get sample counts for current miRNA
  miRNA_counts <- miRNA_expr_line[-1]
  
  miRNA_counts <- data.frame(sample = as.character(names(miRNA_counts)), "miRNA_counts" = as.numeric(unname(miRNA_counts)))
  # compute circRNA expression vs. miRNA expression
  joined_counts <- merge(miRNA_counts, circRNA_counts, by="sample")
  
  # analyse circRNA/miRNA ratio
  mean_circRNA_counts <- mean(joined_counts$circRNA_counts)
  mean_miRNA_counts <- mean(joined_counts$miRNA_counts)
  circRNA_miRNA_ratio <- mean_circRNA_counts/mean_miRNA_counts
  
  # compute number of miRNA binding sites on circRNA
  #mirna <-miRNA
  #binding_sites <- nrow(bindsitDT[miRNA == mirna & Target == circRNA])
  
  # pair <- pairBindSites[miRNA == mirna & Target == circRNA]
  binding_sites <- pairBindSites[circRNA, mirna]
  
  # compute circRNA-miRNA correlation for all samples
  cor_res <- cor.test(joined_counts$miRNA_counts, joined_counts$circRNA_counts,  method = "pearson", use = "complete.obs")
  corr_R <- as.numeric(as.character(cor_res$estimate))
  corr_pval <- as.numeric(as.character(cor_res$p.value))
  
  # compute linear regression
  regression_model <- lm(miRNA_counts~circRNA_counts, data = joined_counts)
  intercept <- summary(regression_model)$coefficients[1,1]
  intercept_pval <- summary(regression_model)$coefficients[1,4]
  slope <- summary(regression_model)$coefficients[2,1]
  slope_pval <- summary(regression_model)$coefficients[2,4]
  adj_r_squared <- summary(regression_model)$adj.r.squared
  
  # compute residuals sum of squares
  # normalize counts for residuals sum of squares
  normalized_counts <- joined_counts[,c("circRNA_counts", "miRNA_counts")]
  min_circRNA_counts <- min(normalized_counts$circRNA_counts)
  max_circRNA_counts <- max(normalized_counts$circRNA_counts)
  normalized_counts[,"circRNA_counts"] <- (normalized_counts[,"circRNA_counts"] - min_circRNA_counts)/(max_circRNA_counts - min_circRNA_counts)
  min_miRNA_counts <- min(normalized_counts$miRNA_counts)
  max_miRNA_counts <- max(normalized_counts$miRNA_counts)
  normalized_counts[,"miRNA_counts"] <- (normalized_counts[,"miRNA_counts"] - min_miRNA_counts)/(max_miRNA_counts - min_miRNA_counts)
  norm_reg_model <- lm(miRNA_counts~circRNA_counts, data = normalized_counts)
  RSS_norm <- sum(norm_reg_model$residuals^2)
  
  res <- data.frame(circRNA = as.character(circRNA), miRNA = as.character(mirna), 
                    circRNA_miRNA_ratio = as.numeric(circRNA_miRNA_ratio), 
                    miRNA_binding_sites = as.numeric(binding_sites), 
                    pearson_R = as.numeric(corr_R), corr_pval = as.numeric(corr_pval), 
                    RSS_norm = RSS_norm, intercept = intercept, 
                    intercept_pval = intercept_pval, slope = slope, 
                    slope_pval = slope_pval, adj_r_squared = adj_r_squared)
  # write correlation info in file
  #write.table(res, file=paste0("filtered_circRNA_miRNA_correlation_libSizeEstNorm_directwritten.tsv"), sep = "\t", quote = F, row.names = F, append = T, col.names = F)
  
  return(res)
}

circRNA_for_row <- function(circRNA_expr_line){
  # get coordinations of current circRNA
  chr <- as.character(circRNA_expr_line[1])
  start <- as.numeric(as.character(circRNA_expr_line[2]))
  end <- as.numeric(as.character(circRNA_expr_line[3]))
  strand <- as.character(circRNA_expr_line[7])
  circRNA <- paste(chr,":", start, "|", end, sep="")
  
  
  # annotate if possible
  if(annotation){
    an = circRNA_expression[circRNA,"CirBaseID"]
    if (an != "None") {
      message("processing: ", circRNA, " (", an, ")")
      circRNA = an
    }
  } else {
    message("processing: ", circRNA)
  }
  # extract pure counts only
  circRNA_counts <- circRNA_expr_line[samples]
  
  # get sample counts for current circRNA
  circRNA_counts <- data.frame(sample = as.character(names(circRNA_counts)), "circRNA_counts" = as.numeric(unname(circRNA_counts)))
  
  res_list <- apply(miRNA_expression, 1, FUN = miRNA_for_row, circRNA, circRNA_counts)
  res_df <- do.call(rbind, res_list)
  return(res_df)
}

correlations_list <- apply(circRNA_expression, MARGIN = 1, circRNA_for_row)
correlations_df <- do.call(rbind, correlations_list)
correlations_df$adj_pval <- p.adjust(correlations_df$corr_pval, method = "BH")
#write.table(correlations_df, file=paste0("filtered_circRNA_miRNA_correlation_v2.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)

