# compute paths
statistics_file <- paste0("sponging_statistics.txt")
plot_folder <- paste0("plots/")
dir.create(file.path(plot_folder))

# get dataset structure, samples, miRNA expression and circRNA expression
dataset <- read.csv("metadata/metadata.csv", header=T, stringsAsFactors = F)
samples <- dataset$filename

miRNA_expression_raw <- read.csv("filtered_expression/mature_counts.csv", header = T, stringsAsFactors = F, check.names = F)
miRNA_expression_raw <- t(miRNA_expression_raw)
circRNA_expression_raw <- read.csv("filtered_expression/circRNA_raw_expression.csv", header = T, stringsAsFactors = F, check.names = F)

miRNA_expression <- read.csv("filtered_expression/miRNA_correlation_input.csv", header = T, stringsAsFactors = F, check.names = F)
circRNA_expression <- read.csv("filtered_expression/circRNA_correlation_input.csv", header = T, stringsAsFactors = F, check.names = F)
rownames(circRNA_expression) <- circRNA_expression$circBaseID

# annotate expression
#circRNA_expression <- annotate(circRNA_expression, circRNA_expression)

# write starting statistics to file
file.create(statistics_file)
cat("Sponging statistics",file=statistics_file,append=TRUE, sep="\n")
cat(paste0("Number of filtered miRNAs used for the correlation analysis: ", nrow(miRNA_expression)," (from ", nrow(miRNA_expression_raw), " unfiltered)"),file=statistics_file,append=TRUE, sep="\n")
cat(paste0("Number of filtered circRNAs used for the correlation analysis: ", nrow(circRNA_expression)," (from ", nrow(circRNA_expression_raw), " unfiltered)"),file=statistics_file,append=TRUE, sep="\n")

# read correlation data
correlations <- data.frame(read.table("correlation/filtered_circRNA_miRNA_correlation.tsv", sep = "\t", stringsAsFactors = F, header = T))
cat(paste0("Initial number of circRNA-miRNA pairs (including NAs in column miRNA_binding_sites): ", nrow(correlations)),file=statistics_file,append=TRUE, sep="\n")

# removing pairs where number of binding sites is NA
correlations_no_NA <- correlations[!is.na(correlations$miRNA_binding_sites),]
cat(paste0("Number of circRNA-miRNA pairs (after removing NAs in column miRNA_binding_sites): ", nrow(correlations_no_NA)),file=statistics_file,append=TRUE, sep="\n")

# define function for plotting correlation distribution
plotCorrelationDistribution <- function(correlations_df, filter_criteria_string, plot_folder, plot_name){
  p <- ggplot(correlations_df, aes(x=pearson_R)) +
    geom_histogram(colour=I("orange"), fill=I("orange"), alpha=I(.2)) +
    labs(title = "Correlation distribution",
         subtitle = paste(nrow(correlations_df) ,"circRNA-miRNA pairs"),
         caption =filter_criteria_string,
         x = "Pearson correlation coefficient R",
         y = "circRNA-miRNA pairs") + xlim(c(-1.1,1.1))
  # add mean line
  y_coord <- max(ggplot_build(p)$data[[1]]$y)/2 # y coordinate for mean label
  p + geom_vline(xintercept = mean(correlations_df$pearson_R), linetype="dashed", 
                 color = "black", size=0.7) + geom_text(data = correlations_df, aes(x=mean(correlations_df$pearson_R), label=paste(round(mean(correlations_df$pearson_R), digits = 3), "\n"), y = y_coord), vjust = 1.25, angle=90)
  ggsave(paste0(plot_folder, "/", plot_name, ".png"), width = 4, height = 3)
}

# unfiltered results
correlations_processed <- correlations_no_NA
n_of_pairs_init <- nrow(correlations_processed)
n_of_pairs <- n_of_pairs_init
cat(paste0("Number of circRNA-miRNA pairs (unfiltered): ", n_of_pairs, " (", round(n_of_pairs/n_of_pairs_init*100,  digits=2), " %)"), file=statistics_file, append=TRUE, sep="\n")
plotCorrelationDistribution(correlations_processed, "unfiltered", plot_folder, paste0("correlation_distribution_unfiltered"))

###### NO BINDING SITES ######
# filter for circRNA-miRNA pairs having no common binding site
correlations_no_bind <- correlations_processed[correlations_processed$miRNA_binding_sites < 1,]
n_of_pairs <- nrow(correlations_no_bind)
cat(paste0("Number of circRNA-miRNA pairs (filter: no binding sites): ", n_of_pairs, " (", round(n_of_pairs/n_of_pairs_init*100,  digits=2), " %)"), file=statistics_file, append=TRUE, sep="\n")
plotCorrelationDistribution(correlations_no_bind, "Filter: no common binding sites", plot_folder, "correlation_distribution_noBindSites")

# filter for circRNA-miRNA pairs having no common binding site for significant p-value < 0.05
adj_pval_filter = 0.05
correlations_no_bind <- data.frame(correlations_no_bind[correlations_no_bind$adj_pval < adj_pval_filter,])
n_of_pairs <- nrow(correlations_no_bind)
cat(paste0("Number of circRNA-miRNA pairs (filter: no binding sites & correlation p-value < ", adj_pval_filter , "): ", n_of_pairs, " (", round(n_of_pairs/n_of_pairs_init*100,  digits=2), " %)"), file=statistics_file, append=TRUE, sep="\n")

###### WITH BINDING SITES ######
# filter for circRNA having mind. 1 binding site from that miRNA
bind_sites_filter = 1
correlations_processed <- correlations_processed[correlations_processed$miRNA_binding_sites >= bind_sites_filter,]
n_of_pairs <- nrow(correlations_processed)
cat(paste0("Number of circRNA-miRNA pairs (filter: number of binding sites >= ", bind_sites_filter,"): ", n_of_pairs, " (", round(n_of_pairs/n_of_pairs_init*100,  digits=2), " %)"), file=statistics_file, append=TRUE, sep="\n")
plotCorrelationDistribution(correlations_processed, paste0("Filter: binding sites > ", bind_sites_filter), plot_folder, paste0("correlation_distribution_minBindSites", bind_sites_filter))

# significant p-value < 0.05 for 1 bind. sites
bind_sites_filter = 1
adj_pval_filter = 0.05
correlations_processed <- correlations_processed[correlations_processed$miRNA_binding_sites >= bind_sites_filter,]
correlations_processed <- data.frame(correlations_processed[correlations_processed$adj_pval < adj_pval_filter,])
n_of_pairs <- nrow(correlations_processed)
cat(paste0("Number of circRNA-miRNA pairs (filter: number of binding sites >= ", bind_sites_filter," & correlation p-value < ", adj_pval_filter , "): ", n_of_pairs, " (", round(n_of_pairs/n_of_pairs_init*100,  digits=2), " %)"), file=statistics_file, append=TRUE, sep="\n")
plotCorrelationDistribution(correlations_processed, paste0("Filter: binding sites > ", bind_sites_filter), plot_folder, paste0("correlation_distribution_minBindSites", bind_sites_filter,"_",adj_pval_filter))

# number of miRNA binding sites > 3
bind_sites_filter = 3
correlations_processed_3 <- correlations_processed[correlations_processed$miRNA_binding_sites > bind_sites_filter,]
n_of_pairs <- nrow(correlations_processed_3)
cat(paste0("Number of circRNA-miRNA pairs (filter: number of binding sites >= ", bind_sites_filter," & correlation p-value < ", adj_pval_filter, "): ", n_of_pairs, " (", round(n_of_pairs/n_of_pairs_init*100,  digits=2), " %)"), file=statistics_file, append=TRUE, sep="\n")
plotCorrelationDistribution(correlations_processed, paste("Filter: adj_pval <", adj_pval_filter, "bind_sites >", bind_sites_filter), plot_folder, paste0("correlation_distribution_minBindSites_", bind_sites_filter, "_adj_pval", adj_pval_filter))

correlations_bind <- correlations_processed


# define function for plotting correlation for specific pair
plotCorrelationForPair <- function(circRNA, miRNA, circRNA_expression_df, miRNA_expression_df, bind_sites, R_value, adjusted_p_value, plot_folder, plot_name){
  name <- paste(circRNA, " VS. ", miRNA, sep = "")
  circRNA_counts <- data.frame(t(circRNA_expression_df[circRNA,samples]))
  colnames(circRNA_counts) <- "circRNA_counts"
  circRNA_counts$sample <- row.names(circRNA_counts)
  circRNA_counts$circRNA_counts <- as.numeric(as.character(circRNA_counts$circRNA_counts))
  
  # get sample counts for current miRNA
  miRNA_counts <- data.frame(t(miRNA_expression_df[miRNA_expression_df$miRNA == miRNA,]))
  miRNA_counts <- as.data.frame(miRNA_counts)
  colnames(miRNA_counts) <- "miRNA_counts"
  miRNA_counts$sample <- row.names(miRNA_counts)
  miRNA_counts <- miRNA_counts[-1, ]
  miRNA_counts$miRNA_counts <- as.numeric(as.character(miRNA_counts$miRNA_counts))
  
  # compute circRNA expression vs. miRNA expression
  joined_counts <- merge(circRNA_counts, miRNA_counts, by="sample")
  
  p <- ggplot(joined_counts, aes(x=circRNA_counts, y=miRNA_counts)) + 
    geom_point(size = 2)+
    geom_smooth(method = "lm", formula = y ~ x) +
    labs(title=name, 
         x ="circRNA counts", 
         y = "miRNA counts", 
         subtitle=paste0("R=",round(R_value, digits = 2),
                         ", bind-sites=", bind_sites,
                         ", p-adj=",round(adjusted_p_value, digits = 8)))
  
  p_labeled <- p + geom_label_repel(data = joined_counts, aes(label=sample), box.padding = 0.35, 
                                    point.padding = 0.5, segment.color = 'black', size = 2, 
                                    segment.size = 0.2, force = 100)
  ggsave(filename = paste0(plot_folder, plot_name,".png"), plot = p,
         width = 6, height = 4)
  ggsave(filename = paste0(plot_folder, plot_name,"_labeled.png"), plot = p_labeled,
         width = 6, height = 4)
  # 
  # if (sample_path != "null") {
  #   sample_structure <- read.table(sample_path, sep = "\t", header=T, stringsAsFactors = F)
  #   joined_counts <- merge(joined_counts, sample_structure, by="sample")
  #   
  #   p_colored <- ggplot(joined_counts, aes(x=circRNA_counts, y=miRNA_counts)) +
  #     geom_point(size = 2, aes(col = group))+
  #     geom_smooth(method = "lm", formula = y ~ x) +
  #     labs(title=name,
  #          x ="circRNA counts",
  #          y = "miRNA counts",
  #          subtitle=paste0("R=",round(R_value, digits = 2),
  #                          ", bind-sites=", bind_sites,
  #                          ", p-adj=",round(adjusted_p_value, digits = 8)))
  #   
  #   ggsave(filename = paste0(plot_folder, plot_name,"_groups.png"), plot = p_colored,
  #          width = 6, height = 4)
#  }
}

# plot top negative correlation
correlations_sign <- correlations_bind
correlations_sign <- correlations_sign[order(correlations_sign$pearson_R),]
top_plots <- list()

library(ggrepel)
for (i in 1:10){
  circRNA_min <- correlations_sign[i,1]
  miRNA_min <- correlations_sign[i,2]
  bind_sites <- correlations_sign[i,"miRNA_binding_sites"]
  R_value <- correlations_sign[i,"pearson_R"]
  adjusted_p_value <- correlations_sign[i,"adj_pval"]
  plotCorrelationForPair(circRNA_min, miRNA_min, 
                         circRNA_expression, miRNA_expression, 
                         bind_sites, R_value, adjusted_p_value, 
                         plot_folder, paste0("correlation_pair_", circRNA_min, "_", miRNA_min))
}