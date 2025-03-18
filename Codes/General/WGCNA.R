library(DESeq2)
library(WGCNA)

#Load the metadata
details <- read.csv("meta/sample_details_1.csv", header = TRUE)
details$filename <- paste0("P1_", details$PATNO, "_", details$EVENT_ID)
rownames(details) <- details$filename
details$PATNO <- as.factor(details$PATNO)
details$EVENT_ID <- factor(details$EVENT_ID, 
                           levels = c("BL", "V06", "V08"))


#Getting circRNAs counts 
gene_mtx <- read.csv("data/gene_count_matrix.csv", row.names = 1, header = TRUE)

#Selecting samples without genetic mutations associated with PD risk
details_up <- subset(details, details$genetic == "LRRK2-/SNCA-/GBA-" & details$APPRDX == "1")


#gene matrix with only counts for final selected dataset
gene_mtx_up <- gene_mtx[, details_up$filename]
dim(gene_mtx_up)

set.seed(123)
library("genefilter")
rv <- rowVars(gene_mtx_up)
summary(rv)

q75 <- quantile(rowVars(gene_mtx_up), .75)

gene_mtx_up_final <- gene_mtx_up[rv > q75, ]
summary(rowVars(gene_mtx_up_final))

#Obtaining the VST matrix
dds <- DESeqDataSetFromMatrix(countData = gene_mtx_up_final,
                              colData = details_up,
                              design = ~ PATNO + EVENT_ID)
dds

#Performing Variance Stabilizing Transformation to obtain normalized counts
vsd <- varianceStabilizingTransformation(dds)
norm_counts <- data.frame(t(assay(vsd)))
colnames(norm_counts) <- rownames(dds)

#save selected genes
selected_genes <- data.frame(ID = rownames(dds), ID_1 = rownames(dds))
selected_genes <- selected_genes %>%
  separate(col = ID_1, into = c("ENSEMBL", "Gene"), sep = "\\|")

#write.csv(selected_genes, file = "WGCNA/WGCNA_selected_genes.csv")

#Checking samples with missing values and outliers
datExpr0 <- norm_counts

gsg <- goodSamplesGenes(datExpr0,verbose=3)
gsg$allOK

#Clustering the samples
sampleTree <- hclust(dist(datExpr0), method="average")

#Plot the sample tree: Open a graphic output window of size 12 by 9 inches #The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,19) 


#pdf(file="Plots/sampleClustering.pdf",width=12,height=9)
par(cex=0.6)
par(mar=c(0,4,2,0)) 
plot(sampleTree, main="Sample clustering to detect outliers", sub="", xlab="", cex.lab= 1.5, cex.axis= 1.5, cex.main= 1.5, cex = 0.6)

#Plot aline to show the cut
abline(h= 80,col= "red")

#Determine cluster under the line
clust <- cutreeStatic(sampleTree, cutHeight = 80, minSize = 10)
table(clust) 

#clust 1 contains the samples we want to keep. 
keepSamples <- (clust==1) 
datExpr <- datExpr0[keepSamples,]
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

#Load the clinical metadata
traitData <- read.csv("meta/WGCNA_meta.csv", header = TRUE)
dim(traitData) 
names(traitData) 

#remove columns that hold information we do not need. 
allTraits <- traitData[,-c(1:2, 6:11,13,20)]
dim(allTraits)
names(allTraits) 

#Form a data frame analogous to expression data that will hold the clinical traits.
PD_Samples <- rownames(datExpr)
traitRows <- match(PD_Samples, allTraits$filename)
datTraits <- allTraits[traitRows, -12]
collectGarbage()

#Re-cluster samples
sampleTree2 <- hclust(dist(datExpr), method="average") 

#Convert traits to a color representation: white means low, red means high,grey means missing entry 
traitColors <- numbers2colors(datTraits, signed=FALSE) 

#Plot the sample dendrogram and the colors underneath. 
plotDendroAndColors(sampleTree2, traitColors, groupLabels=names(datTraits), cex.dendroLabels = 0.5, main="Sample dendrogram and trait heatmap")

#Choose a set of soft-thresholding powers 
powers <- c(c(1:10), seq(from = 12, to = 40, by = 2)) 

#Call the network topology analysis function 
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose=5, networkType ="signed")

#Plot the results: 
sizeGrWindow(9,5)
par(mfrow=c(1,2))
cex1 = 0.9

#Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     xlab = "Soft Threshold (power)", ylab ="Scale Free Topology Model Fit, signed R^2",
     type ="n", main = paste("Scale independence"))

text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, cex = cex1, col = "red")

#this line corresponds to using an R^2 cut-off of h
abline(h = 0.90,col = "red")

#Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1],
     sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n", 
     main = paste("Mean connectivity")) 

text(sft$fitIndices[,1],
     sft$fitIndices[,5],
     labels = powers, cex = cex1, col = "red")

#Module construction
#Using the soft threshold power to create adjacency matrix 
softPower <- 10
adjacency <- adjacency(datExpr, power = softPower)

#Converting the adjacency matrix to the topological overlap matrix (TOM) similarity matrix
TOM <- TOMsimilarity(adjacency, TOMType = "signed")

#Converting the similarity matrix to the dissimilarity
TOM.dissimilarity <- 1-TOM

#Performing hierarchial clustering analysis
geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average") 

#plotting the dendrogram
sizeGrWindow(19,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", 
     labels = FALSE, hang = 0.04)

Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)

table(Modules) #returns a table of the counts of factor levels in an object. In this case how many genes are assigned to each created module. 

set.seed(123)
ModuleColors <- labels2colors(Modules) #assigns each module number a color
table(ModuleColors) #returns the counts for each color (aka the number of genes within each module)

#plots the gene dendrogram with the module colors
plotDendroAndColors(geneTree, ModuleColors,"Module",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

#Identify module eigengenes on expression data
MElist <- moduleEigengenes(datExpr, colors = ModuleColors) 
MEs <- MElist$eigengenes 
head(MEs)

#Module Trait association
# Define numbers of genes and samples
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

module.trait.correlation <- cor(MEs, datTraits, use = "p")

module.trait.Pvalue <- corPvalueStudent(module.trait.correlation, nSamples) #calculate the p-value associated with the correlation

# Will display correlations and their p-values
textMatrix <- paste(signif(module.trait.correlation, 2), "\n(",
                    signif(module.trait.Pvalue, 1), ")", sep = "");
dim(textMatrix) <- dim(module.trait.correlation)
par(mar = c(6, 8.5, 3, 1))

#Rename the module colors (Remove ME)
colors <- names(MEs)
final_colors <- substr(colors, 3, nchar(colors))

# Display the correlation values within a heatmap plot
sizeGrWindow(25,25)
labeledHeatmap(Matrix = module.trait.correlation,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-Trait relationships"),
               cex.lab.x = 0.8,
               cex.lab.y = 0.6)
