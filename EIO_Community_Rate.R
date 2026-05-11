################################################################################
## Data import and compositional data transformations
################################################################################
library(zCompositions)
library(compositions)
library(dplyr)
library(DESeq2)
library(phyloseq)
library(ggplot2)

# Read Bacteria ASV table (noRares)
ASVbact <- read.csv("./Guting/BacteriaASV.csv")

# Assign ASV column to row names
rownames(ASVbact) <- as.character(unlist(ASVbact[, 1]))
ASVbact <- ASVbact[, -1]
# Imputation of zero values using Geometric Bayesian multiplicative replacement
bact.asv.ZeroRepl <- cmultRepl(t(ASVbact), label = 0, method = "CZM", output = "p-counts")

# Centred log ratio transformation
bact.asv.clrTrans <- apply(t(bact.asv.ZeroRepl), 2, function(x) {
  log(x) - mean(log(x))
})

# Transpose for WGCNA
bact.asv.clean = t(bact.asv.clrTrans)

# ASV-level Clustering 
BACT_ASV_tax <- read.csv("./Guting/BacteriaASVtaxonomy.csv") 
rownames(BACT_ASV_tax) <- as.character(unlist(BACT_ASV_tax[, 1]))
BACT_ASV_tax <- as.matrix(BACT_ASV_tax[, -1])

# Combine BACT 16S data into phyloseq object
ASVbact = otu_table(bact.asv.clrTrans, taxa_are_rows = TRUE)
TAXbact = tax_table(BACT_ASV_tax)
BACTseq_asv = phyloseq(ASVbact, TAXbact)
BACTseq_asv <- subset_taxa(BACTseq_asv, Order!="Chloroplast")

bact.asv.clean <- t(as(otu_table(BACTseq_asv), "matrix"))

bactSelect_asv <- subset_taxa(BACTseq_asv)

select.asv.clean = t(as(otu_table(bactSelect_asv), "matrix"))
write.csv(select.asv.clean, file = "./Guting/select.asv.clean.csv", row.names = TRUE)

# Read EIO Metadata
## Samples MUST be in the same order as in your OTU/ASV tables
EIO <- read.csv("./Guting/N2OMetadata.csv", header = TRUE, sep = ",")
rownames(EIO) <- EIO$Sample_ID

# Select relevant subset of variables for PLSR and WGCNA
EIO.trim <- EIO %>% 
  dplyr::select(O2Mol, N2OMol, delta_N2O, NO2, NH4_Prod, NO3_Prod, NO2_Prod, Urea_Prod)

# Subset Bacteria / Metadata to match samples with Archaeal reads (n=18)
common <- intersect(rownames(EIO.trim), rownames(select.asv.clean))
EIO.trim.common <- subset(EIO.trim, rownames(EIO.trim) %in% common)

bact.asv.clean.common <- subset(bact.asv.clean, rownames(bact.asv.clean) 
                                %in% common)
select.asv.clean.common <- subset(select.asv.clean, rownames(select.asv.clean) 
                                %in% common)

EIO.select.asv.comb <- cbind(select.asv.clean.common)
EIO.asv.comb <- cbind(bact.asv.clean.common)

write.csv(EIO.asv.comb, file = "./Guting/EIO_asv_comb.csv", row.names = TRUE)


################################################################################
## Weighted Gene Correlational network analyses
################################################################################
library(WGCNA)

data = EIO.asv.comb
data <- as.matrix(data)
# Take a quick look at what is in the data set:
dim(data);
names(data);

# Assign data to new working variable to streamline code
datExpr0 = as.data.frame(data);

# Check for genes and samples with too many missing values:
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
# If the last statement returns TRUE, all genes have passed the cuts. 
# If not, we remove the offending genes and samples from the data:
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", 
                     paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:", 
                     paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
# Cluster samples to check for obvious outliers.
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))

pdf("./Guting/sample_clustering_plot.pdf", width = 10, height = 8) 
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", 
     xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

# No outliers apparent so no removal steps necessary
datExpr = datExpr0
# Assign trimmed metadata to sample traits dataframe
datTraits = EIO.trim.common
datTraits <- as.matrix(datTraits)


datTraits.scale <- as.data.frame(scale(datTraits, center = TRUE, scale = TRUE))

dim(datTraits.scale)
names(datTraits.scale)

# We now have the expression data in the variable datExpr, and the scaled 
# environmental traits in the variable datTraits. Before we continue with 
# network construction and module detection, we visualize how the sample traits
# relate to the sample dendrogram.

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, 
# grey means missing entry
traitColors = numbers2colors(datTraits.scale, signed = TRUE);
# Plot the sample dendrogram and the colors underneath.
pdf("./Guting/Sample_dendrogram_and_trait_heatmap.pdf", width = 10, height = 8)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits.scale),
                    main = "Sample dendrogram and trait heatmap")  
dev.off()
# In the plot, white means a low value, red a high value, and 
# grey a missing entry.The last step is to save the relevant expression and 
# trait data for use in the next steps of the tutorial.
save(datExpr, datTraits.scale, file = "./Guting/WGCNA_dataInput.RData")

# Load the data saved in the first part
lnames <- load(file = "./Guting/WGCNA_dataInput.RData")


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


# Set the soft power threshold
softPower = 8;
adjacency = adjacency(datExpr, power = softPower, type = "signed");

# Turn adjacency into topological overlap 
TOM = TOMsimilarity(adjacency, TOMType = "signed");
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
pdf("./Guting/Gene_clustering_on_TOM-based_dissimilarity.pdf", width = 10, height = 8)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off()
#  Set the minimum module size relatively high: 20 taxa
minModuleSize = 80;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")


# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# Set similarity threshold for combining modules
MEDissThres = 0.50
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, 
                          verbose = 3)
# The merged module colors
mergedColors = merge$colors;
table(mergedColors)
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;


sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "./Guting/WGCNA-networkConstruction-auto.RData")

write.csv(MEs, file = "./Guting/MEs.csv", row.names = TRUE)
module_data <- data.frame(moduleLabels, moduleColors)
write.csv(module_data, file = "./Guting/ModuleLabels_Colors.csv", row.names = TRUE)

pdf("./Guting/GeneTree.pdf", width = 10, height = 8)
plot(geneTree, main = "Gene Tree", xlab = "Samples", sub = "Clustered Genes")
dev.off()


# Load network data saved in the second part.
lnames = load(file = "./Guting/WGCNA-networkConstruction-auto.RData");
lnames

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

pdf("./Guting/Module-trait_relationships.pdf", width = 10, height = 8)
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()
# Extract the taxa-specific correlation coefficients for each sample trait

# Isolate correlation coefficients for Urea->N2O
trait.urea = as.data.frame(datTraits$Urea_Prod);
names(trait.urea) = "Trait"
# names (colors) of the modules
modNames = substring(names(MEs), 4)
taxaModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(taxaModuleMembership), 
                                          nSamples));

names(taxaModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
taxaTraitSignificance = as.data.frame(cor(datExpr, trait.urea, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(taxaTraitSignificance), 
                                          nSamples));
names(taxaTraitSignificance) = paste("GS.", names(trait.urea), sep="");
names(GSPvalue) = paste("p.GS.", names(trait.urea), sep="");

module.member.traits.urea <- cbind(taxaModuleMembership, taxaTraitSignificance, 
                                    GSPvalue, moduleColors)

# Isolate coefficients for NO3->N2O
trait.NO3Prod = as.data.frame(datTraits$NO3_Prod);
names(trait.NO3Prod) = "Trait"
# names (colors) of the modules
modNames = substring(names(MEs), 4)
taxaModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(taxaModuleMembership), nSamples));

names(taxaModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
taxaTraitSignificance = as.data.frame(cor(datExpr, trait.NO3Prod, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(taxaTraitSignificance), nSamples));
names(taxaTraitSignificance) = paste("GS.", names(trait.NO3Prod), sep="");
names(GSPvalue) = paste("p.GS.", names(trait.NO3Prod), sep="");

module.member.traits.NO3Prod <- cbind(taxaModuleMembership, taxaTraitSignificance,
                                     GSPvalue, moduleColors)

# Isolate coefficients for NH4->N2O
trait.AmmProd = as.data.frame(datTraits$NH4_Prod);
names(trait.AmmProd) = "Trait"
# names (colors) of the modules
modNames = substring(names(MEs), 4)
taxaModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(taxaModuleMembership), nSamples));

names(taxaModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
taxaTraitSignificance = as.data.frame(cor(datExpr, trait.AmmProd, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(taxaTraitSignificance), nSamples));
names(taxaTraitSignificance) = paste("GS.", names(trait.AmmProd), sep="");
names(GSPvalue) = paste("p.GS.", names(trait.AmmProd), sep="");

module.member.traits.AmmProd <- cbind(taxaModuleMembership, taxaTraitSignificance, 
                                     GSPvalue, moduleColors)
# Isolate coefficients for NO2->N2O
trait.NO2Prod = as.data.frame(datTraits$NO2_Prod);
names(trait.NO2Prod) = "Trait"
# names (colors) of the modules
modNames = substring(names(MEs), 4)
taxaModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(taxaModuleMembership), nSamples));

names(taxaModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
taxaTraitSignificance = as.data.frame(cor(datExpr, trait.NO2Prod, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(taxaTraitSignificance), nSamples));
names(taxaTraitSignificance) = paste("GS.", names(trait.NO2Prod), sep="");
names(GSPvalue) = paste("p.GS.", names(trait.NO2Prod), sep="");

module.member.traits.NO2Prod <- cbind(taxaModuleMembership, taxaTraitSignificance,
                                     GSPvalue, moduleColors)

# Isolate coefficients for deltaN2O
trait.delta = as.data.frame(datTraits$delta_N2O);
names(trait.delta) = "Trait"
# names (colors) of the modules
modNames = substring(names(MEs), 4)
taxaModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(taxaModuleMembership), nSamples));

names(taxaModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
taxaTraitSignificance = as.data.frame(cor(datExpr, trait.delta, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(taxaTraitSignificance), nSamples));
names(taxaTraitSignificance) = paste("GS.", names(trait.delta), sep="");
names(GSPvalue) = paste("p.GS.", names(trait.delta), sep="");

module.member.traits.delta <- cbind(taxaModuleMembership, taxaTraitSignificance, 
                                    GSPvalue, moduleColors)

# Now we are going to construct a singal data frame that contains correlation
# coefficients, p-values, connectivity measures and taxonomic information for 
# each ASV contained in the WGCNA analyses
Alldegrees1=intramodularConnectivity(adjacency, mergedColors)
head(Alldegrees1)

connectivity <- cbind(Alldegrees1, mergedColors,taxaModuleMembership)
connectivity <- tibble::rownames_to_column(connectivity, "ASV")
connectivity <- connectivity %>% 
  rename(Subnetwork = mergedColors)
# Rename moduleColors to Subnetworks
connectivity$mergedColors[connectivity$mergedColors == "brown"] <- "SNET1"
connectivity$mergedColors[connectivity$mergedColors == "blue"] <- "SNET2"
connectivity$mergedColors[connectivity$mergedColors == "turquoise"] <- "SNET3"
connectivity$mergedColors[connectivity$mergedColors == "yellow"] <- "SNET4"

connectivity.traits <- cbind(connectivity, 
                             module.member.traits.urea$GS.Trait,
                             module.member.traits.urea$p.GS.Trait,
                             module.member.traits.delta$GS.Trait,
                             module.member.traits.delta$p.GS.Trait,
                             module.member.traits.NO3Prod$GS.Trait,
                             module.member.traits.NO3Prod$p.GS.Trait,
                             module.member.traits.AmmProd$GS.Trait,
                             module.member.traits.AmmProd$p.GS.Trait,
                             module.member.traits.NO2Prod$GS.Trait,
                             module.member.traits.NO2Prod$p.GS.Trait)

connectivity.traits <- connectivity.traits %>% 
  dplyr::rename(
    "GS.urea" = "module.member.traits.urea$GS.Trait",
    "p.urea" = "module.member.traits.urea$p.GS.Trait",
    "GS.Delta" =  "module.member.traits.delta$GS.Trait",
    "p.Delta" =  "module.member.traits.delta$p.GS.Trait",
    "GS.NO3Prod" = "module.member.traits.NO3Prod$GS.Trait",
    "p.NO3Prod" = "module.member.traits.NO3Prod$p.GS.Trait",
    "GS.AmmProd" = "module.member.traits.AmmProd$GS.Trait",
    "p.AmmProd" = "module.member.traits.AmmProd$p.GS.Trait",
    "GS.NO2Prod" = "module.member.traits.NO2Prod$GS.Trait",
    "p.NO2Prod" = "module.member.traits.NO2Prod$p.GS.Trait"
    )

# Include taxonomy information in data frame
# Match datasets by rows
tax <- as.data.frame(rbind(BACT_ASV_tax))
tax <- tibble::rownames_to_column(tax, "ASV")

## Create Data Vector of taxonomies for heatmap labelling
tax <-rbind(BACT_ASV_tax)
tax <- tibble::rownames_to_column(as.data.frame(tax), "ASV")

connectivity.traits.tax <- merge(connectivity.traits, tax[, c("ASV", "Phylum", "Class","Order",
        "Family", "Genus", "Tax")], by="ASV")
rownames(connectivity.traits.tax) <- as.character(unlist(connectivity.traits.tax[, 1]))
write.csv(connectivity.traits.tax, "./Guting/TaxaTraitRelationships.csv")

snet1.col <- "brown"  
snet2.col <- "blue" 
snet3.col <- "turquoise"  
snet4.col <- "yellow"  

# Select all SAR324 and Rhodobacteraceae reads


sar324 <- subset_taxa(BACTseq_asv, Order=="SAR324_clade")
sar324 = as(otu_table(sar324), "matrix")

Rhodobacteraceae <- subset_taxa(BACTseq_asv, Family=="Rhodobacteraceae")
Rhodobacteraceae = as(otu_table(Rhodobacteraceae), "matrix")

names1 = rownames(sar324)
names2 = rownames(Rhodobacteraceae)


library(grid)
# Create connectivity plots
#Urea-N2O
plot1 <- ggplot(connectivity.traits.tax) +
  geom_smooth(aes(MMblue, GS.urea), method = "lm", color="darkgrey", se=FALSE) +
  geom_point(aes(MMblue, GS.urea, colour=mergedColors,
                 size=kWithin), alpha = 0.6) +
  scale_color_brewer(labels = c("SNET1", "SNET2", "SNET3", "SNET4"),
                     name = "mergedColors", palette = "Dark2",
                     direction = -1) +
  scale_size_continuous(range = c(0.2, 10),
                        name=expression(bold("K"["in"]))) +
  labs(y = expression("ASV importance (Urea)"), x = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.title = element_text(face = "plain"),
        legend.position.inside = c(0.85, 0.15),   # 修改为 legend.position.inside
        legend.background = element_rect(color="black"),
        legend.key.size = unit(2.5, "mm"))
plot1
ggsave("./Guting/Kin.pdf", plot = plot1 + theme(legend.position = "right"), 
       width = 7, height = 5.2, units = "in")
plot1 <- plot1 +  
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names1), ],
             aes(MMblue, GS.urea), shape = 2, size = 3) +
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names2), ],
             aes(MMblue, GS.urea), shape = 4, size = 3) +
  guides(color="legend", size = "none")
plot1

ggsave("./Guting/urea.pdf", plot = plot1 + theme(legend.position = "right"), 
       width = 7, height = 5.2, units = "in")

#Ammun-N2O
plot2 <- ggplot(connectivity.traits.tax) +
  geom_smooth(aes(MMblue, GS.AmmProd), method = "lm", color="darkgrey", se=FALSE) +
  geom_point(aes(MMblue, GS.AmmProd, colour=mergedColors,
                 size=kWithin), alpha = 0.6) +
  scale_color_brewer(labels = c("SNET1", "SNET2", "SNET3", "SNET4"),
                     name = "mergedColors", palette = "Dark2",
                     direction = -1) +
  scale_size_continuous(range = c(0.2, 10),
                        name=expression(bold("K"["in"]))) +
  labs(y = expression("ASV importance (Ammun)"), x = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.title = element_text(face = "plain"),
        legend.position.inside = c(0.85, 0.15),   # 修改为 legend.position.inside
        legend.background = element_rect(color="black"),
        legend.key.size = unit(2.5, "mm"))
plot2

plot2 <- plot2 +  
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names1), ],
             aes(MMblue, GS.AmmProd), shape = 2, size = 3) +
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names2), ],
             aes(MMblue, GS.AmmProd), shape = 4, size = 3) +
  guides(color="legend", size = "none")
plot2

ggsave("./Guting/Ammun.pdf", plot = plot2 + theme(legend.position = "right"), 
       width = 7, height = 5.2, units = "in")

#NO2-N2O
plot3 <- ggplot(connectivity.traits.tax) +
  geom_smooth(aes(MMblue, GS.NO2Prod), method = "lm", color="darkgrey", se=FALSE) +
  geom_point(aes(MMblue, GS.NO2Prod, colour=mergedColors,
                 size=kWithin), alpha = 0.6) +
  scale_color_brewer(labels = c("SNET1", "SNET2", "SNET3", "SNET4"),
                     name = "mergedColors", palette = "Dark2",
                     direction = -1) +
  scale_size_continuous(range = c(0.2, 10),
                        name=expression(bold("K"["in"]))) +
  labs(y = expression("ASV importance (NO2)"), x = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.title = element_text(face = "plain"),
        legend.position.inside = c(0.85, 0.15),   # 修改为 legend.position.inside
        legend.background = element_rect(color="black"),
        legend.key.size = unit(2.5, "mm"))
plot3

plot3 <- plot3 +  
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names1), ],
             aes(MMblue, GS.NO2Prod), shape = 2, size = 3) +
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names2), ],
             aes(MMblue, GS.NO2Prod), shape = 4, size = 3) +
  guides(color="legend", size = "none")
plot3

ggsave("./Guting/NO2.pdf", plot = plot3 + theme(legend.position = "right"), 
       width = 7, height = 5.2, units = "in")
#NO3-N2O
plot4 <- ggplot(connectivity.traits.tax) +
  geom_smooth(aes(MMyellow, GS.NO3Prod), method = "lm", color="darkgrey", se=FALSE) +
  geom_point(aes(MMyellow, GS.NO3Prod, colour=mergedColors,
                 size=kWithin), alpha = 0.6) +
  scale_color_brewer(labels = c("SNET1", "SNET2", "SNET3", "SNET4"),
                     name = "mergedColors", palette = "Dark2",
                     direction = -1) +
  scale_size_continuous(range = c(0.2, 10),
                        name=expression(bold("K"["in"]))) +
  labs(y = expression("ASV importance (NO3)"), x = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.title = element_text(face = "plain"),
        legend.position.inside = c(0.85, 0.15),   # 修改为 legend.position.inside
        legend.background = element_rect(color="black"),
        legend.key.size = unit(2.5, "mm"))
plot4

plot4 <- plot4 +  
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names1), ],
             aes(MMyellow, GS.NO3Prod), shape = 2, size = 3) +
  geom_point(data=connectivity.traits.tax[(rownames(connectivity.traits.tax) %in% names2), ],
             aes(MMyellow, GS.NO3Prod), shape = 4, size = 3) +
  guides(color="legend", size = "none")
plot4

ggsave("./Guting/NO3.pdf", plot = plot4 + theme(legend.position = "right"), 
       width = 7, height = 5.2, units = "in")

