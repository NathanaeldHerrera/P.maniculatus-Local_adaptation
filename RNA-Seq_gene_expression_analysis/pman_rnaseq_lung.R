### Load Packages ###########################
rm(list=ls())
message = FALSE # hide package loading messages
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(gprofiler2)
library(plotrix)
library(edgeR) 
library(DESeq2)
library(WGCNA)
library(lme4)
library(plyr)
library(ggsignif)
#==============================================================================
############ Read in raw count data #####################################
run_date <- format(Sys.Date(), "%d%b%Y")
counts0 <- read.csv("./data_files/lung/Pman_readcounts_LUNG_CountMM.csv",header=T) #load the counts.text file from featureCounts

counts <- counts0[, -c(2:6)] # remove the unncessary columns
rownames(counts) <- counts$Geneid #make the rownames the gene id's
counts$Geneid <- NULL #remove the gene id column, since you only want data columns in your matrix

#write.csv(counts, file="./data_files/lung/Pman_LUNG_counts_removeCol.csv")

### filter samples and low reads for PCA plot #####################################
# import metadata
id <- read.csv("./data_files/lung/Pman_rnaseq_LUNG_metadata.csv",header=T)
head(id)
head(counts)
dim(id)
#STOP: make sure the column names in the counts dataset are ordered by sample ID
colnames(counts) == id$code # should see all "TRUE"s 

dim(counts) # verify how many genes (first value) and samples (second value) we have data for

#remove 7 outlier samples identified from technical PCA and clustering tree for LUNG

id <- id[id$code!=c("AC1KNO"),]
id <- id[id$code!=c("CC2KNF"),]
id <- id[id$code!=c("AA4KNV"),]
id <- id[id$code!=c("CC1KNN"),]
id <- id[id$code!=c("CC4MEX"),]
id <- id[id$code!=c("CC3KNN"),]
id <- id[id$code!=c("AA1KNO"),]

dim(id)

# Subset count data to remove outlier individuals for LUNG
counts <- subset(counts, select=-c(AC1KNO,CC2KNF,AA4KNV,CC1KNN,CC4MEX,CC3KNN,AA1KNO))
dim(counts)
#colnames(counts) == id$code # should see all "TRUE"s

## Filter the data
#filter out genes that have <20 reads on average across individuals
counts$mean = rowMeans(counts) 
keep_counts = subset(counts, mean >= 20) 
#dim(keep_counts)

keep_counts$mean = NULL #clean up dataset
dim(keep_counts)
#==============================================================================
# use edgeR to normalize across libraries. 
id$code == colnames(keep_counts) # double check that the files are ordered in the same way
treatment <- id$treatment 
population <- id$population

table <- data.frame(population)
population <- factor(id$population)
treatment <- factor(id$treatment)

table <- cbind(table, treatment)
table$treatment = as.factor(table$treatment)
table$population <- relevel(factor(table$population), ref="Lowlander")
design <- model.matrix(~population*treatment, data=table)
colnames(design) <- c("(Intercept)","population","treatment","pop*treat") # change column names to be simpler

# do normalization 
y <- DGEList(counts=keep_counts, group=population) # make a DGE list
y <- calcNormFactors(y) # normalize
y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

#==============================================================================
# Generate MDS plot
pch <- c(0,1,2,16) #pch <- c(0,1,2,15,16)
colors<- (c("dark green","light green","dark blue","light blue"))

plotMDS(y,col=colors[population],pch=12,labels=colnames(y)) #exploration plots
legend("topleft",legend=levels(population),pch=pch,col=colors,ncol=2,cex=0.8)

pdf(file=paste("./results/lung/Rplot_MDS_n43_",run_date,".pdf",sep=""),width=8,height=8)
plotMDS(y,col=colors[group],pch=12,labels=colnames(y)) #exploration plots
legend("topright",legend=levels(group),pch=pch,col=colors,ncol=2,cex=0.8)
dev.off()
#==============================================================================
# CPM normalize the data
lung.norm <- cpm(y,log=TRUE,prior.count=2,normalized.lib.sizes=TRUE) #cpm normalized and log transformed expression dataset
lung.norm1 <- lung.norm
#transpose expression data for further analysis
Expr0 = as.data.frame(t(lung.norm1));
rownames(Expr0) = colnames(lung.norm1)
#check for genes and samples with too many missing values
gsg = goodSamplesGenes(Expr0, verbose = 0);
gsg$allOK
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(Expr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(Expr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  Expr0 = Expr0[gsg$goodSamples, gsg$goodGenes]
}
Expr = Expr0

#cluster the samples to check for outliers
sampleTree = hclust(dist(Expr), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)

par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
#abline(h = 180, col = "red");

pdf(file = paste("./results/lung/sampleClustering_n41_",run_date,".pdf",sep=""), width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()

#==============================================================================
# Choose a set of soft-thresholding powers
powers = c(c(1:12), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(Expr, powerVector = powers, networkType = "signed", verbose = 0)
# Plot the results:
pdf(paste("results/lung/signed/lung_beta_unsigned_plot_n41_",run_date,".pdf",sep=""), h=7, w=7)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#==============================================================================
# Construct a gene network, using a soft threshold power of 8, based on inflection point in plot.
Net = blockwiseModules(Expr, power =8, maxBlockSize = 18000, # set power according to threshold power plot
                       TOMType = "signed", networkType = "signed",
                       minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "ExprTOM",
                       #loadTOMs = TRUE,
                       verbose = 0,
                       randomSeed=46693) #set seed so that results are fully reproducible

table(Net$colors)
moduleLabels = Net$colors
moduleColors = labels2colors(Net$colors)
MEs = Net$MEs;
geneTree = Net$dendrograms[[1]];
table(moduleColors)
dim(table(moduleColors))

write.csv(table(moduleLabels), file=paste("results/lung/signed/lung_modules_signed",run_date,".csv",sep=""))

#==============================================================================
## Top hub genes from each module
hub_genes <-chooseTopHubInEachModule(Expr, moduleLabels,power=8)
write.csv(hub_genes, file=paste("results/lung/signed/lung_signed_hubGenes_",run_date,".csv",sep=""))

# Convert labels to colors for plotting
moduleColors = labels2colors(Net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(Net$dendrograms[[1]], moduleColors[Net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

pdf(file=paste("results/lung/signed/Rplot_moduleDendrogram_lung_signed",run_date,".pdf",sep=""),h=9,w=9)
plotDendroAndColors(Net$dendrograms[[1]], moduleColors[Net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#==============================================================================
#Analysis of phenotypic correlations with module expression in WGCNA. 
# load trait data
traitData <- id
# Form a data frame analogous to expression data that will hold the phenotypic traits.
Samples = rownames(Expr);
traitRows = match(Samples,traitData$code)# make sure IDs are in the same order as Expr dataset
Traits0 = traitData[traitRows, -1]; # removes Code column
rownames(Traits0) <-  na.omit(traitData[traitRows, 1]);
Traits0$code<- rownames(Traits0)
Traits <- Traits0[, 5]

collectGarbage();

# Trait Module Associations
# Define numbers of genes and samples
nGenes = ncol(Expr);
nSamples = nrow(Expr)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(Expr, moduleLabels)$eigengenes
MEs = orderMEs(MEs0) 
MEs1 <- MEs
moduleTraitCor = cor(MEs1, Traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
corr.table <- as.data.frame(cbind(moduleTraitCor, moduleTraitPvalue))
colnames(corr.table) <- c("cor", "p")
corr.table$p.adjust <- p.adjust(corr.table$p, method="fdr")
corr.table <- corr.table[order(corr.table$p),]

write.csv(corr.table, file=paste("./results/lung/signed/Pman.lung.corrTable_n41_",run_date,".csv",sep=""))

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(Expr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
# Create a dataset containing all gene-specific information
genes=names(Expr)
geneInfo0 = data.frame(Gene = genes,
                       moduleLabel = moduleLabels)
# Add module membership information in the chosen order
oldNames = names(geneInfo0)
geneInfo0 = data.frame(geneInfo0, geneModuleMembership,
                       MMPvalue);
names(geneInfo0) = c(oldNames, paste("MM.", modNames, sep=""),
                     paste("p.MM.", modNames, sep=""))

# Order the genes in the geneInfo variable first by module number, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleLabel);
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file=paste("results/lung/signed/pman.lung.norm_moduleExpression_signed_genes_",run_date,".csv",sep=""))

#module-trait relationship plot

pdf(paste("results/lung/signed/pman_lung_module_trait_relationships_plot_",run_date,".pdf",sep=""), h=11, w=8.5)
#sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(12, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = ("sao2"),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.45,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

# and a version for right here
sizeGrWindow(8,12)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(12, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = ("SaO2"),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.3,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#==============================================================================
# Line plots for treatment vs population effects 
ME3 <- ME_info[ , c("code","population","treatment","ME3")]

L03 <- ddply(ME3, c("population", "treatment"), summarise,
             N    = sum(!is.na(ME3)),
             mean = mean(ME3, na.rm=TRUE),
             sd   = sd(ME3, na.rm=TRUE),
             se   = sd / sqrt(N)
)

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.1) # move them .05 to the left and right

# A much nicer plot
ggplot(L03, aes(x=treatment, y=mean, colour=population, group=population)) + 
  scale_x_discrete(limits=c("Control","Acclimated")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(aes(shape = population), position=pd, fill="white", size=3) +
  scale_shape_manual(values = c(1, 16)) +
  xlab("L03") +
  ylab("Module Expression") +
  scale_color_manual(values=c("black","black")) +
  expand_limits(y=0) +                        # Expand y range
  theme_bw() +
#  theme(legend.position="none")             # for supplement we can exclude 
  theme(legend.justification=c(1,1),
        legend.position=c(1,1))               # Position legend in bottom right
# PDF 
pdf(file=paste("results/lung/signed/plots/eigengene_expr_linePlot_L03_",run_date,".pdf",sep=""),h=3.6,w=4.25) # h=3.6,w=4.25 for main fig; h=2.1,w=2.75 supplement
ggplot(L03, aes(x=treatment, y=mean, colour=population, group=population)) + 
  scale_x_discrete(limits=c("Control","Acclimated")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(aes(shape = population), position=pd, fill="white", size=3) +
  scale_shape_manual(values = c(1, 16)) +
  xlab("L03") +
  ylab("Module Expression") +
  scale_color_manual(values=c("black","black")) +
  expand_limits(y=0) +                        # Expand y range
  theme_bw() +
#  theme(legend.position="none")             # for supplement we can exclude 
  theme(legend.justification=c(1,1),
        legend.position=c(1,1))               # Position legend in bottom right
dev.off()
# Output scatterplots of module expression with phenotypes.
# We can view a plot here and make sure it looks like what we want for the output plots below.
ggplot(ME_info, aes(x=ME3, y=sao2, shape=population, color=treatment)) + scale_shape_manual(values = c(17, 16)) + 
  geom_point(size=3) + geom_smooth(method = lm, aes(group=1))  + theme_bw() +
  scale_color_manual(values=c("black","grey")) + 
  theme(axis.text.x = element_text(color = "black", size = 10, angle = 90, hjust = .5, vjust = .5, face = "plain"), axis.text.y = element_text(color = "black", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"))

#PDF plot The first one with a legend

pdf(file=paste("results/lung/signed/module_expression_scatter_plots/phenotype_lung_standarized_L03_",run_date,".pdf",sep=""),h=11, w=8.5)
ggplot(ME_info, aes(x=ME3, y=rv_standarized, shape=population, color=treatment)) + scale_shape_manual(values = c(17, 16)) + 
  geom_point(size=3) + geom_smooth(method = lm, aes(group=1))  + theme_bw() +
  scale_color_manual(values=c("black","grey")) + 
  theme(axis.text.x = element_text(color = "black", size = 10, angle = 90, hjust = .5, vjust = .5, face = "plain"), axis.text.y = element_text(color = "black", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"))
dev.off()
#==============================================================================
# Test for gene enrichment from WGCNA modules that had significant treatment effect in ANOVA.
all_genes <- geneInfo["Gene"]

# Run GO enrichment using gProfiler package
for (num in c("1","2","4","5","7","8","10","11","12","13","15","16","17","19","20","21","23","32","35","38","39")){
  genes <- geneInfo[geneInfo$moduleLabel==num,]["Gene"]
  genes_GO <- gost(as.vector(genes$Gene), organism = "pmbairdii",
                   ordered_query = F, significant = T, exclude_iea = F,
                   user_threshold=0.05,
                   correction_method = "g_SCS", custom_bg = as.vector(all_genes$Gene), 
                   domain_scope = "annotated", numeric_ns = "",sources=c("GO","KEGG","REAC","TF","HP"))
  write.csv(apply(genes_GO$result,2,as.character),file=paste("results/lung/signed/GO_analysis/new_072324/",num,"_module_GO_",run_date,".csv",sep=""))
}
# module 3, 6, ,9 , 14, 18, 22, 24,  25, 26, 27, 29, 30, 31, 33, 34,  36, and 37 have no significant GO terms
#==============================================================================
# Is there a difference in Male vs Female lung_standardized weight?

id %>%
  group_by(sex) %>%
  get_summary_stats(sao2, type = "mean_sd")

# We can plot this as a box plot (does seem to be a difference)
bxp <- ggboxplot(
  id, x = "sex", y = "sao2", 
  ylab = "sao2", xlab = "sex", add = "jitter"
)
bxp
boxplot

#Nicer plot  
ggplot(id, aes(x=sex, y=sao2)) +
  geom_boxplot() + geom_dotplot(binaxis='y', stackdir='center',
                                position=position_dodge(1), dotsize = 0.5) +
  geom_signif(comparisons = list(c("M", "F")), 
              map_signif_level=TRUE)

# Any outliers? No.
id %>%
  group_by(sex) %>%
  identify_outliers(sao2)

# Compute Shapiro wilk test by groups
id %>%
  group_by(sex) %>%
  shapiro_test(sao2)
# Draw a qq plot by group
ggqqplot(id, x = "sao2", facet.by = "sex")

stat.test <- id %>% 
  t_test(sao2 ~ sex) %>%
  add_significance()
stat.test

stat.test2 <- id %>%
  t_test(sao2 ~ sex, var.equal = TRUE) %>%
  add_significance()
stat.test2
