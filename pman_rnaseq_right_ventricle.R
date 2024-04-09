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
counts0 <- read.csv("./data_files/rv/Pman_readcounts_RV_CountMM.csv",header=T) #load the counts.text file from featureCounts

counts <- counts0[, -c(2:6)] # remove the unncessary columns
rownames(counts) <- counts$Geneid #make the rownames the gene id's
counts$Geneid <- NULL #remove the gene id column, since you only want data columns in your matrix

#write.csv(counts, file="Pman_RV_counts_removeCol.csv")

#==============================================================================
## filter samples and low reads for PCA plot
# import metadata
id <- read.csv("./data_files/rv/Pman_rnaseq_RV_metadata_2.csv",header=T)
head(id)
head(counts)
dim(id)
#STOP: make sure the column names in the counts dataset are ordered by sample ID
colnames(counts) == id$code # should see all "TRUE"s 

dim(counts) # verify how many genes (first value) and samples (second value) we have data for

#remove 3 outlier samples identified from technical PCA and clustering tree
id <- id[id$code!=c("CC3MEBF"),]
id <- id[id$code!=c("KNG2304G"),]
id <- id[id$code!=c("AA4MEBF"),]
dim(id)
# Subsetcount data to remove outlier individuals
counts <- subset(counts, select=-c(CC3MEBF,KNG2304G,AA4MEBF))
#counts0 <- subset(counts0, select=-c(CC3MEBF,KNG2304G,AA4MEBF))
#dim(counts0)
dim(counts)
#colnames(counts) == id$code # should see all "TRUE"s

#==============================================================================
#filter data
counts$mean = rowMeans(counts) #rowMeans takes means of each row
keep_counts = subset(counts, mean >= 20) #filter out genes that have <20 reads on average across individuals
#dim(keep_counts)

keep_counts$mean = NULL #clean up dataset
dim(keep_counts)

#==============================================================================
# Next, we set up our GLM model and use edgeR to normalize across libraries. 
id$code == colnames(keep_counts) # double check that the files are ordered in the same way, other than the "mean" column in the counts df
treatment <- id$treatment #this is the first factor for the GLM model
population <- id$population

#table <- data.frame(treatment)
table <- data.frame(population)
population <- factor(id$population)
treatment <- factor(id$treatment)

table <- cbind(table, treatment)
table$treatment = as.factor(table$treatment)
table$treatment <- relevel(factor(table$treatment), ref="Control")#set Lowlander as the reference 
table$population <- relevel(factor(table$population), ref="Lowlander")
design <- model.matrix(~population*treatment, data=table)
colnames(design) <- c("(Intercept)","population","treatment","pop*treat") # change column names to be simpler

# do normalization 
y <- DGEList(counts=keep_counts, group=treatment) # make a DGE list
dim(y)
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

pdf(file=paste("./results/rv/Rplot_MDS_n43_",run_date,".pdf",sep=""),width=8,height=8)
plotMDS(y,col=colors[group],pch=12,labels=colnames(y)) #exploration plots
legend("topright",legend=levels(group),pch=pch,col=colors,ncol=2,cex=0.8)
dev.off()

#==============================================================================
rv.norm <- cpm(y,log=TRUE,prior.count=2,normalized.lib.sizes=TRUE) #cpm normalized and log transformed expression dataset
rv.norm1 <- rv.norm
#transpose expression data for further analysis
Expr0 = as.data.frame(t(rv.norm1));
rownames(Expr0) = colnames(rv.norm1)
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

pdf(file = paste("./results/rv/sampleClustering_n42_",run_date,".pdf",sep=""), width = 12, height = 9);
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
pdf(paste("results/rv/signed/rv_beta_unsigned_plot_n42_",run_date,".pdf",sep=""), h=7, w=7)
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
# Construct a gene network, using a soft threshold power of 12, based on inflection point in plot. 
Net = blockwiseModules(Expr, power =12, maxBlockSize = 18000, # set power according to threshold power plot
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

write.csv(table(moduleLabels), file=paste("results/rv/signed/rv_modules_signed",run_date,".csv",sep=""))

#==============================================================================
## Top hub genes from each module
hub_genes <-chooseTopHubInEachModule(Expr, moduleLabels,power=12)
write.csv(hub_genes, file=paste("results/rv/signed/rv_modules_signed_hubGenes_",run_date,".csv",sep=""))

# Convert labels to colors for plotting
moduleColors = labels2colors(Net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(Net$dendrograms[[1]], moduleColors[Net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

pdf(file=paste("results/rv/signed/Rplot_moduleDendrogram_rv_signed",run_date,".pdf",sep=""),h=9,w=9)
plotDendroAndColors(Net$dendrograms[[1]], moduleColors[Net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#==============================================================================
# Analysis of phenotypic correlations with module expression in WGCNA. 
# load trait data
traitData <- id
# Form a data frame analogous to expression data that will hold the phenotypic traits.
row_names_df_to_remove<-("6")
#Expr <- Expr[!(row.names(Expr) %in% row_names_df_to_remove),]
traitData <- traitData[!(row.names(traitData) %in% row_names_df_to_remove),]
Samples = rownames(Expr);
traitRows = match(Samples,traitData$code)# make sure IDs are in the same order as Expr dataset
Traits0 = traitData[traitRows, -1]; # removes Code column
rownames(Traits0) <-  na.omit(traitData[traitRows, 1]);
Traits0$code<- rownames(Traits0)
Traits <- Traits0[, 4]

collectGarbage();

#==============================================================================
# Trait Module Associations
# Define numbers of genes and samples
nGenes = ncol(Expr);
nSamples = nrow(Expr)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(Expr, moduleLabels)$eigengenes
MEs = orderMEs(MEs0) #the rownames of this dataset are equal to Expr
MEs1 <- MEs
#match(rownames(Traits), rownames(Expr))
moduleTraitCor = cor(MEs1, Traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
corr.table <- as.data.frame(cbind(moduleTraitCor, moduleTraitPvalue))
colnames(corr.table) <- c("cor", "p")
corr.table$p.adjust <- p.adjust(corr.table$p, method="fdr")
corr.table <- corr.table[order(corr.table$p),]

write.csv(corr.table, file=paste("./results/rv/signed/Pman.rv.norm_signed_CorrTable_n41_",run_date,".csv",sep=""))

# names (numbers) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(Expr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
# Create a dataset containing all gene-specific information
genes=names(Expr)
geneInfo0 = data.frame(Gene = genes,
                       moduleLabels = moduleLabels)
# Add module membership information in the chosen order
oldNames = names(geneInfo0)
geneInfo0 = data.frame(geneInfo0, geneModuleMembership,
                       MMPvalue);
names(geneInfo0) = c(oldNames, paste("MM.", modNames, sep=""),
                     paste("p.MM.", modNames, sep=""))

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleLabels);
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file=paste("results/rv/pman.rv.norm_moduleExpression_signed_genes_",run_date,".csv",sep=""))

#code to generate plot of module-trait relationships

pdf(paste("results/rv/signed/pman_rv_module_trait_relationships_plot_signed",run_date,".pdf",sep=""), h=11, w=8.5)
#sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(12, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = ("RV mass"),
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
               xLabels = ("RV mass"),
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
# ANOVAs on population and treatment effect
# code from Plachetzki et al 2014 (https://bitbucket.org/plachetzki/plachetzki-et-al.-sicb-2014/src/master/)
ME = MEs
ME <- removeGreyME(MEs, greyMEName = paste(moduleColor.getMEprefix(), "0", sep=""))
ME_t <- ME
ME_t$code <- row.names(ME_t)
ME_info <- merge(ME_t,id,by="code")
#ME_anova <- ME_info %>% select(-one_of('code','population','treatment','sex','sao2','hct'))

mods <- colnames(ME)
modNames = names(ME_info)

anovas<-list()
anovasP<-list()
for (i in 1:ncol(ME_info)){
  rME<-rank(ME_info[,i])
  anova<-summary(aov(rME ~ ME_info$population*ME_info$treatment))
  anovasP[[i]]<-c(anova[[1]][["Pr(>F)"]][1:3]); names(anovasP)[i]<-modNames[i]
  anovas[i]<-anova; names(anovas)[i]<-modNames[i]
}
FDR_pvals<-matrix(round(p.adjust(unlist(anovasP),method="fdr"),4),42,3,byrow=T,dimnames=list(names(ME_info),c("population","treatment", "population*treatment")) )            

anova<-summary(aov(ME_info$ME36 ~ ME_info$population*ME_info$treatment))

write.table(FDR_pvals,file=paste("results/rv/signed/table_rv_wgcna_signed_anova_results_",run_date,".txt",sep=""),sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)

#==============================================================================
# Box plots for module expression of population differences and treatment differences
# Plot for right here
ggplot(ME_info, aes(x=population, y=ME3)) +
  geom_boxplot() + geom_dotplot(binaxis='y', stackdir='center',
                                position=position_dodge(1), dotsize = 0.5) +
  geom_signif(comparisons = list(c("Lowlander", "Highlander")), 
              map_signif_level=TRUE,)
# Save as a pdf 
pdf(file=paste("results/rv/signed/module_expression_boxplots/population/eigengene_expr_ME10_",run_date,".pdf",sep=""),h=9,w=9)
ggplot(ME_info, aes(x=population, y=ME10)) + 
  geom_boxplot() + geom_dotplot(binaxis='y', stackdir='center',
                                position=position_dodge(1), dotsize = 0.5) + 
  geom_signif(comparisons = list(c("Lowlander", "Highlander")), 
              map_signif_level=TRUE)
dev.off()

# By treatment 
ggplot(ME_info, aes(x=treatment, y=ME37)) +
  geom_boxplot() + geom_dotplot(binaxis='y', stackdir='center',
                                position=position_dodge(1), dotsize = 0.5) +
  geom_signif(comparisons = list(c("Control", "Acclimated")), 
              map_signif_level=TRUE)
# Save as a pdf 
pdf(file=paste("results/rv/signed/module_expression_boxplots/treatment/eigengene_expr_ME25_",run_date,".pdf",sep=""),h=9,w=9)
ggplot(ME_info, aes(x=treatment, y=ME25)) +
  geom_boxplot() + geom_dotplot(binaxis='y', stackdir='center',
                                position=position_dodge(1), dotsize = 0.5) +
  geom_signif(comparisons = list(c("Control", "Acclimated")), 
              map_signif_level=TRUE)
dev.off()

#==============================================================================
# Line plots for treatment vs population effects 
ME4 <- ME_info[ , c("code","population","treatment","ME4")]

ME4 <- ddply(ME4, c("population", "treatment"), summarise,
              N    = sum(!is.na(ME4)),
              mean = mean(ME4, na.rm=TRUE),
              sd   = sd(ME4, na.rm=TRUE),
              se   = sd / sqrt(N)
)

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.1) # move them .05 to the left and right

# A much nicer plot
ggplot(ME4, aes(x=treatment, y=mean, colour=population, group=population)) + 
  scale_x_discrete(limits=c("Control","Acclimated")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(aes(shape = population), position=pd, fill="white", size=4) +
  xlab("Treatment") +
  ylab("mean expression") +
  scale_color_manual(values=c("black","grey")) +
  ggtitle("Module Expression for the right Ventricle") +
  expand_limits(y=0) +                        # Expand y range
  theme_bw() +
  theme(legend.justification=c(1,1),
        legend.position=c(1,1))               # Position legend in bottom right


pdf(file=paste("results/rv/signed/module_expression_line_plots/eigengene_expr_linePlot_ME10_",run_date,".pdf",sep=""),h=9,w=9)
ggplot(ME10, aes(x=treatment, y=mean, colour=population, group=population)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(aes(shape = population), position=pd, fill="white", size=4) +
  xlab("Treatment") +
  ylab("mean expression") +
  scale_color_manual(values=c("black","grey")) +
  ggtitle("ME10 Module Expression Right ventricle") +
  expand_limits(y=0) +                        # Expand y range
  theme_bw() +
  theme(legend.justification=c(1,1),
        legend.position=c(1,1))               # Position legend in bottom right
dev.off()

#==============================================================================
# Output scatterplots of module expression with phenotypes.
# We can view a plot here and make sure it looks like what we want for the output plots below.
ggplot(ME_info, aes(x=ME37, y=rv_standarized, shape=treatment, color=population)) + geom_point(size=5) + geom_smooth(method = lm, aes(group=1))  + theme_bw() + scale_color_manual(values=c("black","grey")) + theme(legend.position = "right",panel.grid = element_blank())

#PDF plot The first one with a legend
pdf(file=paste("results/rv/signed/module_expression_scatter_plots/phenotype_rv_standarized_ME37_",run_date,".pdf",sep=""),h=11, w=8.5)
ggplot(ME_info, aes(x=ME37, y=rv_standarized, shape=treatment, color=population)) + geom_point(size=5) + geom_smooth(method = lm, aes(group=1))  + theme_bw() + scale_color_manual(values=c("black","grey")) + theme(legend.position = "right",panel.grid = element_blank())
dev.off()

pdf(file=paste("results/rv/signed/module_expression_scatter_plots/phenotype_rv_standarized_ME33_",run_date,".pdf",sep=""),h=11, w=8.5)
ggplot(ME_info, aes(x=ME33, y=rv_standarized, shape=treatment, color=population)) + geom_point(size=5) + geom_smooth(method = lm, aes(group=1))  + theme_bw() + scale_color_manual(values=c("black","grey")) + theme(legend.position = "none",panel.grid = element_blank())
dev.off()

pdf(file=paste("results/rv/signed/module_expression_scatter_plots/phenotype_rv_standarized_ME3_",run_date,".pdf",sep=""),h=11, w=8.5)
ggplot(ME_info, aes(x=ME3, y=rv_standarized, shape=treatment, color=population)) + geom_point(size=5) + geom_smooth(method = lm, aes(group=1))  + theme_bw() + scale_color_manual(values=c("black","grey")) + theme(legend.position = "none",panel.grid = element_blank())
dev.off()

pdf(file=paste("results/rv/signed/module_expression_scatter_plots/phenotype_rv_standarized_ME5_",run_date,".pdf",sep=""),h=11, w=8.5)
ggplot(ME_info, aes(x=ME5, y=rv_standarized, shape=treatment, color=population)) + geom_point(size=5) + geom_smooth(method = lm, aes(group=1))  + theme_bw() + scale_color_manual(values=c("black","grey")) + theme(legend.position = "none",panel.grid = element_blank())
dev.off()

pdf(file=paste("results/rv/signed/module_expression_scatter_plots/phenotype_rv_standarized_ME25_",run_date,".pdf",sep=""),h=11, w=8.5)
ggplot(ME_info, aes(x=ME25, y=rv_standarized, shape=treatment, color=population)) + geom_point(size=5) + geom_smooth(method = lm, aes(group=1))  + theme_bw() + scale_color_manual(values=c("black","grey")) + theme(legend.position = "none",panel.grid = element_blank())
dev.off()

pdf(file=paste("results/rv/signed/module_expression_scatter_plots/phenotype_rv_standarized_ME32_",run_date,".pdf",sep=""),h=11, w=8.5)
ggplot(ME_info, aes(x=ME32, y=rv_standarized, shape=treatment, color=population)) + geom_point(size=5) + geom_smooth(method = lm, aes(group=1))  + theme_bw() + scale_color_manual(values=c("black","grey")) + theme(legend.position = "none",panel.grid = element_blank())
dev.off()

pdf(file=paste("results/rv/signed/module_expression_scatter_plots/phenotype_rv_standarized_ME8_",run_date,".pdf",sep=""),h=11, w=8.5)
ggplot(ME_info, aes(x=ME8, y=rv_standarized, shape=treatment, color=population)) + geom_point(size=5) + geom_smooth(method = lm, aes(group=1))  + theme_bw() + scale_color_manual(values=c("black","grey")) + theme(legend.position = "none",panel.grid = element_blank())
dev.off()

pdf(file=paste("results/rv/signed/module_expression_scatter_plots/phenotype_rv_standarized_ME4_",run_date,".pdf",sep=""),h=11, w=8.5)
ggplot(ME_info, aes(x=ME4, y=rv_standarized, shape=treatment, color=population)) + geom_point(size=5) + geom_smooth(method = lm, aes(group=1))  + theme_bw() + scale_color_manual(values=c("black","grey")) + theme(legend.position = "none",panel.grid = element_blank())
dev.off()

pdf(file=paste("results/rv/signed/module_expression_scatter_plots/phenotype_rv_standarized_ME15_",run_date,".pdf",sep=""),h=11, w=8.5)
ggplot(ME_info, aes(x=ME15, y=rv_standarized, shape=treatment, color=population)) + geom_point(size=5) + geom_smooth(method = lm, aes(group=1))  + theme_bw() + scale_color_manual(values=c("black","grey")) + theme(legend.position = "none",panel.grid = element_blank())
dev.off()

pdf(file=paste("results/rv/signed/module_expression_scatter_plots/phenotype_rv_standarized_ME7_",run_date,".pdf",sep=""),h=11, w=8.5)
ggplot(ME_info, aes(x=ME7, y=rv_standarized, shape=treatment, color=population)) + geom_point(size=5) + geom_smooth(method = lm, aes(group=1))  + theme_bw() + scale_color_manual(values=c("black","grey")) + theme(legend.position = "none",panel.grid = element_blank())
dev.off()

pdf(file=paste("results/rv/signed/module_expression_scatter_plots/phenotype_rv_standarized_ME10_",run_date,".pdf",sep=""),h=11, w=8.5)
ggplot(ME_info, aes(x=ME10, y=rv_standarized, shape=treatment, color=population)) + geom_point(size=5) + geom_smooth(method = lm, aes(group=1))  + theme_bw() + scale_color_manual(values=c("black","grey")) + theme(legend.position = "none",panel.grid = element_blank())
dev.off()

#==============================================================================
# GO enrichment WGCNA modules
# output gene info
all_genes <- geneInfo["Gene"]

# Run GO enrichment using gProfiler package
for (num in c("3","4","5","7","8","10","15","25","32","37")){
  genes <- geneInfo[geneInfo$moduleLabel==num,]["Gene"]
  genes_GO <- gost(as.vector(genes$Gene), organism = "pmbairdii",
                   ordered_query = F, significant = T, exclude_iea = F,
                   user_threshold=0.05,
                   correction_method = "g_SCS",
                   domain_scope = "annotated", numeric_ns = "",sources=c("GO","KEGG","REAC","TF","HP"))
  write.csv(apply(genes_GO$result,2,as.character),file=paste("results/rv/signed/GO_analysis/",num,"_module_GO_",run_date,".csv",sep=""))
}
# module 25 has no significant GO terms

#==============================================================================
# Is there a difference in Male vs Female rv_standardized weight?

id %>%
  group_by(sex) %>%
  get_summary_stats(rv_standardized_bm, type = "mean_sd")

# We can plot this as a box plot (does seem to be a difference)
bxp <- ggboxplot(
  id, x = "sex", y = "rv_standarized", 
  ylab = "rv_standarized", xlab = "sex", add = "jitter"
)
bxp

#Nicer plot  
ggplot(id, aes(x=sex, y=rv_standarized)) +
  geom_boxplot() + geom_dotplot(binaxis='y', stackdir='center',
                                position=position_dodge(1), dotsize = 0.5) +
  geom_signif(comparisons = list(c("M", "F")), 
              map_signif_level=TRUE)

# Any outliers? No.
id %>%
  group_by(sex) %>%
  identify_outliers(rv_standardized_bm)

# Compute Shapiro wilk test by groups
id %>%
  group_by(sex) %>%
  shapiro_test(rv_standardized_bm)
# Draw a qq plot by group
ggqqplot(id, x = "rv_standardized_bm", facet.by = "sex")

stat.test <- id %>% 
  t_test(rv_standardized_bm ~ sex) %>%
  add_significance()
stat.test

stat.test2 <- id %>%
  t_test(rv_standardized_bm ~ sex, var.equal = TRUE) %>%
  add_significance()
stat.test2
