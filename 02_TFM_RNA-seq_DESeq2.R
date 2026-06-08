#### E. CICUTARIUM HERBIVORY ####

## Setting work environment

setwd("D:/Ruben_backup/Postdoc/2022_AlonsoLabPostdoc/Analyses/EcicEcaz_2022/Ecic_rnaseq_herb/results/deseq2/")

library(stringr) # extract patterns in characters

library(DESeq2) # perform differential gene expression analysis
library(edgeR) # generate PCA/MDS

library(reshape2) # Formatting tables
library(ggplot2) # Plotting graphs
library(gridExtra) # Arrange plots

library(RColorBrewer) # generate color palettes
library(gplots) # use of heatmap.2()
library(ggvenn) # for Venn diagrams

## Define functions

# Function for Venn diagrams with two sets

VenniVidiVinci2 <- function(deseqOut1, deseqOut2) {
  vennList <- list(degs1 = rownames(deseqOut1[which(deseqOut1$padj < 0.05), ]), 
                   degs2 = rownames(deseqOut2[which(deseqOut2$padj < 0.05), ]))
  ggvenn(vennList, 
         fill_color = c("lightblue", "tomato"), 
         stroke_size = 0.5, 
         set_name_size = 4, 
         auto_scale = T, 
         show_percentage = F)
}

# Function for Venn diagrams with three sets

VenniVidiVinci3 <- function(deseqOut1, deseqOut2, deseqOut3) {
  vennList <- list(degs1 = rownames(deseqOut1[which(deseqOut1$padj < 0.05), ]), 
                   degs2 = rownames(deseqOut2[which(deseqOut2$padj < 0.05), ]),
                   degs3 = rownames(deseqOut3[which(deseqOut3$padj < 0.05), ]))
  ggvenn(vennList, 
         fill_color = c("lightgreen", "lightblue", "tomato"), 
         stroke_size = 0.5, 
         set_name_size = 4, 
         show_percentage = F)
}

# Function to calculate hypergeometric tests

autoPhyper <- function(deseqOut1, deseqOut2) {
  deglist1 <- rownames(deseqOut1[which(deseqOut1$padj <= 0.05),])
  group1 <- length(deglist1) # number of DEGs in group 1
  deglist2 <- rownames(deseqOut2[which(deseqOut2$padj <= 0.05),])
  group2 <- length(deglist2) # number of DEGs in group 2
  overlap <- length(intersect(deglist1, deglist2)) # number of overlapping DEGs between groups
  totalgenes <- length(unique(c(rownames(deseqOut1), rownames(deseqOut2)))) # all analyzed genes
  HGtest <- phyper(overlap, group2, totalgenes - group2, group1, lower.tail = FALSE, log.p = FALSE) # Hyper geometric test
  # The order of the values in the phyper() function is q, m, n, k:
  #q = size of overlap - 1
  #m = number of DEGs in experiment comparison 1
  #n = total number of genes (62517) - m
  #k = number of DEGs in comparison experiment 2
  print(paste("Hypergeometric test p-value: ", HGtest)) # print the p-value
}


## Prepare input data

# Load counts
counts <- read.csv("../featureCounts/read-counts-greater100.csv", sep = ";", header = T, row.names = 1)

# Generate experimental design table 
samples <- colnames(counts)
herbivory <- as.factor(str_sub(samples, -3, -3))
azacytidine <- as.factor(str_sub(samples, -1, -1))
individual <- as.factor(str_sub(samples, 5, 6))
group <- as.factor(str_sub(samples, -3, -1))
herbivory <- relevel(herbivory, "C")
azacytidine <- relevel(azacytidine, "C")
group <- relevel(group, "C_C")

table <- data.frame(sampleName = samples, fileName = samples, condition = herbivory)
table$azacytidine <- azacytidine

# The experiment involved both exposure to azacytidine and herbivory treatments
# Thus, the counts are being fit to a model with both variables.
# In this study, we are focusing on the herbivory treatments, so the azacytidine 
# contrasts will be ignored despite being included in the model.

# Fitting the counts to a model with variables herbivory and azacytidine
expDesign <- model.matrix(~ herbivory + azacytidine + herbivory:azacytidine) 
ddsMatrix <- DESeqDataSetFromMatrix(countData = counts, colData = table, design = expDesign)

# Perform differential gene expression analysis
dds <- ddsMatrix[rowSums(counts(ddsMatrix)) > 1, ]
rld <- rlog(dds, blind = FALSE)
dds <- DESeq(ddsMatrix)
normCounts <- counts(dds, normalized = TRUE)
# In case we need to export the normalized counts
# write.table(normCounts, file = "ecic-herb_normcounts_deseq2.txt",
#             sep = "\t", row.names = rownames(normCounts), col.names = colnames(normCounts), 
#             quote = F, append = F)

# Boxplot comparing counts before and after normalization (all treatments) ()
# par(mfcol = c(2,1))
# boxplot(log(counts,2), col = as.factor(group)) # Raw counts
# boxplot(log(normCounts, 2), col = as.factor(group)) # Norm counts

## data from "control" samples (i.e., no azacytidine+ samples) only
sampleTable <- data.frame(dds$sampleName, dds$condition, dds$azacytidine)
sampleTable.herbivory <- sampleTable[which(sampleTable$dds.azacytidine == "C"), ]
counts.herbivory <- counts[, sampleTable.herbivory$dds.sampleName]
normCounts.herbivory <- normCounts[, sampleTable.herbivory$dds.sampleName]
normCounts.herbivory <- normCounts.herbivory[,c(18,1:17)]

# Boxplot comparing counts before and after normalization (only samples without azacytidine treatment)
par(mfcol = c(1,2))
boxplot(log(counts.herbivory[,c(18,1:17)],2), col = as.factor(sampleTable.herbivory$dds.condition[c(18,1:17)])) # Raw counts
boxplot(log(normCounts.herbivory, 2), col = as.factor(sampleTable.herbivory$dds.condition[c(18,1:17)])) # Norm counts

# PCA (Only samples with no azacytidine)
sampleColor <- rep(c("white","green","red"), each = 6)
plotMDS(counts.herbivory[,c(18,1:17)], 
        gene.selection = "common", 
        top = 500, 
        dim.plot = c(1,2), 
        labels = NULL,
        pch = 21,
        bg = sampleColor)


## Differential gene expression comparisons 

# Treatment glossary
# C_C = UN unharmed (different from C/control, which means no treatment with azacytidine)
# H_C = MH mechanical herbivory
# S_C = HS herbivory with saliva

# Using the contrast vector: the index represents these groups of samples in object expDesign:
# 1 = (Intercept) 
# 2 = herbivoryH 
# 3 = herbivoryS 
# 4 = azacytidineA 
# 5 = herbivoryH:azacytidineA 
# 6 = herbivoryS:azacytidineA
# The reference group is herviboryC. Thus, if we want to:
# compare UN vs MH (herbivoryC vs herbivoryH): contrast = c(0,-1,0,0,0,0) # logFC < 0: Up-regulated herbivoryH
# compare UN vs HS (herbivoryC vs herbivoryS): contrast = c(0,0,-1,0,0,0) # logFC < 0: Up-regulated herbivoryS
# compare MH vs HS (herbivoryH vs herbivoryS): contrast = c(0,1,-1,0,0,0) # logFC < 0: Up-regulated herbivoryS


# Test UN (C_C) vs MH (H_C)
res.CCvsHC <- results(dds, contrast = c(0,-1,0,0,0,0))
data.CCvsHC <- data.frame(res.CCvsHC)
data.CCvsHC <- data.CCvsHC[with(data.CCvsHC, order(padj)), ]
# write.table(res.CCvsHC, file = "03062026_deseq2_results_ecic-CCvsHC.txt", 
#             sep = "\t", row.names = rownames(res.CCvsHC), col.names = colnames(res.CCvsHC), 
#             quote = F, append = F)
dim(data.CCvsHC[which(data.CCvsHC$padj <= 0.05),]) # 3847 DEGs
dim(data.CCvsHC[which(data.CCvsHC$padj <= 0.05 & data.CCvsHC$log2FoldChange > 0),]) # 2013 Down-regulated herbivoryH
dim(data.CCvsHC[which(data.CCvsHC$padj <= 0.05 & data.CCvsHC$log2FoldChange < 0),]) # 1834 Up-regulated herbivoryH

# Test UN (C_C) Vs HS (S_C)
res.CCvsSC <- results(dds, contrast = c(0,0,-1,0,0,0))
data.CCvsSC <- data.frame(res.CCvsSC)
data.CCvsSC <- data.CCvsSC[with(data.CCvsSC, order(padj)), ]
# write.table(res.CCvsSC, file = "03062026_deseq2_results_ecic-CCvsSC.txt", 
#             sep = "\t", row.names = rownames(res.CCvsSC), col.names = colnames(res.CCvsSC), 
#             quote = F, append = F)
dim(data.CCvsSC[which(data.CCvsSC$padj <= 0.05),]) # 1229 DEGs
dim(data.CCvsSC[which(data.CCvsSC$padj <= 0.05 & data.CCvsSC$log2FoldChange > 0),]) # 853 Down-regulated herbivoryS
dim(data.CCvsSC[which(data.CCvsSC$padj <= 0.05 & data.CCvsSC$log2FoldChange < 0),]) # 376 Up-regulated herbivoryS

# Test MH (H_C) Vs HS (S_C)
res.HCvsSC <- results(dds, contrast = c(0,1,-1,0,0,0))
data.HCvsSC <- data.frame(res.HCvsSC)
data.HCvsSC <- data.HCvsSC[with(data.HCvsSC, order(padj)), ]
# write.table(res.HCvsSC, file = "03062026_deseq2_results_ecic-HCvsSC.txt", 
#             sep = "\t", row.names = rownames(res.HCvsSC), col.names = colnames(res.HCvsSC), 
#             quote = F, append = F)
dim(data.HCvsSC[which(data.HCvsSC$padj <= 0.05),]) # 3079 DEGs
dim(data.HCvsSC[which(data.HCvsSC$padj <= 0.05 & data.HCvsSC$log2FoldChange > 0),]) # 1554 Up-regulated herbivoryH
dim(data.HCvsSC[which(data.HCvsSC$padj <= 0.05 & data.HCvsSC$log2FoldChange < 0),]) # 1525 Up-regulated herbivoryS


#### PLOTS ####

# # PCA for CCvsHC DEGs
# DEGs.CCvsHC <- rownames(data.CCvsHC[which(data.CCvsHC$padj <= 0.05),])
# DEGs.normCounts.CCvsHC <- normCounts.herbivory[which(DEGs.CCvsHC %in% rownames(normCounts.herbivory)), ]
# # PCA for CCvsSC DEGs
# DEGs.CCvsSC <- rownames(data.CCvsSC[which(data.CCvsSC$padj <= 0.05),])
# DEGs.normCounts.CCvsSC <- normCounts.herbivory[which(DEGs.CCvsSC %in% rownames(normCounts.herbivory)), ]
# # PCA for HCvsSC DEGs
# DEGs.HCvsSC <- rownames(data.HCvsSC[which(data.HCvsSC$padj <= 0.05),])
# DEGs.normCounts.HCvsSC <- normCounts.herbivory[which(DEGs.HCvsSC %in% rownames(normCounts.herbivory)), ]
# # PCA for all DEGs
# DEGs.all <- unique(c(DEGs.CCvsHC, DEGs.CCvsSC, DEGs.HCvsSC))
# DEGs.normCounts.all <- normCounts.herbivory[which(DEGs.all %in% rownames(normCounts.herbivory)), ]
# 
# par(mfcol = c(2,2))
# plotMDS(DEGs.normCounts.CCvsHC, gene.selection = "common", top = 500, dim.plot = c(1,2), labels = NULL, 
#         pch = 21,
#         bg = sampleColor)
# plotMDS(DEGs.normCounts.CCvsSC, gene.selection = "common", top = 500, dim.plot = c(1,2), labels = NULL, 
#         pch = 21,
#         bg = sampleColor)
# plotMDS(DEGs.normCounts.HCvsSC, gene.selection = "common", top = 500, dim.plot = c(1,2), labels = NULL, 
#         pch = 21,
#         bg = sampleColor)
# plotMDS(DEGs.normCounts.all, 
#         gene.selection = "common", 
#         top = dim(DEGs.normCounts.all)[1], 
#         dim.plot = c(1,2), 
#         labels = NULL, 
#         pch = 21,
#         bg = sampleColor)


## Number of Up-regulated and Down-regulated DEGs with ggplot2

# UN vs. MH
geneExpCH <- cbind(c(
  dim(data.CCvsHC[which(data.CCvsHC$padj <= 0.05 & data.CCvsHC$log2FoldChange < 0),])[1],
  dim(data.CCvsHC[which(data.CCvsHC$padj <= 0.05 & data.CCvsHC$log2FoldChange > 0),])[1]
))
colnames(geneExpCH) <- "UN vs. MH"
rownames(geneExpCH) <- c("Up-regulated", "Down-regulated")

geneExpCH.melt <- melt(geneExpCH)

plotCH.2 <- ggplot(data = geneExpCH.melt, aes(x = Var1, y = value, fill = Var1)) +  
  geom_bar(stat = "identity", position = position_dodge())  + theme_classic() +
  xlab("UN vs. MH") + ylab("# DEGs") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2100)) + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

# UN vs. HS
geneExpCS <- cbind(c(
  dim(data.CCvsSC[which(data.CCvsSC$padj <= 0.05 & data.CCvsSC$log2FoldChange < 0),])[1],
  dim(data.CCvsSC[which(data.CCvsSC$padj <= 0.05 & data.CCvsSC$log2FoldChange > 0),])[1]
))
colnames(geneExpCS) <- "UN vs. HS"
rownames(geneExpCS) <- c("Up-regulated", "Down-regulated")

geneExpCS.melt <- melt(geneExpCS)

plotCS.2 <- ggplot(data = geneExpCS.melt, aes(x = Var1, y = value, fill = Var1)) +  
  geom_bar(stat = "identity", position = position_dodge())  + theme_classic() +
  xlab("UN vs. HS") + ylab("# DEGs") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2100)) + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

# MH vs. HS
geneExpHS <- cbind(c(
  dim(data.HCvsSC[which(data.HCvsSC$padj <= 0.05 & data.HCvsSC$log2FoldChange > 0),])[1], 
  dim(data.HCvsSC[which(data.HCvsSC$padj <= 0.05 & data.HCvsSC$log2FoldChange < 0),])[1]
))
colnames(geneExpHS) <- "MH vs. HS"
rownames(geneExpHS) <- c("Up-reg. in MH", "Up-reg. in HS")

geneExpHS.melt <- melt(geneExpHS)

plotHS.2 <- ggplot(data = geneExpHS.melt, aes(x = Var1, y = value, fill = Var1)) +  
  geom_bar(stat = "identity", position = position_dodge())  + theme_classic() +
  xlab("MH vs. HS") + ylab("# DEGs") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2100)) + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

# Arrange the barplots
grid.arrange(plotCH.2, plotCS.2, plotHS.2, ncol = 3, nrow = 1)


## HEATMAP
# DEGs overlapping the three contrasts

# Set DEG lists (again)
DEGs.CCvsHC <- data.CCvsHC[which(data.CCvsHC$padj <= 0.05),]
DEGs.CCvsSC <- data.CCvsSC[which(data.CCvsSC$padj <= 0.05),]
DEGs.HCvsSC <- data.HCvsSC[which(data.HCvsSC$padj <= 0.05),]

# Set DEGs normalized counts
normCounts.CCvsHC <- normCounts.herbivory[which(rownames(normCounts.herbivory) %in% rownames(DEGs.CCvsHC)), ]
normCounts.CCvsSC <- normCounts.herbivory[which(rownames(normCounts.herbivory) %in% rownames(DEGs.CCvsSC)), ]
normCounts.HCvsSC <- normCounts.herbivory[which(rownames(normCounts.herbivory) %in% rownames(DEGs.HCvsSC)), ]

# Merge DEG normCounts and remove duplicated rows
degs <- rbind(normCounts.CCvsHC, normCounts.CCvsSC, normCounts.HCvsSC)
degs <- degs[!duplicated(degs), ]

coul <- colorRampPalette(brewer.pal(8, "Blues"))(25)

samples.h <- colnames(degs)
genes.h <- rownames(degs)
sampleColor.2 <- rep(c("black","green","red"), each = 6)
group.h <- as.factor(sampleColor.2)

heatmap.2(as.matrix(log2(degs + 1)), 
          trace = "none", 
          col = coul, 
          # scale = "row", 
          labRow = FALSE, 
          colCol = as.vector(group.h))


### Volcano Plots ###

# In case datasets are not loaded
# volcDataH <- read.table('03062024_deseq2_results_ecic-CCvsHC.txt', sep = '\t', header = TRUE)
# volcDataS <- read.table('03062024_deseq2_results_ecic-CCvsSC.txt', sep = '\t', header = TRUE)
# volcDataHS <- read.table('03062024_deseq2_results_ecic-HCvsSC.txt', sep = '\t', header = TRUE)

# HerbivoryH
volcDataH <- res.CCvsHC
volcDataH$diffexpr <- "NO"
volcDataH$diffexpr[volcDataH$log2FoldChange > 0 & volcDataH$padj < 0.05] <- "DOWN"
volcDataH$diffexpr[volcDataH$log2FoldChange < 0 & volcDataH$padj < 0.05] <- "UP"
volcPlotH <- ggplot(data = volcDataH, aes(x = log2FoldChange, y = -log10(padj), colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "UN vs. MH") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 FC") + ylab("-log10 FDR") +
  scale_colour_manual(values = c("steelblue", "grey", "red")) + 
  scale_x_continuous(limits = c(-25, 25)) + scale_y_continuous(limits = c(0, 25))

# HerbivoryS
volcDataS <- res.CCvsSC
volcDataS$diffexpr <- "NO"
volcDataS$diffexpr[volcDataS$log2FoldChange > 0 & volcDataS$padj < 0.05] <- "DOWN"
volcDataS$diffexpr[volcDataS$log2FoldChange < 0 & volcDataS$padj < 0.05] <- "UP"
volcPlotS <- ggplot(data = volcDataS, aes(x = log2FoldChange, y = -log10(padj), colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "UN vs. HS") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 FC") + ylab("-log10 FDR") +
  scale_colour_manual(values = c("steelblue", "grey", "red")) + 
  scale_x_continuous(limits = c(-25, 25)) + scale_y_continuous(limits = c(0, 25))

# Within Herbivory
volcDataHS <- res.HCvsSC
volcDataHS$diffexpr <- "NO"
volcDataHS$diffexpr[volcDataHS$log2FoldChange > 0 & volcDataHS$padj < 0.05] <- "MH"
volcDataHS$diffexpr[volcDataHS$log2FoldChange < 0 & volcDataHS$padj < 0.05] <- "HS"
volcPlotHS <- ggplot(data = volcDataHS, aes(x = log2FoldChange, y = -log10(padj), colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "MH vs. HS") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 FC") + ylab("-log10 FDR") +
  scale_colour_manual(values = c("steelblue", "grey", "red")) + 
  scale_x_continuous(limits = c(-25, 25)) + scale_y_continuous(limits = c(0, 25))


grid.arrange(volcPlotH, volcPlotS, volcPlotHS, ncol = 3, nrow = 1)

### MA PLOTS ###

MAPlotH <- ggplot(data = volcDataH, aes(x = log2(baseMean), y = log2FoldChange, colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "UN vs. MH") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 CPM") + ylab("log2 FC") +
  scale_colour_manual(values = c("steelblue", "grey", "red"))

MAPlotS <- ggplot(data = volcDataS, aes(x = log2(baseMean), y = log2FoldChange, colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "UN vs. HS") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 CPM") + ylab("log2 FC") +
  scale_colour_manual(values = c("steelblue", "grey", "red"))

MAPlotHS <- ggplot(data = volcDataHS, aes(x = log2(baseMean), y = log2FoldChange, colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "MH vs. HS") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 CPM") + ylab("log2 FC") +
  scale_colour_manual(values = c("red", "steelblue", "grey")) # had to rearrange colors due to different labelling

grid.arrange(MAPlotH, MAPlotS, MAPlotHS, ncol = 3, nrow = 1)

#### OVERLAPPING DEGS ####

## Intersections
# DEGs between CCvsHC and CCvsSC: 331
length(intersect(as.vector(rownames(data.CCvsHC[which(data.CCvsHC$padj <= 0.05),])), as.vector(rownames(data.CCvsSC[which(data.CCvsSC$padj <= 0.05),]))))
# DEGs between CCvsHC and HCvsSC: 1084
length(intersect(as.vector(rownames(data.CCvsHC[which(data.CCvsHC$padj <= 0.05),])), as.vector(rownames(data.HCvsSC[which(data.HCvsSC$padj <= 0.05),]))))
# DEGs between CCvsSC and HCvsSC: 375
length(intersect(as.vector(rownames(data.CCvsSC[which(data.CCvsSC$padj <= 0.05),])), as.vector(rownames(data.HCvsSC[which(data.HCvsSC$padj <= 0.05),]))))

## Venn diagrams:

# All DEGs
venn.all <-
VenniVidiVinci3(data.CCvsHC, data.CCvsSC, data.HCvsSC)

# Up-regulated by herbivory MH and HS
venn.Up <-
VenniVidiVinci2(data.CCvsHC[which(data.CCvsHC$log2FoldChange < 0),],
                data.CCvsSC[which(data.CCvsSC$log2FoldChange < 0),])

# Down-regulated by herbivory MH and HS
venn.Down <-
VenniVidiVinci2(data.CCvsHC[which(data.CCvsHC$log2FoldChange > 0),],
                data.CCvsSC[which(data.CCvsSC$log2FoldChange > 0),])

grid.arrange(venn.all, venn.Up, venn.Down, ncol = 3)

## Hypergeometric tests

# UN vs. MH U UN vs. HS
autoPhyper(data.CCvsHC, data.CCvsSC)
#[1] "Hypergeometric test p-value:  1.04374200927257e-121"

# UN vs. MH U MH vs. HS
autoPhyper(data.CCvsHC, data.HCvsSC)
#[1] "Hypergeometric test p-value:  0"

# UN vs. HS U MH vs. HS
autoPhyper(data.CCvsSC, data.HCvsSC)
#[1] "Hypergeometric test p-value:  1.00131749907196e-191"



# #### PLOTS FOR GOIs ####
# 
# # Annotated DEGs in the Control vs Azacitidine comparison
# annotationsDEG <- c("g25155.t1", "g22465.t1", "g81973.t1", "g64755.t1",	
#                     "g169477.t1",	"g121748.t1",	"g15886.t1", "g141517.t1", 
#                     "g148881.t1", "g61285.t1",	"g178124.t1",	"g60462.t1", 
#                     "g88692.t1",	"g95990.t1")
# degs.annot <- degs[, annotationsDEG]
# degs.annot$Herbivory <- degs$Herbivory
# degs.annot$Azacitidine <- degs$Azacitidine
# 
# degs.annot.melt <- melt(degs.annot)
# degs.annot.plot <- ggplot(degs.annot.melt, aes(x = Herbivory, y = log(value), color = Azacitidine)) + 
#   #geom_line(aes(group = treatment)) + 
#   geom_boxplot() + #scale_x_discrete(expand = c(0.2, 0.2)) + 
#   facet_wrap(~ variable, scale = "free")
# 
