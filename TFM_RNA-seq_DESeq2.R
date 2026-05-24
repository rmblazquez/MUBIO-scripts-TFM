#### E. CICUTARIUM HERBIVORY ####

setwd("D:/Ruben_backup/Postdoc/2022_AlonsoLabPostdoc/Analyses/EcicEcaz_2022/Ecic_rnaseq_herb/results/deseq2/")

library(DESeq2)
counts <- read.csv("../featureCounts/read-counts-greater100.csv", sep = ";", header = T, row.names = 1)

library(stringr)
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

#### For herbivory and azacitidine general effects ####

expDesign <- model.matrix(~ herbivory + azacytidine + herbivory:azacytidine) 
ddsMatrix <- DESeqDataSetFromMatrix(countData = counts, colData = table, design = expDesign)

dds <- ddsMatrix[rowSums(counts(ddsMatrix)) > 1, ]
rld <- rlog(dds, blind = FALSE)
dds <- DESeq(ddsMatrix)
normCounts <- counts(dds, normalized = TRUE)
write.table(normCounts, file = "ecic-herb_normcounts_deseq2.txt", 
            sep = "\t", row.names = rownames(normCounts), col.names = colnames(normCounts), 
            quote = F, append = F)

par(mfcol = c(2,1))
boxplot(log(counts,2), col = as.factor(group)) # Raw counts
boxplot(log(normCounts, 2), col = as.factor(group)) # Norm counts

## data from control samples only
sampleTable <- data.frame(dds$sampleName, dds$condition, dds$azacytidine)
sampleTable.herbivory <- sampleTable[which(sampleTable$dds.azacytidine == "C"), ]
counts.herbivory <- counts[, sampleTable.herbivory$dds.sampleName]
normCounts.herbivory <- normCounts[, sampleTable.herbivory$dds.sampleName]

# boxplots
par(mfcol = c(2,1))
boxplot(log(counts.herbivory[,c(18,1:17)],2), col = as.factor(sampleTable.herbivory$dds.condition[c(18,1:17)])) # Raw counts
boxplot(log(normCounts.herbivory[,c(18,1:17)], 2), col = as.factor(sampleTable.herbivory$dds.condition[c(18,1:17)])) # Norm counts

# PCA (Only control treatments, no azacytidine)
sampleColor <- rep(c("white","green","red"), each = 6)
library(edgeR)
plotMDS(counts.herbivory[,c(18,1:17)], 
        gene.selection = "common", 
        top = 500, 
        dim.plot = c(1,2), 
        labels = NULL,
        # labels = colnames(counts.herbivory[,c(18,1:17)]),
        pch = 21,
        bg = sampleColor)

# Treatment glossary
# C_C = UN unharmed (different from C/control, which means no treatment with azacytidine)
# H_C = MH mechanical herbivory
# S_C = HS herbivory with saliva

# Test C_C Vs H_C
res.CCvsHC <- results(dds, contrast = c(0,-1,0,0,0,0))
data.CCvsHC <- data.frame(res.CCvsHC)
data.CCvsHC <- data.CCvsHC[with(data.CCvsHC, order(padj)), ]
write.table(res.CCvsHC, file = "13032024_deseq2_results_ecic-CCvsHC.txt", 
            sep = "\t", row.names = rownames(res.CCvsHC), col.names = colnames(res.CCvsHC), 
            quote = F, append = F)
dim(data.CCvsHC[which(data.CCvsHC$padj <= 0.05),]) # 3847 DEGs
dim(data.CCvsHC[which(data.CCvsHC$padj <= 0.05 & data.CCvsHC$log2FoldChange > 0),]) # 2013 Down-regulated herbivoryH
dim(data.CCvsHC[which(data.CCvsHC$padj <= 0.05 & data.CCvsHC$log2FoldChange < 0),]) # 1834 Up-regulated herbivoryH

# Test C_C Vs S_C
res.CCvsSC <- results(dds, contrast = c(0,0,-1,0,0,0))
data.CCvsSC <- data.frame(res.CCvsSC)
data.CCvsSC <- data.CCvsSC[with(data.CCvsSC, order(padj)), ]
write.table(res.CCvsSC, file = "13032024_deseq2_results_ecic-CCvsSC.txt", 
            sep = "\t", row.names = rownames(res.CCvsSC), col.names = colnames(res.CCvsSC), 
            quote = F, append = F)
dim(data.CCvsSC[which(data.CCvsSC$padj <= 0.05),]) # 1229 DEGs
dim(data.CCvsSC[which(data.CCvsSC$padj <= 0.05 & data.CCvsSC$log2FoldChange > 0),]) # 853 Down-regulated herbivoryS
dim(data.CCvsSC[which(data.CCvsSC$padj <= 0.05 & data.CCvsSC$log2FoldChange < 0),]) # 376 Up-regulated herbivoryS

# Test H_C Vs S_C
res.HCvsSC <- results(dds, contrast = c(0,1,-1,0,0,0))
data.HCvsSC <- data.frame(res.HCvsSC)
data.HCvsSC <- data.HCvsSC[with(data.HCvsSC, order(padj)), ]
write.table(res.HCvsSC, file = "13032024_deseq2_results_ecic-HCvsSC.txt", 
            sep = "\t", row.names = rownames(res.HCvsSC), col.names = colnames(res.HCvsSC), 
            quote = F, append = F)
dim(data.HCvsSC[which(data.HCvsSC$padj <= 0.05),]) # 3079 DEGs
dim(data.HCvsSC[which(data.HCvsSC$padj <= 0.05 & data.HCvsSC$log2FoldChange > 0),]) # 1554 Up-regulated herbivoryH
dim(data.HCvsSC[which(data.HCvsSC$padj <= 0.05 & data.HCvsSC$log2FoldChange < 0),]) # 1525 Up-regulated herbivoryS


#### PLOTS ####

# PCA for CCvsHC DEGs
DEGs.CCvsHC <- rownames(data.CCvsHC[which(data.CCvsHC$padj <= 0.05),])
DEGs.normCounts.CCvsHC <- normCounts[which(DEGs.CCvsHC %in% rownames(normCounts)), ]
# PCA for CCvsSC DEGs
DEGs.CCvsSC <- rownames(data.CCvsSC[which(data.CCvsSC$padj <= 0.05),])
DEGs.normCounts.CCvsSC <- normCounts[which(DEGs.CCvsSC %in% rownames(normCounts)), ]
# PCA for HCvsSC DEGs
DEGs.HCvsSC <- rownames(data.HCvsSC[which(data.HCvsSC$padj <= 0.05),])
DEGs.normCounts.HCvsSC <- normCounts[which(DEGs.HCvsSC %in% rownames(normCounts)), ]

par(mfcol = c(2,2))
plotMDS(DEGs.normCounts.CCvsHC, gene.selection = "common", top = 500, dim.plot = c(1,2), labels = NULL, 
        pch = 21,
        bg = sampleColor)
plotMDS(DEGs.normCounts.CCvsSC, gene.selection = "common", top = 500, dim.plot = c(1,2), labels = NULL, 
        pch = 21,
        bg = sampleColor)
plotMDS(DEGs.normCounts.HCvsSC, gene.selection = "common", top = 500, dim.plot = c(1,2), labels = NULL, 
        pch = 21,
        bg = sampleColor)

# PCA for all DEGs
DEGs.all <- unique(c(DEGs.CCvsHC, DEGs.CCvsSC, DEGs.HCvsSC))
DEGs.normCounts.all <- normCounts[which(DEGs.all %in% rownames(normCounts)), ]
plotMDS(DEGs.normCounts.all, 
        gene.selection = "common", 
        top = dim(DEGs.normCounts.all)[1], 
        dim.plot = c(1,2), 
        labels = NULL, 
        pch = 21,
        bg = sampleColor)


### Number of Up-regulated and Down-regulated DEGs with ggplot2 ###

library(reshape2) # Formatting tables
library(ggplot2) # Plotting graphs
library(gridExtra) # Arrange plots

# HerbivoryH
geneExpCH <- cbind(c(1834, 2013))
colnames(geneExpCH) <- "herbivoryH"
rownames(geneExpCH) <- c("Up-regulated", "Down-regulated")

geneExpCH.melt <- melt(geneExpCH)

plotCH.2 <- ggplot(data = geneExpCH.melt, aes(x = Var1, y = value, fill = Var1)) +  
  geom_bar(stat = "identity", position = position_dodge())  + theme_classic() +
  xlab("Herbivory") + ylab("# DEGs") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + theme(legend.position = "none")

# HerbivoryS
geneExpCS <- cbind(c(376, 853))
colnames(geneExpCS) <- "herbivoryS"
rownames(geneExpCS) <- c("Up-regulated", "Down-regulated")

geneExpCS.melt <- melt(geneExpCS)

plotCS.2 <- ggplot(data = geneExpCS.melt, aes(x = Var1, y = value, fill = Var1)) +  
  geom_bar(stat = "identity", position = position_dodge())  + theme_classic() +
  xlab("Saliva") + ylab("# DEGs") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + theme(legend.position = "none")

# Within Herbivory
geneExpHS <- cbind(c(1554, 1525))
colnames(geneExpHS) <- "withinHerbivory"
rownames(geneExpHS) <- c("HerbivoryH.Up", "HerbivoryS.Up")

geneExpHS.melt <- melt(geneExpHS)

plotHS.2 <- ggplot(data = geneExpHS.melt, aes(x = Var1, y = value, fill = Var1)) +  
  geom_bar(stat = "identity", position = position_dodge())  + theme_classic() +
  xlab("Herbivory H vs S") + ylab("# DEGs") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + theme(legend.position = "none")

grid.arrange(plotCH.2, plotCS.2, plotHS.2, ncol = 3, nrow = 1)


### HEATMAP ###
# DEGs overlapping the three contrasts

# Set DEG lists
DEGs.CCvsHC <- data.CCvsHC[which(data.CCvsHC$padj <= 0.05),]
DEGs.CCvsSC <- data.CCvsSC[which(data.CCvsSC$padj <= 0.05),]
DEGs.HCvsSC <- data.HCvsSC[which(data.HCvsSC$padj <= 0.05),]

# Set DEGs normalized counts
normCounts.CCvsHC <- normCounts[which(rownames(normCounts) %in% rownames(DEGs.CCvsHC)), ]
normCounts.CCvsSC <- normCounts[which(rownames(normCounts) %in% rownames(DEGs.CCvsSC)), ]
normCounts.HCvsSC <- normCounts[which(rownames(normCounts) %in% rownames(DEGs.HCvsSC)), ]

# Merge DEG normCounts and remove duplicated rows
degs <- rbind(normCounts.CCvsHC, normCounts.CCvsSC, normCounts.HCvsSC)
degs <- degs[!duplicated(degs), ]

library(RColorBrewer)
coul <- colorRampPalette(brewer.pal(8, "Blues"))(25)

samples.h <- colnames(degs)
genes.h <- rownames(degs)
group.h <- as.factor(herbivory)
levels(group.h) <- c("white","red","green")

library(gplots)
heatmap.2(as.matrix(scale(degs)), 
          trace = "none", 
          col = coul, 
          scale = "row", 
          labRow = FALSE, 
          colCol = as.vector(group.h))


### Volcano Plots ###

library(ggplot2)
library("RColorBrewer")

# HerbivoryH
volcDataH <- read.table('13032024_deseq2_results_ecic-CCvsHC.txt', sep = '\t', header = TRUE)
volcDataH$diffexpr <- "NO"
volcDataH$diffexpr[volcDataH$log2FoldChange > 0 & volcDataH$padj < 0.05] <- "DOWN"
volcDataH$diffexpr[volcDataH$log2FoldChange < 0 & volcDataH$padj < 0.05] <- "UP"
volcPlotH <- ggplot(data = volcDataH, aes(x = log2FoldChange, y = -log10(padj), colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "Control vs HerbivoryH") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 fold change") + ylab("-log10 q-value") +
  scale_colour_manual(values = c("steelblue", "grey", "red")) + 
  scale_x_continuous(limits = c(-25, 25)) + scale_y_continuous(limits = c(0, 25))

# HerbivoryS
volcDataS <- read.table('13032024_deseq2_results_ecic-CCvsSC.txt', sep = '\t', header = TRUE)
volcDataS$diffexpr <- "NO"
volcDataS$diffexpr[volcDataS$log2FoldChange > 0 & volcDataS$padj < 0.05] <- "DOWN"
volcDataS$diffexpr[volcDataS$log2FoldChange < 0 & volcDataS$padj < 0.05] <- "UP"
volcPlotS <- ggplot(data = volcDataS, aes(x = log2FoldChange, y = -log10(padj), colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "Control vs HerbivoryS") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 fold change") + ylab("-log10 q-value") +
  scale_colour_manual(values = c("steelblue", "grey", "red")) + 
  scale_x_continuous(limits = c(-25, 25)) + scale_y_continuous(limits = c(0, 25))

# Within Herbivory
volcDataHS <- read.table('13032024_deseq2_results_ecic-HCvsSC.txt', sep = '\t', header = TRUE)
volcDataHS$diffexpr <- "NO"
volcDataHS$diffexpr[volcDataHS$log2FoldChange > 0 & volcDataHS$padj < 0.05] <- "HERBIVORY"
volcDataHS$diffexpr[volcDataHS$log2FoldChange < 0 & volcDataHS$padj < 0.05] <- "SALIVA"
volcPlotHS <- ggplot(data = volcDataHS, aes(x = log2FoldChange, y = -log10(padj), colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "HerbivoryH vs HerbivoryS") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 fold change") + ylab("-log10 q-value") +
  scale_colour_manual(values = c("steelblue", "grey", "red")) + 
  scale_x_continuous(limits = c(-25, 25)) + scale_y_continuous(limits = c(0, 25))


library(gridExtra)
grid.arrange(volcPlotH, volcPlotS, volcPlotHS, ncol = 3, nrow = 1)

### MA PLOTS ###

MAPlotH <- ggplot(data = volcDataH, aes(x = log2(baseMean), y = log2FoldChange, colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "Control vs HerbivoryH") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 CPM") + ylab("log2 fold change") +
  scale_colour_manual(values = c("steelblue", "grey", "red"))

MAPlotS <- ggplot(data = volcDataS, aes(x = log2(baseMean), y = log2FoldChange, colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "Control vs HerbivoryS") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 CPM") + ylab("log2 fold change") +
  scale_colour_manual(values = c("steelblue", "grey", "red"))

MAPlotHS <- ggplot(data = volcDataHS, aes(x = log2(baseMean), y = log2FoldChange, colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "HerbivoryH vs HerbivoryS") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 CPM") + ylab("log2 fold change") +
  scale_colour_manual(values = c("steelblue", "grey", "red"))

library(gridExtra)
grid.arrange(MAPlotH, MAPlotS, MAPlotHS, ncol = 3, nrow = 1)

#### HYPERGEOMETRIC TEST ####

## Intersections
# DEGs between CCvsHC and CCvsSC: 331
length(intersect(as.vector(rownames(data.CCvsHC[which(data.CCvsHC$padj <= 0.05),])), as.vector(rownames(data.CCvsSC[which(data.CCvsSC$padj <= 0.05),]))))
# DEGs between CCvsHC and HCvsSC: 1084
length(intersect(as.vector(rownames(data.CCvsHC[which(data.CCvsHC$padj <= 0.05),])), as.vector(rownames(data.HCvsSC[which(data.HCvsSC$padj <= 0.05),]))))
# DEGs between CCvsSC and HCvsSC: 375
length(intersect(as.vector(rownames(data.CCvsSC[which(data.CCvsSC$padj <= 0.05),])), as.vector(rownames(data.HCvsSC[which(data.HCvsSC$padj <= 0.05),]))))

## Formal hypergeometric test
# CCvsHC U CCvsSC
phyper(330, 3847, 58670, 1229, lower.tail = FALSE, log.p = FALSE) # P-value < 0.001 (6.34611e-121) ***
# CCvsHC U HCvsSC
phyper(1083, 3847, 58670, 3079, lower.tail = FALSE, log.p = FALSE) # P-value < 0.001 (0) ***
# CCvsSC U HCvsSC
phyper(374, 1229, 61288, 3079, lower.tail = FALSE, log.p = FALSE) # P-value < 0.001 (9.556417e-191) ***

#q = size of overlap - 1
#m = number of DEGs in experiment comparison 1
#n = total number of genes (62517) - m
#k = number of DEGs in comparison experiment 2


#### PLOTS FOR GOIs ####

# Annotated DEGs in the Control vs Azacitidine comparison
annotationsDEG <- c("g25155.t1", "g22465.t1", "g81973.t1", "g64755.t1",	
                    "g169477.t1",	"g121748.t1",	"g15886.t1", "g141517.t1", 
                    "g148881.t1", "g61285.t1",	"g178124.t1",	"g60462.t1", 
                    "g88692.t1",	"g95990.t1")
degs.annot <- degs[, annotationsDEG]
degs.annot$Herbivory <- degs$Herbivory
degs.annot$Azacitidine <- degs$Azacitidine

library(reshape2)
degs.annot.melt <- melt(degs.annot)
library(ggplot2)
degs.annot.plot <- ggplot(degs.annot.melt, aes(x = Herbivory, y = log(value), color = Azacitidine)) + 
  #geom_line(aes(group = treatment)) + 
  geom_boxplot() + #scale_x_discrete(expand = c(0.2, 0.2)) + 
  facet_wrap(~ variable, scale = "free")

