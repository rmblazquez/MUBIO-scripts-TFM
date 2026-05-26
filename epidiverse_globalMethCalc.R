#######################################
#### Global methylation calculator ####
#######################################

## Set path (in EBD Genomics-B)
setwd("/home/user/epidiv/ecic_herbivory_test/wgbs/wgbs/bedGraph")

## List with all files (bedGraph)

bedgraph <- list.files(recursive = TRUE, pattern = "*.bedGraph")
bedgraph.CpG <- bedgraph[37:54]
bedgraph.CHG <- bedgraph[1:18]
bedgraph.CHH <- bedgraph[19:36]

## Calculate percentage of methylated Cs per sample

# Generate empty  list and data frame
reports <- list()
globalMeth <- c()

# CpG
for (i in 1:length(bedgraph.CpG)) {
  reports[[i]] <- read.csv(bedgraph.CpG[i], sep = '\t', skip = 1, header = FALSE)
  # Filter positions with 1 or less reads
  reports[[i]] <- reports[[i]][(reports[[i]][,5] + reports[[i]][,6] > 1), ]
  # Calculate total percentage of methylated Cs
  globalMeth.sample <- c(bedgraph.CpG[i], sum(reports[[i]][,5]) / (sum(reports[[i]][,5]) + sum(reports[[i]][,6])))
  globalMeth <- rbind(globalMeth, globalMeth.sample)
  print(paste0("Sample ", bedgraph.CpG[i], " processed!"))  
}

gc()

# CHG
for (i in 1:length(bedgraph.CHG)) {
  reports[[i]] <- read.csv(bedgraph.CHG[i], sep = '\t', skip = 1, header = FALSE)
  reports[[i]] <- reports[[i]][(reports[[i]][,5] + reports[[i]][,6] > 1), ]
  globalMeth.sample <- c(bedgraph.CHG[i], sum(reports[[i]][,5]) / (sum(reports[[i]][,5]) + sum(reports[[i]][,6])))
  globalMeth <- rbind(globalMeth, globalMeth.sample)
  print(paste0("Sample ", bedgraph.CHG[i], " processed!"))  
}

gc()

# CHH
for (i in 1:length(bedgraph.CHH)) {
  reports[[i]] <- read.csv(bedgraph.CHH[i], sep = '\t', skip = 1, header = FALSE)
  reports[[i]] <- reports[[i]][(reports[[i]][,5] + reports[[i]][,6] > 1), ]
  globalMeth.sample <- c(bedgraph.CHH[i], sum(reports[[i]][,5]) / (sum(reports[[i]][,5]) + sum(reports[[i]][,6])))
  globalMeth <- rbind(globalMeth, globalMeth.sample)
  print(paste0("Sample ", bedgraph.CHH[i], " processed!"))
}

gc()

## Summarize methylation % table and edit metadata

# Edit column names and variable format
library(stringr)
colnames(globalMeth) <- c("sample", "globalMethylation")
globalMeth$globalMethylation <- round(as.numeric(globalMeth$globalMethylation) * 100, digits = 2)

# Add context and treatment variables
globalMeth$context <- as.factor(substr(globalMeth$sample, 15, 17))
globalMeth$context <- relevel(globalMeth$context, "CpG")
globalMeth$treatment <- as.factor(substr(globalMeth$sample, 36, 36))
levels(globalMeth$treatment) <- c("UN", "MH", "HS")

# Remove path and file extension
globalMeth$sample <- str_replace(globalMeth$sample, "C[pH][HG]/", "")
globalMeth$sample <- str_replace(globalMeth$sample, ".bedGraph", "")

write.table(globalMeth, file = "globalMeth_report.txt", sep = "\t", quote = FALSE)

## Plot the results

# Generate a faceted boxplot
library(ggplot2)
ggplot(data = globalMeth, aes(x = treatment, y = globalMethylation, fill = treatment)) + 
  geom_boxplot() + 
  facet_wrap(~ globalMeth$context, scales = "free") +
  theme_light()

## Statistical analysis

# Descriptive stats
for (i in 1:length(levels(globalMeth$context))) {
  print(
    paste0("Mean + SD Methylation percentage for ", 
            levels(globalMeth$context)[i], 
            " context is ", 
            round(as.numeric(mean(globalMeth$globalMethylation[which(globalMeth$context == levels(globalMeth$context)[i])])), digits = 2),
            " + ",
            round(as.numeric(sd(globalMeth$globalMethylation[which(globalMeth$context == levels(globalMeth$context)[i])])), digits = 2)
    )
  )
}

# Analysis of variance
methANOVA <- aov(globalMethylation ~ treatment + context + treatment * context, data = globalMeth)
summary(methANOVA)
methTukey <- TukeyHSD(methANOVA)

methTukey$`treatment:context`[c(1, 2, 9), ] # within CpG
methTukey$`treatment:context`[c(22, 23, 27), ] # within CHG
methTukey$`treatment:context`[34:36, ] # within CHH
