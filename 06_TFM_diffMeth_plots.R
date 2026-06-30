####################
## Plots for DMRs ##
####################

library(ggvenn) # for Venn diagrams

library(ggplot2) # for plots
library(stringr) # for character manipulation
library(gridExtra) # for figure compositions

# Set Work Directory
setwd("D:/Ruben_backup/Postdoc/2022_AlonsoLabPostdoc/Analyses/EcicEcaz_2022/Ecic_WGBS_herb/ecic_herbivory_test/dmrs/dmrs")

#### Import DMR datasets ####

DMR.list <- list.files(recursive = TRUE, pattern = "metilene_C[pH][HG].[CH]_vs_[HS].txt")

DMR.tables <- list()
for (i in 1:length(DMR.list)) {
  DMR.tables[[i]] <- read.csv(DMR.list[[i]], sep = '\t', header = F)
  colnames(DMR.tables[[i]]) <- c("chr", "start", "stop", "q.value", "diffMeth", "cov", "p_MWU", "p_2D_KS", "mean_g1", "mean_g2")
  DMR.tables[[i]]$region <- paste0(DMR.tables[[i]]$chr, ":", DMR.tables[[i]]$start, "-", DMR.tables[[i]]$stop)
}

#### Venn diagrams ####

# Function for Venn diagrams with two sets (from WGBS analysis)
VenniVidiVinci2 <- function(metileneOut1, metileneOut2) {
  vennList <- list(dmrs1 = metileneOut1[which(metileneOut1$q.value < 0.05), ]$region, 
                   dmrs2 = metileneOut2[which(metileneOut2$q.value < 0.05), ]$region)
  ggvenn(vennList, 
         fill_color = c("purple", "yellow"), 
         stroke_size = 0.5, 
         set_name_size = 4, 
         auto_scale = T, 
         show_percentage = F)
}

# Function for Venn diagrams with three sets (from WGBS analysis)
VenniVidiVinci3 <- function(metileneOut1, metileneOut2, metileneOut3) {
  vennList <- list(dmrs1 = metileneOut1[which(metileneOut1$q.value < 0.05), ]$region, 
                   dmrs2 = metileneOut2[which(metileneOut2$q.value < 0.05), ]$region,
                   dmrs3 = metileneOut3[which(metileneOut3$q.value < 0.05), ]$region)
  ggvenn(vennList, 
         fill_color = c("orange", "purple", "yellow"), 
         stroke_size = 0.5, 
         set_name_size = 4, 
         show_percentage = F)
}

# Plot array
grid.arrange(
  VenniVidiVinci3(DMR.tables[[7]], DMR.tables[[8]], DMR.tables[[9]]), # CpG
  VenniVidiVinci3(DMR.tables[[1]], DMR.tables[[2]], DMR.tables[[3]]), # CHG
  VenniVidiVinci3(DMR.tables[[4]], DMR.tables[[5]], DMR.tables[[6]]), # CHH
  ncol = 3, nrow = 1)

#### Volcano Plots ####

# Function for editing group names
groupEdit <- function(group) {
  if (group == "C") {
    group <- "UN"
  } else if (group == "H") {
    group <- "MH"
  } else {
    group <- "HS"
  }
}

# Generate Volcano plot list
volcPlot <- list()
for (i in 1:length(DMR.list)) {
  context <- substr(DMR.list[[i]], 21, 23)
  groupA <- substr(DMR.list[[i]], 25, 25)
  groupA <- groupEdit(groupA)
  groupB <- substr(DMR.list[[i]], 30, 30)
  groupB <- groupEdit(groupB)
  volcData <- DMR.tables[[i]]
  volcData$methylation <- "NO"
  volcData$methylation[volcData$diffMeth > 0 & volcData$q.value < 0.05] <- "HIPO"
  volcData$methylation[volcData$diffMeth < 0 & volcData$q.value < 0.05] <- "HIPER"
  volcPlot[[i]] <- ggplot(data = volcData, 
                     aes(x = diffMeth, 
                         y = -log10(q.value), 
                         colour = methylation)) +	
    geom_point(alpha = 0.4, size = 1.75) + 
    theme_minimal() + labs(title = paste0(groupA, " vs. ", groupB, " (", context, ")")) +
    theme(legend.position = "none", 
          plot.title = element_text(hjust = 0.5), 
          axis.line = element_line(colour = "black", 
                                   size = 1, 
                                   linetype = "solid")) + 
    xlab("diffMeth") + ylab("-log10 FDR") +
    scale_colour_manual(values = 
                          if (length(table(volcData$methylation)) == 2) {
                            c("orange", "grey")
                          } else {
                            c("purple", "orange", "grey")
                          }
                        ) + 
    scale_x_continuous(limits = c(-0.66, 0.66)) + 
    scale_y_continuous(limits = c(0, 20))
}

# Plot array
grid.arrange(volcPlot[[7]], volcPlot[[8]], volcPlot[[9]], # CpG
             volcPlot[[1]], volcPlot[[2]], volcPlot[[3]], # CHG
             volcPlot[[4]], volcPlot[[5]], volcPlot[[6]], # CHH
             ncol = 3, nrow = 3)
