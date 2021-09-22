#1) Clear memory
{
  rm(list=ls())
}

#2) Load packages
{
  library(parallel)
  library(tidyverse)
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  
  setwd("/media/felix/Data/Felix_Data/Science/PhD/12_RangeExpansionGenomics")
}


#4) Load subsetted data containing the peaks
{
  #Load peak data
  load(file="2_data/R_DataFiles/PoPoolationData/PeakData.RData")
  
  #load significant gene data
  load(file="2_data/R_DataFiles/PoPoolationData/PeakSummaryGenes.RData")
  genes.sign <- data.both.peaks$genes
  peaks.sign <- data.both.peaks$peaks
  
  genes.sign$'Gene function' <- genes.sign$geneFunction
  peaks.sign$`Gene function` <- peaks.sign$geneFunction
  #4b) Subset data to identify peaks
  {
    #Plot data with the significant genes
    x <- ggplot(cmh.data.peaks, aes(x=coordinate, y=P_P, group=peak)) + geom_point(size=0.5, show.legend = F) + geom_hline(yintercept=8, colour="red", size=.5, alpha=0.5) + 
      geom_line() + geom_segment(inherit.aes=F, data=filter(genes.sign, maxSign>=8), mapping=aes(x = coordinate, y = maxSign + 4, xend = coordinate, 
                                                                                                 yend = maxSign + 1, colour=`Gene function`),
                   arrow = arrow(length = unit(0.3, "cm")), size=1.2) +
      xlab("Genome coordinate") + ylab(expression(paste(-log[10], " of p-value"))) +
      theme_light() + theme(axis.text=element_text(size=18), legend.text=element_text(size=12),legend.title=element_text(size=24), strip.text.x=element_text(30),
                            axis.title=element_text(size=16), strip.text = element_text(size=24), legend.position = "right")
    x
    data.x <- filter(genes.sign, maxSign>=8) %>% arrange(coordinate)
    ggsave(x, file="4_output/2_CMH_test_GO_Enrichment/CMH_data_Peaks.png", device = "png", dpi = 320, units = "mm", width=300, height = 150)
    
    #Plot data with most significant signal per peak
    y <- ggplot(cmh.data.peaks, aes(x=coordinate, y=P_P, group=peak)) + geom_point(size=0.5, show.legend = F) + geom_hline(yintercept=8, colour="red", size=0.5) + 
      geom_line() + geom_segment(inherit.aes=F, data=peaks.sign, mapping=aes(x = coordinate, y = P_P + 4, xend = coordinate, yend = P_P + 1, colour=`Gene function`),
                                 arrow = arrow(length = unit(0.3, "cm")))  + 
      xlab("Genome coordinate") + ylab(expression(paste(-log[10], " of p-value"))) +
      theme_light() + theme(axis.text=element_text(size=18), legend.text=element_text(size=12),legend.title=element_text(size=16), strip.text.x=element_text(30),
                            axis.title=element_text(size=24), strip.text = element_text(size=24), legend.position = "right")
    y
    ggsave(y, file="4_output/2_CMH_test_GO_Enrichment/CMH_data_Peaks_nc.png", device = "png", dpi = 320, units = "mm", width=300, height = 150)
    peaks.sign$significance <- peaks.sign$P_P
    data.y <- peaks.sign %>% arrange(coordinate) %>% select(-position, -geneID_num, -START, -STOP, -distance, -P_P)
    data.y
  }
  
  
  
}
