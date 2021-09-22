#1) Clear memory
{
  rm(list=ls())
}

#2) Load packages
{
  library(parallel)
  library(tidyverse)
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  setwd("/media/felix/Elements/12_RangeExpansionGenomics")
}

#3) Load data
{
  #Load the cmh data
  load("2_data/R_DataFiles/PoPoolationData/RawDataset.RData")
  
  #Order the data by coordinate
  cmh.data <- arrange(cmh.data, coordinate)
  
  #Plot and save raw data
  #Plot filtered data to verify
  y <- ggplot(cmh.data, aes(x=coordinate, y=P_P)) + geom_point(size=0.5, show.legend = F) + geom_hline(yintercept=8, colour="red", size=2) + 
    xlab("Genome coordinate") + ylab(expression(paste(-log[10], "of p-value"))) + 
    theme_light() + theme(axis.text=element_text(size=18), legend.text=element_text(size=24),legend.title=element_text(size=24), strip.text.x=element_text(30),
                          axis.title=element_text(size=24), strip.text = element_text(size=24)) 
  ggsave(y, file="4_output/2_CMH_test_GO_Enrichment/CMH_data_raw.png", device = "png", dpi = 320, units = "mm", width=575, height = 250)
}

#4) Create subset data for a) all values significant after bonferonni and b) containing the peaks
{
  #4a) Subset bonferonni data
  {
    #This means, data where p-value < 0.001 / 100.000 or P_P >=8
    cmh.data.bonferroni <- filter(cmh.data, P_P>=8)
  }
  
  #4b) Subset data to identify peaks
  {
    
    #Add a variable that shows the distance of the variant site to the previous one
    cmh.data$distance <- 0
    cmh.data[2:97690, "distance"] <- cmh.data[2:97690, "coordinate"] - cmh.data[1:97689, "coordinate"]
    
    #Show a histogram of the distances
    hist(log(cmh.data$distance), breaks = 200)
    
    #The expected average distance between 2 snps is 1054.4 bp
    #One set of close variants (peak) is defined as all variants with a maximum distance to the next variant of half this number (527 bp)
    
    cmh.data$peak <- 1
    for ( i in 2:nrow(cmh.data)){
      cmh.data[i, "peak"] <- ifelse(cmh.data[i, "distance"]-cmh.data[i-1, "distance"]<=54, cmh.data[i-1, "peak"], cmh.data[i-1, "peak"]+1)
    }
    
    #A real peak is defined as a peak which has at least one variant with a significant p-value (P_P>=8) 
    # and at least a length of X variants
    
    #Add peak length and max P_P to the data 
    datasum <- cmh.data %>% group_by(peak) %>% summarize(length=n(), maxp=max(P_P), largeP=sum(P_P>=8))
    datasum$percSign <- datasum$largeP/datasum$length
    hist(datasum$length, breaks = 200)
    
    #Filter data containing only peaks with less than 10% significant observations
    filteredsum <- filter(datasum, percSign>=0.05)
    cmh.data.peaks <- filter(cmh.data, peak %in% filteredsum$peak)
    
    #Filter out peaks with fewer than 20 observations per peak
    filteredsum <- filter(filteredsum, length>=20)
    cmh.data.peaks <-filter(cmh.data, peak %in% filteredsum$peak)
    
    #Plot filtered data to verify
    x <- ggplot(cmh.data.peaks, aes(x=coordinate, y=P_P)) + geom_point(size=0.5, show.legend = F) + geom_hline(yintercept=8, colour="red", size=2) + 
      xlab("Genome coordinate") + ylab(expression(paste(-log[10], " of p-value"))) +
      theme_light() + theme(axis.text=element_text(size=18), legend.text=element_text(size=24),legend.title=element_text(size=24), strip.text.x=element_text(30),
                            axis.title=element_text(size=24), strip.text = element_text(size=24))
    ggsave(x, file="4_output/2_CMH_test_GO_Enrichment/CMH_data_Peaks.png", device = "png", dpi = 320, units = "mm", width=575, height = 250)
  }
  
  #Save the files to later use with GOWINDA
  save(cmh.data.peaks, file="2_data/R_DataFiles/PoPoolationData/PeakData.RData")
  write_tsv(cmh.data.peaks, path="2_data/R_DataFiles/PoPoolationData/PeakData.gwas", col_names = F)
  
  #Save the highly significant changes
  save(cmh.data.bonferroni, file="2_data/R_DataFiles/PoPoolationData/BonferonniData.RData")
  write_tsv(cmh.data.bonferroni, path="2_data/R_DataFiles/PoPoolationData/BonferonniData.gwas", col_names=F)
  
  #Create a combined plot
  library(ggpubr)
  x2 <- ggarrange(y, x,
                 labels = c("A", "B"),
                 ncol = 1, nrow = 2)
  ggsave(x2, file="4_output/2_CMH_test_GO_Enrichment/CMH_data_All.png", device = "png", dpi = 320, units = "mm", width=285, height = 250)
  
}
