#1) Clear memory
{
  rm(list=ls())
}

#2) Load packages and set working directory
{
  library(parallel)
  library(tidyverse)
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  setwd("/media/felix/Elements/12_RangeExpansionGenomics")
}

#3) Read in data
{
  load("2_data/R_DataFiles/SlidingWindowAnalysis/SlidingWindow.20kb.AFC.RData")
}

#4) Find genes/snps that are common between at least 12 of the 16 populations
{
  #Summarize the data
  summarized.peaks <- df.genes.in.peaks %>% group_by(coordinate) %>% summarize(CHROM=unique(CHROM), POS=unique(POS), AFC_min=min(AFC_raw), AFC_max=max(AFC_raw),
                                                                         populations=n())
  
  #Add a variable to check if direction of change is the same
  summarized.peaks$const.change <- ifelse(summarized.peaks$AFC_min*summarized.peaks$AFC_max>0, T, F)
  
  #Filter variants where change of direction is the same and variants occur in at least 12 of the total populations
  summarized.peaks.filtered <- filter(summarized.peaks, const.change==T, populations >= 12)
  
  #Condition only true for 1 variant, so no use in continuing this analysis
}
