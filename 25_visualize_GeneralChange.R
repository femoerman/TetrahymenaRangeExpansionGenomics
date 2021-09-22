#1) Clear memory
{
  rm(list=ls())
}

#2) Load packages
{
  library(parallel)
  library(tidyverse)
  library(nlme)
  library(MuMIn)
  '%!in%' <- function(x,y)!('%in%'(x,y))
}

#3) Set working directory
{
  setwd("/media/felix/Elements/12_RangeExpansionGenomics")
  load("2_data/R_DataFiles/AFC_Analysis/EvolvedData_AFC_func.RData")
}

#4) Create a summary of the data for every site as row, and filter to keep sites with AFC in the same direction for all sites (with change being significant)
{
  #Get filtered data with only significant changes, and only sites that changed more than 0.3 (exclude clonal selection)
  df.sign <- filter(data.evolved.func, p_BinomProb>=10, AFC_abs>=0.3)
  
  #Create a summary by position of the variant
  summarized.sites <- df.sign %>% group_by(coordinate) %>% summarize(CHROM=unique(CHROM), POS=unique(POS), min.change=min(AFC_raw), 
                                                                     max.change=max(AFC_raw), populations=n(), averagep=mean(p_BinomProb))
  
  #Add variable that tells if change happened in same direction
  summarized.sites$const.change <- ifelse(summarized.sites$min.change*summarized.sites$max.change>0, T, F)
  
  #Filter the summarized dataframe with the folowing conditions
  {
    #1: We observe change in the same direction for all variants
    summarized.sites.filtered <- filter(summarized.sites, const.change==T)
    
    #2: We observe change in at least 75% of the populations
    summarized.sites.filtered <- filter(summarized.sites.filtered, populations>=12)
  }
  
  #Output a list with the positions where we see change happening
  output.list <- select(summarized.sites.filtered, CHROM, POS, populations, coordinate)
  output.list.all <- select(summarized.sites, CHROM, POS, populations, coordinate, averagep)
}

#5) Add a variable to shoz if change is general, and visualize data
{
  #Add variable
  output.list.all$general <- ifelse(output.list.all$coordinate %in% output.list$coordinate, "General", "Not general")
  
  #Plot output
  x <- ggplot(output.list.all, aes(x=coordinate, y=populations, colour=general)) + geom_point()+ 
    xlab("Genome coordinate") + ylab("Number of populations with significant AFC-change") + 
    theme_light() + theme(axis.text=element_text(size=12), legend.text=element_text(size=12),legend.title=element_text(size=16), strip.text.x=element_text(24),
                          axis.title=element_text(size=16), strip.text = element_text(size=16),
                          axis.text.x=element_text(angle=90, hjust=1)) + 
    scale_color_manual(values=c("#d7191c", "#2c7bb6"))
  x
  ggsave(x, path = "4_output/2_CMH_test_GO_Enrichment/", filename="Counts_general.png", device="png", height=170, width=250, dpi=320, units = "mm")
}