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
  setwd("/media/felix/DataDrive2/Documenten/PhD/12_RangeExpansionGenomics")
}

#4) Load data on mutations, and filter to keep significant sites with AFC>0.3 only (exclude clonal selection)
{
  load("2_data/R_DataFiles/MutationAnalysis/AllMutations_Func.RData")
  Mutations.all.func.filt <- filter(Mutations.all.func, alt.freq>0.1)
  
  #Calculate the percentage of snp versus indel
  Mutations.all.func %>% group_by(TYPE) %>% summarize(count=n())
}

#5) Check how common mutations are in all genes (i.e. in how many populations they occur)
{
  gradient <- c("EVO4",  "EVO5",  "EVO8",  "EVO9", "EVO11", "EVO14", "EVO18", "EVO20")
  
  Mutations.all.freq <- Mutations.all.func %>% group_by(geneID) %>% 
    summarise(populations = length(unique(sample)), samples = list(unique(sample)), geneFunction=unique(geneFunction)) %>% 
    arrange(desc(populations))
  Mutations.all.freq <- data.frame(Mutations.all.freq)
  Mutations.all.freq$gradientCount <- 0
  for (i in 1:nrow(Mutations.all.freq)){
    Mutations.all.freq[i, "gradientCount"] <- ifelse(length(unlist(Mutations.all.freq[i, "samples"]))==0, 0, 
     length(intersect(gradient, unlist(Mutations.all.freq[i, "samples"]))))
  }
  Mutations.all.freq$uniformCount <- Mutations.all.freq$populations -  Mutations.all.freq$gradientCount
  Mutations.all.freq$difference <- abs(Mutations.all.freq$gradientCount - Mutations.all.freq$uniformCount)
  Mutations.all.freq <- select(Mutations.all.freq, -samples)
}

#6) Make histogram of frequencies for both groups, and create list with top candidates for general and gradient specific adaptation
{
  hist(Mutations.all.freq$populations)
  hist(Mutations.all.freq$gradientCount - Mutations.all.freq$uniformCount)
  
  #Get top candidates for general adaptation
  general <- Mutations.all.freq %>% filter(populations >= 12)
  general <- filter(general, geneID != ".")
  print(select(general, -gradientCount, -uniformCount, -difference), max.print = 66)
  sum.gen <- select(general, -gradientCount, -uniformCount, -difference)
  
  #Differently represented in the two treatments
  abiotic <- Mutations.all.freq %>% arrange(desc(difference)) %>% filter(difference>=4)
  print(select(abiotic, -populations), max.print = 66)
  sum.grad <- select(abiotic, -populations)
}