#1) Clear memory
{
  rm(list=ls())
  setwd("/media/felix/Elements/12_RangeExpansionGenomics")
}

#2) Load packages
{
  library(tidyverse)
  library(parallel)
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  #setwd("/media/felix/Elements/12_RangeExpansionGenomics")
  setwd("/media/felix/DataDrive2/Documenten/PhD/12_RangeExpansionGenomics")
}

#3) Load CMH data and calculate genome coordinate value
{
  #Load data
  cmh.data <- read_delim("2_data/PopoolationData/Alldata_try2.cmh.gwas", 
                                 "\t", escape_double = FALSE, trim_ws = TRUE)
  
  #Calculate -log10 of p-value
  cmh.data$P_P <- -log10(cmh.data$P)
  
  #Calculate chromosome coordinate
  #Add position (coordinate)
  {
    #Read in scaffold data
    scaffolds <- read_delim("2_data/TetThermRefGenome/GCF_000189635.1_JCVI-TTA1-2.2_ScaffoldLength.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)
    scaffolds$name <- scaffolds$`Sequence-Name`
    scaffolds <- scaffolds[order(scaffolds$name),] 
    scaffolds$start <- 1
    for (i in 2:nrow(scaffolds)){
      scaffolds[i, "start"] <- scaffolds[i-1, "start"] + scaffolds[i-1, "Sequence-Length"]
    }
    scaffolds <- as.data.frame(scaffolds)
    rownames(scaffolds) <- scaffolds$`RefSeq-Accn`
    
    #Calculate genome coordinates
    cmh.data$coordinate <- scaffolds[cmh.data$CHR, "start"] + cmh.data$BP
    
    #Add a colour value for the chromosomes
    cmh.data$colour <- "yellow"
    for (i in 2:nrow(cmh.data)){
      cmh.data[i, "colour"] <- ifelse(
        cmh.data[i, "CHR"] == cmh.data[i-1, "CHR"],
        cmh.data[i-1, "colour"],
        ifelse(cmh.data[i-1, "colour"]=="yellow", "blue", "yellow")
      )
      
    }
  }
  
  #Save the whole dataset for future use
  save(cmh.data, file="2_data/R_DataFiles/PoPoolationData/RawDataset.RData")
  
}

#4) Visualize the results of the CMH test
{
  ggplot(cmh.data, aes(x=coordinate, y=P_P, colour=colour)) + geom_point(size=0.5, show.legend = F) + geom_hline(yintercept=8, colour="red", size=2) + 
    xlab("Genome coordinate") + ylab(expression(paste(-log[10], "of p-value"))) + 
    scale_colour_manual(breaks = c("yellow", "blue"), values=c("gold3", "blue4")) + 
    theme_light() + theme(axis.text=element_text(size=18), legend.text=element_text(size=24),legend.title=element_text(size=24), strip.text.x=element_text(30),
                          axis.title=element_text(size=24), strip.text = element_text(size=24)) 
  ggsave(file="CMH_data.pdf", device = "pdf", dpi = 320, units = "mm", width=575, height = 250)
}
