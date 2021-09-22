#1) Clear memory
{
  rm(list=ls())
}

#2) Load packages
{
  library(tidyverse)
  library(parallel)
  '%!in%' <- function(x,y)!('%in%'(x,y))
}

#3) Set working directory
{
  setwd("/media/felix/Elements/12_RangeExpansionGenomics")
}

#4) Read in data
{
  #4.1) Gene ontology data
  {
    #Load in data
    Gene.Ontology <- read_delim("2_data/TetThermRefGenome/uniprot-taxonomy_312017.tab", 
                                "\t", escape_double = FALSE, trim_ws = TRUE)
    
    #Subset data to get entries containing gene ID's (TTHERM_...)
    include <- unlist(lapply(Gene.Ontology$`Gene names`, grepl, pattern="TTHERM"))
    Gene.Ontology <- Gene.Ontology[include, ]
    
    #Extract Gene ID from gene names
    {
      #Process the problematic cases here (only 1 TTHERM present)
      data.problem <- filter(Gene.Ontology, !startsWith(`Gene names`, prefix="TTHERM"))
      
      #Separately process the one case where there is more than 1 TTHERM id
      special.case <- rbind(data.problem[5, ], data.problem[5, ])
      special.case$`Gene names` <- unlist(str_split(data.problem[5, "Gene names"], pattern = ";"))
      special.case$geneID <- unlist(str_split(special.case$`Gene names`, pattern = " "))[c(2,5)]
      
      #Process the other problematic data here
      data.problem <- data.problem[-5, ]
      data.problem$geneID <- ""
      for (i in 1:nrow(data.problem)){
        temp <- unlist(str_split(data.problem[i, "Gene names"], pattern=" "))
        ID <- temp[length(temp)]
        data.problem[i, "geneID"] <- ID
      }
      
      #Add geneID to regular data (subsetted to exclude special cases)
      Gene.Ontology <- filter(Gene.Ontology, startsWith(`Gene names`, prefix="TTHERM"))
      Gene.Ontology$geneID <- Gene.Ontology$`Gene names`
      
      #Bring all data back together
      Gene.Ontology.cleaned <- rbind(special.case, data.problem, Gene.Ontology)
    }
    
    #Make variable for numerical geneID 
    Gene.Ontology.cleaned$geneID_num <- as.numeric(sub(sub(Gene.Ontology.cleaned$geneID, pattern="TTHERM_", replacement=""), pattern="A", replacement=""))
    
    #Check entries with problems 
    temp1 <- filter(Gene.Ontology.cleaned, is.na(geneID_num))
    temp2 <- filter(Gene.Ontology.cleaned, is.na(geneID_num))
    temp1$`Gene names` <- unlist(str_split(temp1$`Gene names`, pattern = " "))[seq(from=1, to=36, by=2)]
    temp2$`Gene names` <- unlist(str_split(temp2$`Gene names`, pattern = " "))[seq(from=2, to=36, by=2)]
    temp.all <- rbind(temp1, temp2)
    temp.all$geneID <- temp.all$`Gene names`
    temp.all$geneID_num <- as.numeric(sub(sub(temp.all$geneID, pattern="TTHERM_", replacement=""), pattern="A", replacement=""))
    
    #Bring together with the rest of the data
    Gene.Ontology.cleaned <- rbind(filter(Gene.Ontology.cleaned, !is.na(geneID_num)), temp.all)
    
    #Check if all entries have a unique gene ID
    length(unique(Gene.Ontology.cleaned$geneID))
    
    #Remove one duplicated entry
    Gene.Ontology.cleaned <- Gene.Ontology.cleaned[!duplicated(Gene.Ontology.cleaned$geneID), ]
    
    #Check if numeric gene ontology is unique
    uniques <- Gene.Ontology.cleaned[duplicated(Gene.Ontology.cleaned$geneID_num), "geneID_num"]
    non_unique <- filter(Gene.Ontology.cleaned, geneID_num %in% uniques$geneID_num) %>% arrange(geneID_num)
    non_unique <- non_unique[seq(from=1, to=324, by=2), ]
    
    #Remove duplicated and add single entries to dataset
    Gene.Ontology.cleaned <- rbind(filter(Gene.Ontology.cleaned, geneID_num %!in% uniques$geneID_num), non_unique)
    
    #Save the cleaned gene ontology file
    save(Gene.Ontology.cleaned, file="2_data/R_DataFiles/Gene.Ontology.Cleaned.RData")

  }
}
