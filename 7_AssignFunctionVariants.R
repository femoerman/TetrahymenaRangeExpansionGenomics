#1) Clear memory
{
  rm(list=ls())
}

#2) Load packages
{
  library(tidyverse)
  library(parallel)
}

#3) Set working directory
{
  setwd("/media/felix/Elements/12_RangeExpansionGenomics")
}

#4) Load and process data
{
  #Load cleaned gene ontology data
  load(file="2_data/R_DataFiles/Gene.Ontology.Cleaned.RData")
  Gene.Ontology.cleaned <- as.data.frame(Gene.Ontology.cleaned)
  
  #Change row names of Gene ontology file to the numeric geneID
  rownames(Gene.Ontology.cleaned) <- Gene.Ontology.cleaned$geneID_num
  
  #Load data with mutations (biallellic)
  load(file="2_data/R_DataFiles/MutationAnalysis/BiallellicMutations.RData")
  
  #Load data of all mutations
  load(file="2_data/R_DataFiles/MutationAnalysis/AllMutations.RData")
  
}

#5) Assign function to biallellic mutations
{
  #Separate Biallellic sites in coding/noncoding block
  coding <- filter(Mutations.biallellic, geneID_num != -1)
  noncoding <- filter(Mutations.biallellic, geneID_num == -1)
  
  #Add gene function (protein name) and gene ontology to Biallellic Mutations data frame
  coding$geneFunction <- Gene.Ontology.cleaned[as.character(coding$geneID_num), "Protein names"]
  coding$GO <- Gene.Ontology.cleaned[as.character(coding$geneID_num), "Gene ontology (GO)"]
  coding$GO_biol <- Gene.Ontology.cleaned[as.character(coding$geneID_num), "Gene ontology (biological process)"]
  coding$GO_mol <- Gene.Ontology.cleaned[as.character(coding$geneID_num), "Gene ontology (molecular function)"]
  coding$GO_cell <- Gene.Ontology.cleaned[as.character(coding$geneID_num), "Gene ontology (cellular component)"]
  
  #Add empty entries to noncoding block
  noncoding$geneFunction <- ""
  noncoding$GO <- ""
  noncoding$GO_biol <- ""
  noncoding$GO_mol <- ""
  noncoding$GO_cell <- ""
  
  #Bind together
  Mutations.biallellic.func <- rbind(coding, noncoding)
  
  #Save output
  save(Mutations.biallellic.func, file = "2_data/R_DataFiles/MutationAnalysis/BiallellicMutations_Func.RData")
}

#6) Assign function to all mutations
{
  #Separate Biallellic sites in coding/noncoding block
  coding <- filter(Mutations.all, geneID_num != -1)
  noncoding <- filter(Mutations.all, geneID_num == -1)
  
  #Add gene function (protein name) and gene ontology to Biallellic Mutations data frame
  coding$geneFunction <- Gene.Ontology.cleaned[as.character(coding$geneID_num), "Protein names"]
  coding$GO <- Gene.Ontology.cleaned[as.character(coding$geneID_num), "Gene ontology (GO)"]
  coding$GO_biol <- Gene.Ontology.cleaned[as.character(coding$geneID_num), "Gene ontology (biological process)"]
  coding$GO_mol <- Gene.Ontology.cleaned[as.character(coding$geneID_num), "Gene ontology (molecular function)"]
  coding$GO_cell <- Gene.Ontology.cleaned[as.character(coding$geneID_num), "Gene ontology (cellular component)"]
  
  #Add empty entries to noncoding block
  noncoding$geneFunction <- ""
  noncoding$GO <- ""
  noncoding$GO_biol <- ""
  noncoding$GO_mol <- ""
  noncoding$GO_cell <- ""
  
  #Bind together
  Mutations.all.func <- rbind(coding, noncoding)
  
  #Save output
  save(Mutations.all.func, file = "2_data/R_DataFiles/MutationAnalysis/AllMutations_Func.RData")
}

#7) Assign function to the AFC data
{
  #Load the data for the AFC
  load("2_data/R_DataFiles/AFC_Analysis/EvolvedData_AFC.RData")
  
  #Separate Biallellic sites in coding/noncoding block
  coding <- filter(data.evolved, !is.na(geneID_num))
  noncoding <- filter(data.evolved, is.na(geneID_num))
  
  #Add gene function (protein name) and gene ontology to Biallellic Mutations data frame
  coding$geneFunction <- Gene.Ontology.cleaned[as.character(coding$geneID_num), "Protein names"]
  coding$GO <- Gene.Ontology.cleaned[as.character(coding$geneID_num), "Gene ontology (GO)"]
  coding$GO_biol <- Gene.Ontology.cleaned[as.character(coding$geneID_num), "Gene ontology (biological process)"]
  coding$GO_mol <- Gene.Ontology.cleaned[as.character(coding$geneID_num), "Gene ontology (molecular function)"]
  coding$GO_cell <- Gene.Ontology.cleaned[as.character(coding$geneID_num), "Gene ontology (cellular component)"]
  
  #Add empty entries to noncoding block
  noncoding$geneFunction <- ""
  noncoding$GO <- ""
  noncoding$GO_biol <- ""
  noncoding$GO_mol <- ""
  noncoding$GO_cell <- ""
  
  #Bind together
  data.evolved.func <- rbind(coding, noncoding)
  
  #Save output
  save(data.evolved.func, file = "2_data/R_DataFiles/AFC_Analysis/EvolvedData_AFC_func.RData")
}