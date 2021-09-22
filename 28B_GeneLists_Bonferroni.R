#1) Clear memory
{
  rm(list=ls())
}

#2) Load packages
{
  library(parallel)
  library(tidyverse)
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  setwd("/media/felix/DataDrive2/Documenten/PhD/12_RangeExpansionGenomics")
}

#3) Load data bonferroni dataset and data with assigned function 
{
  #Import data
  load("/media/felix/DataDrive2/Documenten/PhD/12_RangeExpansionGenomics/2_data/R_DataFiles/PoPoolationData/BonferonniData.RData")
  
  #Load data with gene function of variants
  load("2_data/R_DataFiles/AFC_Analysis/EvolvedData_AFC_func.RData")
}

#4) Assign geneID's to the bonferroni data
{
  #Save the positions of the peaks data
  temp <- select(cmh.data.bonferroni, CHR, BP, P_P, coordinate)
  colnames(temp) <- c("#CHROM", "START", "P_P", "coordinate")
  temp$STOP <- temp$START+1
  temp <- select(temp, "#CHROM", START, STOP, P_P, coordinate)
  write_tsv((as.data.frame(temp)), path = "2_data/R_DataFiles/PoPoolationData/BONF_POS.txt", col_names = F)
  
  #Run intersect function of bedtools to identify gene ID's
  {
    #Define location of gtf file
    gtf.file <- paste(getwd(), "2_data/TetThermRefGenome/GCF_000189635.1_JCVI-TTA1-2.2_genomic_genesOnly.gtf", sep="/")
    
    #Define the location of the position file
    vcf.file <- paste(getwd(), "2_data/R_DataFiles/PoPoolationData/BONF_POS.txt", sep="/")
    
    #Define the output.file
    output.file <- paste(getwd(), "2_data/R_DataFiles/PoPoolationData/BONF_POS_geneID.txt", sep="/")
    
    #Perform system command to create annotation file
    command <- paste("bedtools intersect -wao -a ", vcf.file, " -b ", gtf.file, " > ", output.file, sep="")
    # "bedtools intersect -wao -a ", vcf.file, " -b ", gtf.file, " > ", output.file, sep="")
    system(command)
  }
  
}

#5) Clear memory and read geneID data
{
  #Load in data
  rm(list=ls())
  peak_geneID <- read_delim("2_data/R_DataFiles/PoPoolationData/BONF_POS_geneID.txt", 
                            "\t", escape_double = FALSE, col_names = FALSE, 
                            trim_ws = TRUE)
  
  #Filter out locations without annotations
  data.filt <- filter(peak_geneID, X9 != -1)
  unique(data.filt$X6)
  unique(peak_geneID$X6)
  
  #Remove unnecessary columns
  data.filt <- select(data.filt, -X15, -X13, -X12, -X11, -X10, -X9, -X8, -X7)
  
  #Isolate geneID
  geneIDS <- str_split(data.filt$X14, patter='"')
  #Get the geneID
  fun1 <- function(lst, n){
    sapply(lst, `[`, n)
  }
  data.filt$geneID <- fun1(geneIDS, 4)
  data.filt$geneID <- replace_na(data.filt$geneID, ".")
  data.filt <- select(data.filt, -X14)
  temp <- unlist(str_split(data.filt$geneID, "_"))
  data.filt$geneID_num <-  as.numeric(temp[temp !="TTHERM"])
  
  #Make the colnames more informative
  colnames(data.filt) <- c("#CHROM", "START", "STOP", "P_P", "coordinate", "position", "geneID", "geneID_num")
  data.filt <- mutate(data.filt, geneID_num=as.numeric(substr(geneID, start=8, stop=nchar(geneID))))
  
  #Assign the biological function to the genes
  {
    #Load gene ontology data
    load(file="2_data/R_DataFiles/Gene.Ontology.Cleaned.RData")
    Gene.Ontology.cleaned <- as.data.frame(Gene.Ontology.cleaned)
    
    #Change row names of Gene ontology file to the numeric geneID
    rownames(Gene.Ontology.cleaned) <- Gene.Ontology.cleaned$geneID_num
    
    #Assign function to the genes
    #Add gene function (protein name) and gene ontology to Biallellic Mutations data frame
    data.filt$geneFunction <- Gene.Ontology.cleaned[as.character(data.filt$geneID_num), "Protein names"]
  }
}

#6) Create a summary of the most important genes (i.e. those with most significant effects, and most common)
{
  #Get unique entries
  temp <- data.filt[!duplicated(data.filt$geneID), ]
  "NW_002475923.1" %in% temp
  common <- temp %>% group_by(geneFunction) %>% summarise(count=n()) %>% arrange(desc(count))
  
  
  highSign <- data.filt %>% arrange(desc(P_P))
}

#7) Create a nice list for Latex
{
  gene.list <- temp %>% select(geneID, geneFunction)
  gene.list$geneFunction <- paste(gene.list$geneFunction, " \\")
  write_delim(gene.list, delim = " & ", path="4_output/2_CMH_test_GO_Enrichment/GeneListBonferroni.txt", col_names = F)
}