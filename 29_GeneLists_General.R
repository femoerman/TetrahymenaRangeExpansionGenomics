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

#3) Read and prepare data
{
  sites <- read_delim("2_data/GOWINDA_Data/AFC_General/SiteList.txt", 
                      "\t", escape_double = FALSE, col_names = FALSE, 
                      trim_ws = TRUE)
  sites$end <- sites$X2+1
  sites <- select(sites, X1, X2, end, X3, X4)
  
  #Save positions for assigning gene ID
  write_tsv((as.data.frame(sites)), path = "2_data/R_DataFiles/PoPoolationData/GENERAL_POS.txt", col_names = F)
  
  #Run intersect function of bedtools to identify gene ID's
  {
    #Define location of gtf file
    gtf.file <- paste(getwd(), "2_data/TetThermRefGenome/GCF_000189635.1_JCVI-TTA1-2.2_genomic_genesOnly.gtf", sep="/")
    
    #Define the location of the position file
    vcf.file <- paste(getwd(), "2_data/R_DataFiles/PoPoolationData/GENERAL_POS.txt", sep="/")
    
    #Define the output.file
    output.file <- paste(getwd(), "2_data/R_DataFiles/PoPoolationData/GENERAL_POS_geneID.txt", sep="/")
    
    #Perform system command to create annotation file
    command <- paste("bedtools intersect -wao -a ", vcf.file, " -b ", gtf.file, " > ", output.file, sep="")
    system(command)
  }
}

#5) Clear memory and read geneID data
{
  #Load in data
  rm(list=ls())
  peak_geneID <- read_delim("2_data/R_DataFiles/PoPoolationData/GENERAL_POS_geneID.txt", 
                            "\t", escape_double = FALSE, col_names = FALSE, 
                            trim_ws = TRUE)
  
  #Filter out locations without annotations
  data.filt <- filter(peak_geneID, X9 != -1)
  
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
  
  #Make the colnames more informative
  colnames(data.filt) <- c("#CHROM", "START", "STOP", "Populations", "coordinate", "CHROM", "geneID")
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

output <- data.filt %>% select(-START, -STOP, -CHROM, -geneID_num) %>% arrange(desc(Populations))
output2 <- select(output, geneID, geneFunction, Populations)
output2.uniquegene <- output2[!duplicated(output2$geneID), ]
output2 <- arrange(output2, geneID) %>% select(geneID, Populations, geneFunction)
print(output2.uniquegene, n=43)

#Get a list with proteins that have multiple variants
output3 <- output2 %>% group_by(geneID) %>% summarize(geneFunction=unique(geneFunction), count=n())
