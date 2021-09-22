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

#3) Load data Peaks dataset and data with assigned function 
{
  #Import data
  load("/media/felix/DataDrive2/Documenten/PhD/12_RangeExpansionGenomics/2_data/R_DataFiles/PoPoolationData/PeakData.RData")
  
  #Load data with gene function of variants
  load("2_data/R_DataFiles/AFC_Analysis/EvolvedData_AFC_func.RData")
}

#4) Assign geneID's to the PEAK data
{
  #Save the positions of the peaks data
  temp <- select(cmh.data.peaks, CHR, BP, P_P, distance, peak, coordinate)
  colnames(temp) <- c("#CHROM", "START", "P_P", "distance", "peak", "coordinate")
  temp$STOP <- temp$START+1
  temp <- select(temp, "#CHROM", START, STOP, P_P, distance, peak, coordinate)
  write_tsv((as.data.frame(temp)), path = "2_data/R_DataFiles/PoPoolationData/PEAKS_POS.txt", col_names = F)
  
  #Run intersect function of bedtools to identify gene ID's
  {
    #Define location of gtf file
    gtf.file <- paste(getwd(), "2_data/TetThermRefGenome/GCF_000189635.1_JCVI-TTA1-2.2_genomic_genesOnly.gtf", sep="/")
    
    #Define the location of the position file
    vcf.file <- paste(getwd(), "2_data/R_DataFiles/PoPoolationData/PEAKS_POS.txt", sep="/")
    
    #Define the output.file
    output.file <- paste(getwd(), "2_data/R_DataFiles/PoPoolationData/PEAKS_POS_geneID.txt", sep="/")
    
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
  peak_geneID <- read_delim("2_data/R_DataFiles/PoPoolationData/PEAKS_POS_geneID.txt", 
                            "\t", escape_double = FALSE, col_names = FALSE, 
                            trim_ws = TRUE)
  
  #Filter out locations without annotations
  data.filt <- filter(peak_geneID, X11 != -1)
  unique(data.filt$X6)
  unique(peak_geneID$X6)
  
  #Remove unnecessary columns
  data.filt <- select(data.filt, -X17, -X15, -X14, -X13, -X12, -X11, -X10, -X9)
  
  #Isolate geneID
  geneIDS <- str_split(data.filt$X16, patter='"')
  #Get the geneID
  fun1 <- function(lst, n){
    sapply(lst, `[`, n)
  }
  data.filt$geneID <- fun1(geneIDS, 4)
  data.filt$geneID <- replace_na(data.filt$geneID, ".")
  data.filt <- select(data.filt, -X16)
  
  #Make the colnames more informative
  colnames(data.filt) <- c("#CHROM", "START", "STOP", "P_P", "distance", "peak", "coordinate", "position", "geneID")
  
  #Make two subsetted datasets, 1) with all unique geneID's and 2) With the most significant geneID per peak (2 is non-informative)
  {
    data.all.genes <- data.filt %>% group_by(geneID) %>% summarize(maxSign=max(P_P)) %>% arrange(desc(maxSign))
    data.all.genes <- mutate(data.all.genes, geneID_num=as.numeric(substr(geneID, start=8, stop=nchar(geneID))))
    data.all.genes$coordinate <- 0
    for (i in 1:nrow(data.all.genes)){
      temp <- filter(data.filt, P_P==as.numeric(data.all.genes[i, "maxSign"]))
      data.all.genes[i, "coordinate"] <- as.numeric(temp$coordinate)
    }
  }
  
  #Assign the biological function to the genes
  {
    #Load gene ontology data
    load(file="2_data/R_DataFiles/Gene.Ontology.Cleaned.RData")
    Gene.Ontology.cleaned <- as.data.frame(Gene.Ontology.cleaned)
    
    #Change row names of Gene ontology file to the numeric geneID
    rownames(Gene.Ontology.cleaned) <- Gene.Ontology.cleaned$geneID_num
    
    #Assign function to the genes
    #Add gene function (protein name) and gene ontology to Biallellic Mutations data frame
    data.all.genes$geneFunction <- Gene.Ontology.cleaned[as.character(data.all.genes$geneID_num), "Protein names"]
  }
}

#6) Redo previous step, now getting for each peak the most significant gene/region
{
  #Organize per peak
  peaks <- peak_geneID %>% group_by(X6) %>% summarize(maxSign=max(X4))
  
  #Get the relevant information for each of the peaks
  {
    peak_sum <- filter(peak_geneID, X4==as.double(peaks[1, "maxSign"]), X6==as.double(peaks[1, "X6"]))
    
    for (i in 2:nrow(peaks)){
      temp <- filter(peak_geneID, X4==as.double(peaks[i, "maxSign"]), X6==as.double(peaks[i, "X6"]))
      peak_sum <- rbind(peak_sum, temp)
    }
  }
  
  #Select relevant variable
  peak_sum <- select(peak_sum, -X17, -X15, -X14, -X13, -X12, -X11, -X10, -X9)
  
  #Isolate geneID
  geneIDS <- str_split(peak_sum$X16, patter='"')
  #Get the geneID
  fun1 <- function(lst, n){
    sapply(lst, `[`, n)
  }
  peak_sum$geneID <- fun1(geneIDS, 4)
  peak_sum$geneID <- replace_na(peak_sum$geneID, ".")
  peak_sum <- select(peak_sum, -X16)
  
  #Make the colnames more informative
  colnames(peak_sum) <- c("#CHROM", "START", "STOP", "P_P", "distance", "peak", "coordinate", "position", "geneID")
  peak_sum <- mutate(peak_sum, geneID_num=as.numeric(substr(geneID, start=8, stop=nchar(geneID))))
  
  #Assign the biological function to the genes
  {
    #Load gene ontology data
    load(file="2_data/R_DataFiles/Gene.Ontology.Cleaned.RData")
    Gene.Ontology.cleaned <- as.data.frame(Gene.Ontology.cleaned)
    
    #Change row names of Gene ontology file to the numeric geneID
    rownames(Gene.Ontology.cleaned) <- Gene.Ontology.cleaned$geneID_num
    
    #Split in gene and non gene
    coding <- filter(peak_sum, geneID != ".")
    noncoding <- filter(peak_sum, geneID == ".")
    
    #Assign function to the genes
    #Add gene function (protein name) and gene ontology to Biallellic Mutations data frame
    coding$geneFunction <- Gene.Ontology.cleaned[as.character(coding$geneID_num), "Protein names"]
    noncoding$geneFunction <- "."
    
    #Bring back together 
    peak_sum <- rbind(coding, noncoding)
  }
}
peak_sum$geneID <- ifelse(peak_sum$geneID==".", "Non-coding region", peak_sum$geneID)
peak_sum$geneFunction <- ifelse(peak_sum$geneFunction==".", "Non-coding region", peak_sum$geneFunction)
data.both.peaks=list(genes=data.all.genes, peaks=peak_sum)
save(data.both.peaks, file="2_data/R_DataFiles/PoPoolationData/PeakSummaryGenes.RData")
