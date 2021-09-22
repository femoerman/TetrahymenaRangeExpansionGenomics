#1) Clear memory
{
  rm(list=ls())
}

#2) Load packages and set working directory
{
  library(parallel)
  library(tidyverse)
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  #setwd("/media/felix/Elements/12_RangeExpansionGenomics")
  setwd("/media/felix/DataDrive2/Documenten/PhD/12_RangeExpansionGenomics")
}

#3) Load data
{
  load("2_data/R_DataFiles/AFC_Analysis/EvolvedData_AFC_func.RData")
}

#4) Create a summary of the data for every site as row, and filter to keep sites with AFC in the same direction for all sites (with change being significant)
{
  #Get filtered data with only significant changes, and only sites that changed more than 0.3 (exclude clonal selection)
  df.sign <- filter(data.evolved.func, p_BinomProb>=10, AFC_abs>=0.3)
  
  #Create a summary by position of the variant
  summarized.sites <- df.sign %>% group_by(coordinate) %>% summarize(CHROM=unique(CHROM), POS=unique(POS), min.change=min(AFC_raw), 
                                                                     max.change=max(AFC_raw), populations=n())
  
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
  
  #Save this list as tab/delimited file without header for use in GOWINDA
  write_tsv(output.list, path="2_data/GOWINDA_Data/AFC_General/SiteList.txt", col_names = F)
  
  #Also save a list with all possible variant sites
  df.all <- select(summarized.sites, CHROM, POS, populations, coordinate)
  write_tsv(df.all, path="2_data/GOWINDA_Data/AFC_General/AllSites.txt", col_names = F)
}

#5) Perform GOWINDA analysis for this list
{
  #5.1) Define paths
  {
    #Define Gowinda path
    gowinda.path <- paste(getwd(), "5_Packages/Gowinda-1.12.jar", sep="/")
    
    #Define location of subsetted gwas file
    gwas.file <- paste(getwd(), "2_data/GOWINDA_Data/AFC_General/SiteList.txt", sep="/")
    
    #Define location of gene annotation file (.gtf) containing the entries for the exons
    gtf.path <- paste(getwd(),"2_data/TetThermRefGenome/Simplifiedgtf.gtf", sep="/")
    
    #Define the path to the file with the feneset data from GOminer
    geneset.path <-  paste(getwd(), "2_data/GOWINDA_Data/GOList_GOminer.txt", sep="/")
    
    #Define the path to the unfiltered GWAS data file
    unfiltered <- paste(getwd(), "2_data/GOWINDA_Data/AFC_General/AllSites.txt", sep="/")
  }
  
  #5.2) Perform basic analysis
  {
    #Define the output file
    output.file <- paste(getwd(), "2_data/GOWINDA_Data/AFC_General/BasicAnalysis.txt", sep="/")
    
    #Prepare the code to run
    code <- paste("java -Xmx4g -jar ", gowinda.path, " --snp-file ", unfiltered, " --candidate-snp-file ", gwas.file, " --gene-set-file ", geneset.path, " --annotation-file ", gtf.path, " --simulations 100000 --min-significance 1 --gene-definition gene --threads 6 --output-file ", output.file, " --mode gene --min-genes 1 ", sep="")
    
    #Run the code
    system(code)
  }
  
  #5.3) Perform analysis including regulatory regions
  {
    #Define the output file
    output.file <- paste(getwd(), "2_data/GOWINDA_Data/AFC_General/Analysis_with_Regulatory.txt", sep="/")
    
    #Prepare the code to run
    code <- paste("java -Xmx4g -jar ", gowinda.path, " --snp-file ", unfiltered, " --candidate-snp-file ", gwas.file, " --gene-set-file ", geneset.path, " --annotation-file ", gtf.path, " --simulations 100000 --min-significance 1 --gene-definition updownstream2000 --threads 6 --output-file ", output.file, " --mode gene --min-genes 1 ", sep="")
    
    #Run the code
    system(code)
  }
  
  #5.4) Perform high resolution GO enrichment analysis
  {
    #Define the output file
    output.file <- paste(getwd(), "2_data/GOWINDA_Data/AFC_General/Analysis_HighRes.txt", sep="/")
    
    #Prepare the code to run
    code <- paste("java -Xmx4g -jar ", gowinda.path, " --snp-file ", unfiltered, " --candidate-snp-file ", gwas.file, " --gene-set-file ", geneset.path, " --annotation-file ", gtf.path, " --simulations 100000 --min-significance 1 --gene-definition updownstream2000 --threads 6 --output-file ", output.file, " --mode snp --min-genes 1 ", sep="")
    
    #Run the code
    system(code)
  }
}
