#1) Clear memory
{
  rm(list=ls())
  setwd("/media/felix/Elements/12_RangeExpansionGenomics")
}

#2) Load packages
{
  library(tidyverse)
  library(parallel)
}

#3) Create function to query vcf files (with variants only) and run in a parallellized way
{
  #Get a list with all vcf files
  vcf.files <-  paste(getwd(), list.files(path = "2_data/VariantCalling_OnlyDeviating/", 
                                                   recursive = T, pattern = "Mdup.vcf$", full.names = T), sep="/")
  
  #Create function to query files for necessary information
  query.vcf <- function(vcf.file, output.folder){
    
    #Get the name of the file
    temp.name <- unlist(str_split(vcf.file, pattern = "/"))[length(unlist(str_split(vcf.file, pattern = "/")))]
    name <- paste(unlist(str_split(temp.name, pattern=".bam"))[1], ".query", sep="")
  
    #Define the output file
    output.file <- paste(output.folder, name, sep="")
    
    #Issue query command
    command <- paste("bcftools query -H -f '%CHROM\t%POS\t%END\t%REF\t%ALT\t%QUAL\t%TYPE\t%DP\t[%AD]\n' ", vcf.file, " > ", output.file, sep="")
    system(command)
  }
  
  #Define output folder
  output.folder <- paste(getwd(), "2_data/QueriedVCF/VariantsOnly/", sep="/")
  
  #Perform parallelized query
 # mclapply(vcf.files, query.vcf, output.folder, mc.cores=6)
}

#4) Create function to query vcf files (AFC positions) and run in a parallellized way
{
  #Get a list with all vcf files
  vcf.files <-  paste(getwd(), list.files(path = "2_data/VariantCallingComplete/", 
                                          recursive = T, pattern = "Mdup.vcf$", full.names = T), sep="/")
  
  #Create function to query files for necessary information
  query.vcf <- function(vcf.file, output.folder){
    
    #Get the name of the file
    temp.name <- unlist(str_split(vcf.file, pattern = "/"))[length(unlist(str_split(vcf.file, pattern = "/")))]
    name <- paste(unlist(str_split(temp.name, pattern=".bam"))[1], ".query", sep="")
    
    #Compress and index vcf file
    system(paste("bgzip ", vcf.file))
    system(paste("tabix ", vcf.file, ".gz", sep=""))
    
    #Define the output file
    output.file <- paste(output.folder, name, sep="")
    
    #Refer to the positions file
    positions <- paste(getwd(), "2_data/TetThermRefGenome/Positions.txt", sep="/")
    
    #Issue query command
    command <- paste("bcftools query -H -R ",positions, " -f '%CHROM\t%POS\t%END\t%REF\t%ALT\t%QUAL\t%TYPE\t%DP\t[%AD]\n' ", vcf.file, ".gz > ", output.file, sep="")
    system(command)
  }
  
  #Define output folder
  # output.folder <- paste(getwd(), "2_data/QueriedVCF/AllPositions/", sep="/")
  
  #Perform parallelized query
  #mclapply(vcf.files, query.vcf, output.folder, mc.cores=6)
  
  #vcf.file <- vcf.files[1]
}

#5) Create annotation files (variants only)
{
  #Get a list with all vcf files
  vcf.files <-  paste(getwd(), list.files(path = "2_data/VariantCalling_OnlyDeviating/", 
                                          recursive = T, pattern = "Mdup.vcf$", full.names = T), sep="/")
  
  #Filter gtf file to only contain genes
  {
    #Read in
    gtf.source <- read_delim("2_data/TetThermRefGenome/GCF_000189635.1_JCVI-TTA1-2.2_genomic.gtf", 
                             "\t", escape_double = FALSE, col_names = FALSE, 
                             trim_ws = TRUE, skip = 3)
    
    #Filter for gene
    gtf.filter <- filter(gtf.source, X3=="gene")
    
    #Save output
    #write_tsv(gtf.filter,col_names = F, path="2_data/TetThermRefGenome/GCF_000189635.1_JCVI-TTA1-2.2_genomic_genesOnly.gtf")
  }
  
  
  
  #Create function to do intersection
  annotateVCF <- function(vcf.file, gtf.file, output.folder){
    
    #Get the name of the file
    temp.name <- unlist(str_split(vcf.file, pattern = "/"))[length(unlist(str_split(vcf.file, pattern = "/")))]
    name <- paste(unlist(str_split(temp.name, pattern=".bam"))[1], ".annot", sep="")
    
    #Define the output file
    output.file <- paste(output.folder, name, sep="/")
    
    #Perform system command to create annotation file
    command <- paste("bedtools intersect -wao -a ", vcf.file, " -b ", gtf.file, " > ", output.file, sep="")
    system(command)
    
  }
  
  #Define output folder
  output.folder <- paste(getwd(), "2_data/AnnotationsVCF/VariantsOnly", sep="/")
  
  #Define location of gtf file
  gtf.file <- paste(getwd(), "2_data/TetThermRefGenome/GCF_000189635.1_JCVI-TTA1-2.2_genomic_genesOnly.gtf", sep="/")
  
  #Perform parallellized annotation
  #mclapply(vcf.files, annotateVCF, gtf.file, output.folder, mc.cores=6)
}

#5) Create annotation files (AFC positions)
{
  #Get a list with all vcf files
  vcf.files <-  paste(getwd(), list.files(path = "2_data/VariantCallingComplete", 
                                          recursive = T, pattern = "Mdup.vcf.gz$", full.names = T), sep="/")
  
  #Subset vcf files to position
  {
    #Define output folder for subsetting
    output.folder <- paste(getwd(), "2_data/VariantCallingComplete_Subsetted/", sep="/")
    
    #Define file with positions to subset
    regions.file <- paste(getwd(), "2_data/TetThermRefGenome/Positions.txt", sep="/")
    
    #Make function for subsetting
    subsetVCF <- function(vcf.file, output.folder, regions.file){
      
      #Get the name of the file
      temp.name <- unlist(str_split(vcf.file, pattern = "/"))[length(unlist(str_split(vcf.file, pattern = "/")))]
      
      #Define output location
      output.file <- paste(output.folder, temp.name, sep="/")
      
      #Perform subsetting command
      system(paste("bcftools view -R ", regions.file, " -O z ", vcf.file, " -o ", output.file))
    }
    
    #Perform parallellized subsetting
    # mclapply(vcf.files, subsetVCF, output.folder, regions.file, mc.cores=6)
  }
  
  #Define files for intersection
  vcf.files <-  paste(getwd(), list.files(path = "2_data/VariantCallingComplete_Subsetted", 
                                          recursive = T, pattern = "Mdup.vcf.gz$", full.names = T), sep="/")
  
  #Create function to do intersection
  annotateVCF <- function(vcf.file, gtf.file, output.folder){
    
    #Get the name of the file
    temp.name <- unlist(str_split(vcf.file, pattern = "/"))[length(unlist(str_split(vcf.file, pattern = "/")))]
    name <- paste(unlist(str_split(temp.name, pattern=".bam"))[1], ".annot", sep="")
    
    #Define the output file
    output.file <- paste(output.folder, name, sep="/")
    
    #Perform system command to create annotation file
    system(paste("bedtools intersect -wao -a ", vcf.file, " -b ", gtf.file, " > ", output.file, sep=""))
  }
  
  #Define output folder
  output.folder <- paste(getwd(), "/2_data/AnnotationsVCF/AllPositions", sep="")
  
  #Define location of gtf file
  gtf.file <- paste(getwd(), "2_data/TetThermRefGenome/GCF_000189635.1_JCVI-TTA1-2.2_genomic_genesOnly.gtf", sep="/")
  
  #Perform parallellized annotation
  mclapply(vcf.files, annotateVCF, gtf.file, output.folder, mc.cores=6)
}