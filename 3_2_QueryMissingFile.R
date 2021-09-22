#1) Clear memory and load packages and set woking directory
{
  rm(list=ls())
  cores<-6
  
  library(parallel)
  library(tidyverse)
  
  setwd("/media/felix/Elements/12_RangeExpansionGenomics")
}

#Get a list with all vcf files
vcf.files <-  paste(getwd(), list.files(path = "2_data/VariantCallingComplete/", 
                                        recursive = T, pattern = "Mdup.vcf.gz$", full.names = T), sep="/")

vcf.file <- vcf.files[1]

#Define output folder
output.folder <- paste(getwd(), "2_data/QueriedVCF/AllPositions/", sep="/")

#Refer to the positions file
positions <- paste(getwd(), "2_data/TetThermRefGenome/Positions.txt", sep="/")

#Get the name of the file
temp.name <- unlist(str_split(vcf.file, pattern = "/"))[length(unlist(str_split(vcf.file, pattern = "/")))]
name <- paste(unlist(str_split(temp.name, pattern=".bam"))[1], ".query", sep="")

#Define the output file
output.file <- paste(output.folder, name, sep="")

#Issue query command for missing sample
command <- paste("bcftools query -H -R ",positions, " -f '%CHROM\t%POS\t%END\t%REF\t%ALT\t%QUAL\t%TYPE\t%DP\t[%AD]\n' ", vcf.file, " > ", output.file, sep="")
system(command)
