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
  setwd("/media/felix/Elements/12_RangeExpansionGenomics")
}

#4) Load data on mutations, and filter to keep significant sites with AFC>0.3 only (exclude clonal selection)
{
  load("2_data/R_DataFiles/MutationAnalysis/AllMutations_Func.RData")
  Mutations.all.func.filt <- filter(Mutations.all.func, alt.freq>0.1)
}

#5) Subset data to contain gradient/uniform data and top10% data
{
  grad.all <- filter(Mutations.all.func, gradient=="Gradient") %>% select(CHROM, POS)
  unif.all <- filter(Mutations.all.func, gradient=="Uniform") %>% select(CHROM, POS)
  grad.10 <- filter(Mutations.all.func, gradient=="Gradient") %>% filter(alt.freq>= quantile(alt.freq, 0.9)) %>% select(CHROM, POS)
  unif.10 <- filter(Mutations.all.func, gradient=="Uniform") %>% filter(alt.freq>= quantile(alt.freq, 0.9)) %>% select(CHROM, POS)
  all.10 <- filter(Mutations.all.func) %>% filter(alt.freq>= quantile(alt.freq, 0.9)) %>% select(CHROM, POS)
}

#6) Save the lists as txt files for GOWINDA analysis
{
  write_tsv(grad.all, path = "2_data/GOWINDA_Data/Mutations/gradall.txt", col_names = F)
  write_tsv(unif.all, path = "2_data/GOWINDA_Data/Mutations/grada10.txt", col_names = F)
  write_tsv(grad.10, path = "2_data/GOWINDA_Data/Mutations/unifall.txt", col_names = F)
  write_tsv(unif.10, path = "2_data/GOWINDA_Data/Mutations/unif10.txt", col_names = F)
  write_tsv(all.10, path = "2_data/GOWINDA_Data/Mutations/all10.txt", col_names = F)
  
  #Save complete list
  complete.list <- select(Mutations.all.func, CHROM, POS)
  write_tsv(complete.list, path = "2_data/GOWINDA_Data/Mutations/completelist.txt", col_names = F)
}

#7) Perform GOWINDA analysis for mutations
{
  #7.0) Define locations
  {
    #Define Gowinda path
    gowinda.path <- paste(getwd(), "5_Packages/Gowinda-1.12.jar", sep="/")
    
    #Define location of gene annotation file (.gtf) containing the entries for the exons
    gtf.path <- paste(getwd(),"2_data/TetThermRefGenome/Simplifiedgtf.gtf", sep="/")
    
    #Define the path to the file with the geneset data from GOminer
    geneset.path <-  paste(getwd(), "2_data/GOWINDA_Data/GOList_GOminer.txt", sep="/")
    
    #Define the path to the unfiltered data file
    unfiltered <- paste(getwd(), "2_data/GOWINDA_Data/Mutations/completelist.txt", sep="/")
  }
  
  #7.1) Do first analysis (grad.all)
  {
    #Define input file
    gwas.file <- paste(getwd(), "2_data/GOWINDA_Data/Mutations/gradall.txt", sep="/")
    
    #Define the output file
    output.file <- paste(getwd(), "2_data/GOWINDA_Data/Mutations/GradAllFullAnalysis.txt", sep="/")
    
    #Prepare the code to run
    code <- paste("java -Xmx4g -jar ", gowinda.path, " --snp-file ", unfiltered, " --candidate-snp-file ", gwas.file, " --gene-set-file ", geneset.path, " --annotation-file ", gtf.path, " --simulations 100000 --min-significance 1 --gene-definition updownstream2000 --threads 6 --output-file ", output.file, " --mode snp --min-genes 1 ", sep="")
    
    #Run the code
    system(code)
  }
  
  #7.2) Do second analysis (unif.all)
  {
    #Define input file
    gwas.file <- paste(getwd(), "2_data/GOWINDA_Data/Mutations/unifall.txt", sep="/")
    
    #Define the output file
    output.file <- paste(getwd(), "2_data/GOWINDA_Data/Mutations/UnifAllFullAnalysis.txt", sep="/")
    
    #Prepare the code to run
    code <- paste("java -Xmx4g -jar ", gowinda.path, " --snp-file ", unfiltered, " --candidate-snp-file ", gwas.file, " --gene-set-file ", geneset.path, " --annotation-file ", gtf.path, " --simulations 100000 --min-significance 1 --gene-definition updownstream2000 --threads 6 --output-file ", output.file, " --mode snp --min-genes 1 ", sep="")
    
    #Run the code
    system(code)
  }
  
  #7.3) Do third analysis (grad.10)
  {
    #Define input file
    gwas.file <- paste(getwd(), "2_data/GOWINDA_Data/Mutations/grada10.txt", sep="/")
    
    #Define the output file
    output.file <- paste(getwd(), "2_data/GOWINDA_Data/Mutations/Grad10FullAnalysis.txt", sep="/")
    
    #Prepare the code to run
    code <- paste("java -Xmx4g -jar ", gowinda.path, " --snp-file ", unfiltered, " --candidate-snp-file ", gwas.file, " --gene-set-file ", geneset.path, " --annotation-file ", gtf.path, " --simulations 100000 --min-significance 1 --gene-definition updownstream2000 --threads 6 --output-file ", output.file, " --mode snp --min-genes 1 ", sep="")
    
    #Run the code
    system(code)
  }
  
  #7.4) Do fourth analysis (unif.10)
  {
    #Define input file
    gwas.file <- paste(getwd(), "2_data/GOWINDA_Data/Mutations/unif10.txt", sep="/")
    
    #Define the output file
    output.file <- paste(getwd(), "2_data/GOWINDA_Data/Mutations/Unif10FullAnalysis.txt", sep="/")
    
    #Prepare the code to run
    code <- paste("java -Xmx4g -jar ", gowinda.path, " --snp-file ", unfiltered, " --candidate-snp-file ", gwas.file, " --gene-set-file ", geneset.path, " --annotation-file ", gtf.path, " --simulations 100000 --min-significance 1 --gene-definition updownstream2000 --threads 6 --output-file ", output.file, " --mode snp --min-genes 1 ", sep="")
    
    #Run the code
    system(code)
  }
  
  #7.5) Do fifth analysis (all.10)
  {
    #Define input file
    gwas.file <- paste(getwd(), "2_data/GOWINDA_Data/Mutations/all10.txt", sep="/")
    
    #Define the output file
    output.file <- paste(getwd(), "2_data/GOWINDA_Data/Mutations/All10FullAnalysis.txt", sep="/")
    
    #Prepare the code to run
    code <- paste("java -Xmx4g -jar ", gowinda.path, " --snp-file ", unfiltered, " --candidate-snp-file ", gwas.file, " --gene-set-file ", geneset.path, " --annotation-file ", gtf.path, " --simulations 100000 --min-significance 1 --gene-definition updownstream2000 --threads 6 --output-file ", output.file, " --mode snp --min-genes 1 ", sep="")
    
    #Run the code
    system(code)
  }
}
