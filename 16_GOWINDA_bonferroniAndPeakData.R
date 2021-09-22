#1) Clear memory
{
  rm(list=ls())
}

#2) Load packages
{
  library(parallel)
  library(tidyverse)
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  setwd("/media/felix/Elements/12_RangeExpansionGenomics")
}

#3) Prepare the gene ontology file
{
  #3.1) Create a list with the gene id's
  GeneFunctions_Complete <- read_delim("2_data/TetThermRefGenome/GeneFunctions_Complete", 
                                       "\t", escape_double = FALSE, trim_ws = TRUE)
  geneIDs <- as.data.frame(unique(c(GeneFunctions_Complete$`2006_Gene_ID`, GeneFunctions_Complete$`2008_Gene_ID`)))
  colnames(geneIDs) <- "gene ID"
  write_tsv(geneIDs, path="2_data/TetThermRefGenome/geneIDs.txt", col_names = F)
  
  #3.2) For this list, a query is performed at GO miner, with following parameters
  {
    #Input:gene ID list
    #Changed list: gene ID list
    #Taxon ID: 312017 (T. thermophila SB210, associated with the reference genome)
    #All other inout: standard input
  }
  
  #3.3) Transform output of query (.gce file) in format usable by GOWINDA
  path.script <- paste(getwd(), "5_Packages/GoMiner2FuncAssociate.py", sep="/")
  gce.file <- paste(getwd(), "2_data/GOWINDA_Data/GO_Query/geneIDs.txt1474171867.dir/geneIDs.txt.dir/geneIDs.txt.change.gce", sep="/")
  output.file <- paste(getwd(), "2_data/GOWINDA_Data/GOList_GOminer.txt", sep="/")
  system(paste("python ", path.script, " --input ", gce.file, " > ", output.file))
}

#3) Define paths to be used for the analysis
{
  #Define Gowinda path
  gowinda.path <- paste(getwd(), "5_Packages/Gowinda-1.12.jar", sep="/")
  
  #Define location of subsetted gwas file
  gwas.file <- paste(getwd(), "2_data/R_DataFiles/PoPoolationData/BonferonniData.gwas", sep="/")
  
  #Define location of gene annotation file (.gtf) containing the entries for the exons
  gtf.path <- paste(getwd(),"2_data/TetThermRefGenome/Simplifiedgtf.gtf", sep="/")
  
  #Define the path to the file with the feneset data from GOminer
  geneset.path <-  paste(getwd(), "2_data/GOWINDA_Data/GOList_GOminer.txt", sep="/")
  
  #Define the path to the unfiltered GWAS data file
  unfiltered <- paste(getwd(), "2_data/PopoolationData/Alldata_try2.cmh_simplified.gwas", sep="/")
}

#4) Perform GOWINDA analysis for bonferroni corrected data
{
  #4.1) Perform basic analysis
  {
    #Define the output file
    output.file <- paste(getwd(), "2_data/GOWINDA_Data/BonferroniData/BasicAnalysis.txt", sep="/")
    
    #Prepare the code to run
    code <- paste("java -Xmx4g -jar ", gowinda.path, " --snp-file ", unfiltered, " --candidate-snp-file ", gwas.file, " --gene-set-file ", geneset.path, " --annotation-file ", gtf.path, " --simulations 100000 --min-significance 1 --gene-definition gene --threads 6 --output-file ", output.file, " --mode gene --min-genes 1 ", sep="")
    
    #Run the code
    system(code)
  }
  
  #4.2) Perform analysis including regulatory regions
  {
    #Define the output file
    output.file <- paste(getwd(), "2_data/GOWINDA_Data/BonferroniData/Analysis_with_Regulatory.txt", sep="/")
    
    #Prepare the code to run
    code <- paste("java -Xmx4g -jar ", gowinda.path, " --snp-file ", unfiltered, " --candidate-snp-file ", gwas.file, " --gene-set-file ", geneset.path, " --annotation-file ", gtf.path, " --simulations 100000 --min-significance 1 --gene-definition updownstream2000 --threads 6 --output-file ", output.file, " --mode gene --min-genes 1 ", sep="")
    
    #Run the code
    system(code)
  }
  
  #4.3) Perform high resolution GO enrichment analysis
  {
    #Define the output file
    output.file <- paste(getwd(), "2_data/GOWINDA_Data/BonferroniData/Analysis_HighRes.txt", sep="/")
    
    #Prepare the code to run
    code <- paste("java -Xmx4g -jar ", gowinda.path, " --snp-file ", unfiltered, " --candidate-snp-file ", gwas.file, " --gene-set-file ", geneset.path, " --annotation-file ", gtf.path, " --simulations 100000 --min-significance 1 --gene-definition updownstream2000 --threads 6 --output-file ", output.file, " --mode snp --min-genes 1 ", sep="")
   
    #Run the code
    system(code)
  }
}

#5) Perform GOWINDA analysis for peak data
{
  gwas.file <- paste(getwd(), "2_data/R_DataFiles/PoPoolationData/PeakData.gwas", sep="/")
  #5.1) Perform basic analysis
  {
    #Define the output file
    output.file <- paste(getwd(), "2_data/GOWINDA_Data/PeakData/BasicAnalysis.txt", sep="/")
    
    #Prepare the code to run
    code <- paste("java -Xmx4g -jar ", gowinda.path, " --snp-file ", unfiltered, " --candidate-snp-file ", gwas.file, " --gene-set-file ", geneset.path, " --annotation-file ", gtf.path, " --simulations 100000 --min-significance 1 --gene-definition gene --threads 6 --output-file ", output.file, " --mode gene --min-genes 1 ", sep="")
    
    #Run the code
    system(code)
  }
  
  #5.2) Perform analysis including regulatory regions
  {
    #Define the output file
    output.file <- paste(getwd(), "2_data/GOWINDA_Data/PeakData/Analysis_with_Regulatory.txt", sep="/")
    
    #Prepare the code to run
    code <- paste("java -Xmx4g -jar ", gowinda.path, " --snp-file ", unfiltered, " --candidate-snp-file ", gwas.file, " --gene-set-file ", geneset.path, " --annotation-file ", gtf.path, " --simulations 100000 --min-significance 1 --gene-definition updownstream2000 --threads 6 --output-file ", output.file, " --mode gene --min-genes 1 ", sep="")
    
    #Run the code
    system(code)
  }
  
  #5.3) Perform high resolution GO enrichment analysis
  {
    #Define the output file
    output.file <- paste(getwd(), "2_data/GOWINDA_Data/PeakData/Analysis_HighRes.txt", sep="/")
    
    #Prepare the code to run
    code <- paste("java -Xmx4g -jar ", gowinda.path, " --snp-file ", unfiltered, " --candidate-snp-file ", gwas.file, " --gene-set-file ", geneset.path, " --annotation-file ", gtf.path, " --simulations 100000 --min-significance 1 --gene-definition updownstream2000 --threads 6 --output-file ", output.file, " --mode snp --min-genes 1 ", sep="")
    
    #Run the code
    system(code)
  }
}