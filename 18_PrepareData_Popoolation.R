#This analysis prepares data for use with popoolation2

#1) Clear memory and load packages and set woking directory
{
  rm(list=ls())
  cores<-detectCores()-2
  
  library(parallel)
  library(tidyverse)
  
  #setwd("/media/felix/BackupPlus/MoermanF_Data/Science/PhD/12_RangeExpansionGenomics")
}

#2) Load data
{
  
  #List all the merged bam files from evolved populations
  All.files <- paste(getwd(), list.files(path = "2_data/MergedBam/", recursive = T, pattern = "Mdup.bam$", full.names = T), sep="/")[4:19]
  
  #Order them by sample name
  All.files <- c(All.files[13:16], All.files[1:12])
  
  All.files <-paste(unlist(All.files), sep=" ", collapse = " ")
  
  #Do the mpileup of all populations together, as needed for the PoPoolation analysis
  #Name file for outputting the popoolation data (combined vcf file)
  output1 <- paste(getwd(), "2_data/PopoolationData/Alldata_try2.mpileup", sep="/")
  system(paste("samtools mpileup -B ", All.files, " > ", output1))
  
  output2 <- paste(getwd(), "2_data/PopoolationData/Alldata_try2.sync", sep="/")
  system(paste("java -ea -Xmx7g -jar /media/felix/Elements/12_RangeExpansionGenomics/5_Packages/popoolation2_1201/mpileup2sync.jar --input ", output1, " --output ", output2," --fastq-type sanger --min-qual 30 --threads 6", sep=""))
}

#3) Calculate allele frequency differences between samples
{
  snp.freq <- paste(getwd(), "5_Packages/popoolation2_1201/snp-frequency-diff.pl", sep="/")
  file.in <- paste(getwd(), "2_data/PopoolationData/Alldata_try2.sync", sep="/")
  file.out <- paste(getwd(), "2_data/PopoolationData/Alldata_try2", sep="/")
  system(paste("perl ", snp.freq , " --input ",file.in, " --output-prefix ", file.out, " --min-count 2 --min-coverage 30 --max-coverage 300", sep=""))
}

#4) Calculate Fst per site over genome
{
  Fst.calc <- paste(getwd(), "5_Packages/popoolation2_1201/fst-sliding.pl", sep="/")
  file.in <- paste(getwd(), "2_data/PopoolationData/Alldata_try2.sync", sep="/")
  file.out <- paste(getwd(), "2_data/PopoolationData/Alldata_try2.fst", sep="/")
  #system(paste("perl ", Fst.calc, " --input ", file.in, " --output ", file.out, " --suppress-noninformative --min-count 2 --min-coverage 30 --max-coverage 300 --min-covered-fraction 1 --window-size 1 --step-size 1 --pool-size 500", sep=""))
  
  #Create igv file
  output2 <- paste(getwd(), "2_data/PopoolationData/Alldata_try2.fst.igv", sep="/")
  convert <- paste(getwd(), "5_Packages/popoolation2_1201/export/pwc2igv.pl", sep="/")
  system(paste("perl ", convert, " --input ", file.out, " --output ", output2))
  
}

#5) Calculate sliding window Fst
{
  Fst.calc <- paste(getwd(), "5_Packages/popoolation2_1201/fst-sliding.pl", sep="/")
  file.in <- paste(getwd(), "2_data/PopoolationData/Alldata_try2.sync", sep="/")
  file.out <- paste(getwd(), "2_data/PopoolationData/Alldata_try2.fst_slidingWindow", sep="/")
  #system(paste("perl ", Fst.calc, " --input ", file.in, " --output ", file.out, " --suppress-noninformative --min-count 2 --min-coverage 30 --max-coverage 300 --min-covered-fraction 1 --window-size 1000 --step-size 1000 --pool-size 1000", sep=""))
  
  #Create igv file
  output2 <- paste(getwd(), "2_data/PopoolationData/Alldata_try2.fst_slidingWindow.igv", sep="/")
  convert <- paste(getwd(), "5_Packages/popoolation2_1201/export/pwc2igv.pl", sep="/")
  system(paste("perl ", convert, " --input ", file.out, " --output ", output2))
}

#6) Perform Fisher's exact test to compare significance
{
  Fisher <- paste(getwd(), "5_Packages/popoolation2_1201/fisher-test.pl", sep="/")
  file.in <- paste(getwd(), "2_data/PopoolationData/Alldata_try2.sync", sep="/")
  file.out <- paste(getwd(), "2_data/PopoolationData/Alldata_try2.fet", sep="/")
  #system(paste("perl ", Fisher, " --input ", file.in, " --output ", file.out, " --min-count 2 --min-coverage 30 --max-coverage 300 --suppress-noninformative", sep=""))
  
  #Convert to IGV format
  output2 <- paste(getwd(), "2_data/PopoolationData/Alldata_try2.fet.igv", sep="/")
  convert <- paste(getwd(), "5_Packages/popoolation2_1201/export/pwc2igv.pl", sep="/")
  system(paste("perl ", convert, " --input ", file.out, " --output ", output2))
}

#7) Perform CMH test to compare populations in gradient/homogeneous environment
{
  CMH <- paste(getwd(), "5_Packages/popoolation2_1201/cmh-test.pl", sep="/")
  file.in <- paste(getwd(), "2_data/PopoolationData/Alldata_try2.sync", sep="/")
  file.out <- paste(getwd(), "2_data/PopoolationData/Alldata_try2.cmh", sep="/")
  #system(paste("perl ", CMH, " --input ", file.in, " --output ", file.out, " --min-count 16 --min-coverage 30 --max-coverage 300 --population 1-9,2-10,3-11,4-12,5-13,6-14,7-15,8-16", sep=""))
  
  #Convert to IGV gwas format
  output2 <- paste(getwd(), "2_data/PopoolationData/Alldata_try2.cmh.gwas", sep="/")
  convert <- paste(getwd(), "5_Packages/popoolation2_1201/export/cmh2gwas.pl", sep="/")
  system(paste("perl ", convert, " --input ", file.out, " --output ", output2, " --min-pvalue 1.0e-50"))
}

#8) Calculate Fst for genes
{
  genewise <- paste(getwd(), "5_Packages/popoolation2_1201/create-genewise-sync.pl", sep="/")
  file.in <- paste(getwd(), "2_data/PopoolationData/Alldata_try2.sync", sep="/")
  file.out <- paste(getwd(), "2_data/PopoolationData/Alldata_try2_genes.sync", sep="/")
  gtf.file <- paste(getwd(), "2_data/TetThermRefGenome/GCF_000189635.1_JCVI-TTA1-2.2_genomic_genesOnly.gtf", sep="/")
  system(paste("perl ", genewise, " --input ", file.in, " --gtf ", gtf.file,  " --output ", file.out, sep=""))
}