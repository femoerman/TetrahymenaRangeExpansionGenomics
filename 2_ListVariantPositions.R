#1) Clear memory
{
  rm(list=ls())
}

#2) load packages
{
  library(tidyverse)
  library(parallel)
}

#3) Load data of ancestors
{
  #B2086.2
  B2086.2 <- read_delim("2_data/VariantCalling_OnlyDeviating/B2086_ANC_P1_merged.bam_ordered_Mdup.vcf", 
                                                     "\t", escape_double = FALSE, trim_ws = TRUE, 
                                                     skip = 1192)
  B2086.2$sample <- "B2086.2"
  B2086.2 <- B2086.2[, -10]
  
  #CU427.4
  CU427.4 <- read_table2("2_data/VariantCalling_OnlyDeviating/CU427_4ANC_P2_merged.bam_ordered_Mdup.vcf", 
                            skip = 1192)
  CU427.4$sample <- "CU427.4"
  CU427.4 <- CU427.4[, -10]
  
  #CU428.2 
  CU428.2 <- read_table2("2_data/VariantCalling_OnlyDeviating/CU428_2_ANC_P2_merged.bam_ordered_Mdup.vcf", 
  skip = 1192)
  CU428.2$sample <- "CU428.2"
  CU428.2 <- CU428.2[, -10]
  
  #SB3539 
  SB3539 <- read_table2("2_data/VariantCalling_OnlyDeviating/SB3539_ANC_P2_merged.bam_ordered_Mdup.vcf", 
                         skip = 1192)
  SB3539$sample <- "SB3539"
  SB3539 <- SB3539[, -10]
  
  #Concatenate all data
  data.all <- rbind(B2086.2, CU427.4, CU428.2, SB3539)
}

#4) Create a list with all positions to include
{
  #Create position factor
  data.all$position <- paste(data.all$`#CHROM`, data.all$POS)
  
  #Create filtered vector of positions
  positions <- data.all[row.names(unique(data.all[,c("position")])),]
  positions <- select(positions, "#CHROM", "POS")
  
  #Write a file with the positions
  write_tsv(positions, path="2_data/TetThermRefGenome/Positions.txt")
}