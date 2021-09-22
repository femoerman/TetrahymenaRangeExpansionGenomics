#1) Clear memory
{
  rm(list=ls())
}

#2) Load packages and set working directory
{
  library(parallel)
  library(tidyverse)
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  setwd("/media/felix/Elements/12_RangeExpansionGenomics")
}

#3) Load datasets
{
  peaks.gowinda <-read_delim("2_data/GOWINDA_Data/PeakData//Analysis_HighRes.txt", 
                             "\t", escape_double = FALSE, col_names = FALSE, 
                             trim_ws = TRUE)
  colnames(peaks.gowinda) <- c("GO_term", "average_gene_number", "snp_number", "p_value", "FDR", "unique_genes", "most_genes", "total_genes", "Description_GO_term", "gene_IDs")
  
  
  bonferroni.gowinda <-read_delim("2_data/GOWINDA_Data/BonferroniData/Analysis_HighRes.txt", 
                                  "\t", escape_double = FALSE, col_names = FALSE, 
                                  trim_ws = TRUE)
  colnames(bonferroni.gowinda) <- c("GO_term", "average_gene_number", "snp_number", "p_value", "FDR", "unique_genes", "most_genes", "total_genes", "Description_GO_term", "gene_IDs")
  
    
  general.gowinda <- read_delim("2_data/GOWINDA_Data/AFC_General/Analysis_HighRes.txt", 
                                "\t", escape_double = FALSE, col_names = FALSE, 
                                trim_ws = TRUE)
  colnames(general.gowinda) <- c("GO_term", "average_gene_number", "snp_number", "p_value", "FDR", "unique_genes", "most_genes", "total_genes", "Description_GO_term", "gene_IDs")

}

#4) Filter the datasets, to keep only GO terms that have a false detection rate < 5% (FDR<0.05)
{
  peaks.gowinda <- filter(peaks.gowinda, FDR<0.05)
  write_tsv(peaks.gowinda, path = "4_output/2_CMH_test_GO_Enrichment/peaks_gowinda.txt")
  bonferroni.gowinda <- filter(bonferroni.gowinda, FDR<0.05)
  write_tsv(bonferroni.gowinda, path = "4_output/2_CMH_test_GO_Enrichment/bonferroni_gowinda.txt")
  general.gowinda <- filter(general.gowinda, FDR<0.05)
  write_tsv(general.gowinda, path = "4_output/2_CMH_test_GO_Enrichment/general_gowinda.txt")
}

#5) Print 25 most significant GO terms for each list
{
  bonf <- bonferroni.gowinda %>% arrange(FDR) %>% select(Description_GO_term)
  peak <- peaks.gowinda %>% arrange(FDR) %>% select(Description_GO_term)
  gen <- general.gowinda %>% arrange(FDR) %>% select(Description_GO_term)
}