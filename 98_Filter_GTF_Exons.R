#1) Clear memory
{
  rm(list=ls())
}

#2) Load packages
{
  library(tidyverse)
  library(parallel)
  '%!in%' <- function(x,y)!('%in%'(x,y))
}

#3) Set working directory
{
  setwd("/media/felix/Elements/12_RangeExpansionGenomics")
}

#Load in data
Gene.Ontology <- read_delim("2_data/TetThermRefGenome/uniprot-taxonomy_312017.tab", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)
#Read in
gtf.source <- read_delim("2_data/TetThermRefGenome/GCF_000189635.1_JCVI-TTA1-2.2_genomic.gtf", 
                         "\t", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE, skip = 3)

#Filter for gene
gtf.filter <- filter(gtf.source, X3=="exon")

#Save output
write_tsv(gtf.filter,col_names = F, path="2_data/TetThermRefGenome/GCF_000189635.1_JCVI-TTA1-2.2_genomic_exonsOnly.gtf")
