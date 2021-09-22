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

#4) Load datasets
{
  UnifAll <- read_delim("2_data/GOWINDA_Data/Mutations/UnifAllFullAnalysis.txt", 
                                    "\t", escape_double = FALSE, col_names = FALSE, 
                                    trim_ws = TRUE) %>% filter(X5 <= 0.05)
  Unif10 <- read_delim("2_data/GOWINDA_Data/Mutations/Unif10FullAnalysis.txt", 
                        "\t", escape_double = FALSE, col_names = FALSE, 
                        trim_ws = TRUE) %>% filter(X5 <= 0.05)
  GradAll <- read_delim("2_data/GOWINDA_Data/Mutations/GradAllFullAnalysis.txt", 
                        "\t", escape_double = FALSE, col_names = FALSE, 
                        trim_ws = TRUE) %>% filter(X5 <= 0.05)
  Grad10 <- read_delim("2_data/GOWINDA_Data/Mutations/Grad10FullAnalysis.txt", 
                        "\t", escape_double = FALSE, col_names = FALSE, 
                        trim_ws = TRUE) %>% filter(X5 <= 0.05)
  All10 <- read_delim("2_data/GOWINDA_Data/Mutations/All10FullAnalysis.txt", 
                        "\t", escape_double = FALSE, col_names = FALSE, 
                        trim_ws = TRUE) %>% filter(X5 <= 0.05)
}