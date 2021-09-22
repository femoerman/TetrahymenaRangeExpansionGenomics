#1) Clear memory and set number of cores
{
  rm(list = ls())
  library(parallel)
  cores <- detectCores() - 1
  library(tidyverse)
  library(scales)
  library(purrr)
  library(ggfortify)
  library(MASS)
  library(janitor)
}

#2) Set working directory
{
  setwd("/media/felix/Seagate/Documenten/PhD/12_RangeExpansionGenomics")
}

#3) Load the data
{
  load(file="2_data/7_Allel_frequencies_Ancestors/Regular_VariantCalling/AFC_Data_Final_Calculation.RData")
}

#4) Filter data
{
  Evo.data <- Evo.data %>% filter(AFC_abs>=0.3)
}

#5) Assign function to genes
{
  source("3_analysis/Functions/01_AddGeneFunction.R")
  Evo.data.func <- assignFunction(Evo.data)
}

#6) do PCA for gene function
{
  #6.1) prepare data
  {
    #Split data by sample
    datasplit <- split(Evo.data.func , f = Evo.data.func$sample)
    datasplit <- lapply(datasplit, mutate, geneFunction=tolower(geneFunction))
    datasplit <- lapply(datasplit, filter, !is.na(geneFunction))
    datasplit <- lapply(datasplit, filter, geneFunction !="hypothetical protein")
    datasplit <- lapply(datasplit, filter, geneFunction !="conserved hypothetical protein")
    
    #Summarize for each population the number variants for each gene function
    summarizeVariants <- function(df.filtered){
      df.filtered %>% group_by(geneFunction) %>% summarize(number=n())
    }
    summarizedFunctions <- lapply(datasplit, summarizeVariants)
    
    #Add variable containing sample name to all data frames
    output.df <- data.frame(geneFunction=NA, number=NA, sample=NA)
    for (i in names(summarizedFunctions)){
      temp <- as.data.frame(summarizedFunctions[i])
      colnames(temp) <- c("geneFunction", "number")
      temp$sample <- i
      output.df <- rbind(output.df, temp)
    }
    
    output.df <- output.df[2:nrow(output.df), ]
    
    #Spread data frame so that columns correspond to the gene functions
    output.df <- spread(output.df, key=geneFunction, value=number)
    
    #Replace NA's with 0
    output.df <- output.df %>% replace(is.na(.), 0)
  }
  
  
  #6.2) Add treatment info
  {
    #Add treatment info to the data (for gene function)
    output.df$gradient <- ifelse(output.df$sample %in% c("MIX_EVO_11", "MIX_EVO_14", "MIX_EVO_18", "MIX_EVO_20",
                                                                 "MIX_EVO_4", "MIX_EVO_5", "MIX_EVO_8", "MIX_EVO_9"), "Gradient", "Uniform")
    
    output.df$reproduction <-  ifelse(output.df$sample %in% c("MIX_EVO_11", "MIX_EVO_14", "MIX_EVO_23", "MIX_EVO_24",
                                                                      "MIX_EVO_4", "MIX_EVO_5", "MIX_EVO_33", "MIX_EVO_35_P3_AACGAGGCCG-CACGGAA"), "Sexual", "Asexual")
    
    output.df$geneFlow <-  ifelse(output.df$sample %in% c("MIX_EVO_23", "MIX_EVO_24", "MIX_EVO_29", "MIX_EVO_30",
                                                                  "MIX_EVO_4", "MIX_EVO_5", "MIX_EVO_8", "MIX_EVO_9"), "Present", "Absent")
  }
  
  #6.3) Plot a PCA
  {
    
    output.df$treatment <- paste(output.df$gradient, output.df$reproduction, output.df$geneFlow)
    df.function <- output.df[2:106]
    function.pca <- prcomp(df.function, scale=F)
    print(function.pca)
    plot(function.pca)
    autoplot(function.pca, data = output.df, colour = 'gradient', loadings=F, loadings.label=F)
    autoplot(function.pca, data = output.df, colour = 'reproduction', loadings=F, loadings.label=F)
    autoplot(function.pca, data = output.df, colour = 'geneFlow', loadings=F, loadings.label=F)
    autoplot(function.pca, data = output.df, colour = 'treatment', loadings=F, loadings.label=F)
  }
}

#7) do PCA for GO (biological function)
{
  #7.1) prepare data
  {
    #Filter out entries where gene ontology is missing
    data.GO.biol <-lapply(datasplit, filter, !is.na(GO_biol))
    
    #Create an output data frame
    output.GO.biol <- data.frame(GO_biol=NA, number=NA, sample=NA)
    
    #Create a loop to make a list of GO_ID's per population
    for (i in names(data.GO.biol)){
      temp <- as.data.frame(data.GO.biol[i])
      if(i=="MIX_EVO_35_P3_AACGAGGCCG-CACGGAA"){
        var <-"MIX_EVO_35_P3_AACGAGGCCG.CACGGAA.GO_biol"
      } else {var <- paste(i, "GO_biol", sep=".")}
      
      GO_list <- as.data.frame(unlist(str_split((temp[, var]), pattern="; ")))
      colnames(GO_list) <- c("GO_biol")
      GO_summ <- data.frame(GO_list) %>% group_by(GO_biol) %>% summarize(number=n())
      GO_summ$sample <- i
      output.GO.biol <- rbind(output.GO.biol, GO_summ)
    }
    
    #Remove the first entry
    output.GO.biol <- output.GO.biol[2:nrow(output.GO.biol), ]
    
    #Spread the data and replace NA with 0
    output.GO.biol <- spread(output.GO.biol, key=GO_biol, value=number)
    output.GO.biol <- output.GO.biol %>% replace(is.na(.), 0)
  }
  
  #7.2) Add treatment info
  {
    #Add treatment info to the data (for gene function)
    output.GO.biol$gradient <- ifelse(output.GO.biol$sample %in% c("MIX_EVO_11", "MIX_EVO_14", "MIX_EVO_18", "MIX_EVO_20",
                                                         "MIX_EVO_4", "MIX_EVO_5", "MIX_EVO_8", "MIX_EVO_9"), "Gradient", "Uniform")
    
    output.GO.biol$reproduction <-  ifelse(output.GO.biol$sample %in% c("MIX_EVO_11", "MIX_EVO_14", "MIX_EVO_23", "MIX_EVO_24",
                                                              "MIX_EVO_4", "MIX_EVO_5", "MIX_EVO_33", "MIX_EVO_35_P3_AACGAGGCCG-CACGGAA"), "Sexual", "Asexual")
    
    output.GO.biol$geneFlow <-  ifelse(output.GO.biol$sample %in% c("MIX_EVO_23", "MIX_EVO_24", "MIX_EVO_29", "MIX_EVO_30",
                                                          "MIX_EVO_4", "MIX_EVO_5", "MIX_EVO_8", "MIX_EVO_9"), "Present", "Absent")
  }
  
  #6.3) Plot a PCA
  {
    output.GO.biol$treatment <- paste(output.GO.biol$gradient, output.GO.biol$reproduction, output.GO.biol$geneFlow)
    df.function <- output.df[2:43]
    function.pca <- prcomp(df.function, scale=T)
    print(function.pca)
    plot(function.pca)
    autoplot(function.pca, data = output.df, colour = 'gradient', loadings=F, loadings.label=F)
    autoplot(function.pca, data = output.df, colour = 'reproduction', loadings=F, loadings.label=F)
    autoplot(function.pca, data = output.df, colour = 'geneFlow', loadings=F, loadings.label=F)
    autoplot(function.pca, data = output.df, colour = 'treatment', loadings=F, loadings.label=F)
  }
}

#Summarize data by counting the total occurrences per treatment (gradient/uniform) and total occurrences
{
  #Start with the whole dataset, filter out hypothetical proteins and unknown functions
  data.filt <- Evo.data.func %>% filter(geneFunction != "conserved hypothetical protein", geneFunction != "hypothetical protein", !is.na(geneFunction))
  
  #Summarize per treatment group, and count number of samples and occurences
  data.summ <- data.filt %>% group_by(geneFunction, gradient) %>% summarize(occurrences =n(), samples=length(unique(sample)), meanFreq=mean(first.alt.perc))
  
  #Spread out data to compare numbers better
  
  myspread <- function(df, key, value) {
    # quote key
    keyq <- rlang::enquo(key)
    # break value vector into quotes
    valueq <- rlang::enquo(value)
    s <- rlang::quos(!!valueq)
    df %>% gather(variable, value, !!!s) %>%
      unite(temp, !!keyq, variable) %>%
      spread(temp, value)
  }
  
  data.spread <- myspread(data.summ, gradient, c(occurrences, samples, meanFreq))
  data.spread <- data.spread %>% replace(is.na(.), 0)
  data.spread <- mutate(data.spread, samplediff=Gradient_samples-Uniform_samples, countdiff=Gradient_occurrences-Uniform_occurrences, freqdiff=Gradient_meanFreq-Uniform_meanFreq)
}
