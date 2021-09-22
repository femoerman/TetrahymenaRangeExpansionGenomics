#1) Clear memory
{
  rm(list=ls())
}

#2) Load packages
{
  library(parallel)
  library(tidyverse)
  '%!in%' <- function(x,y)!('%in%'(x,y))
}

#3) Set working directory
{
  setwd("/media/felix/Elements/12_RangeExpansionGenomics")
}

#4) Load and prepare data
{
  #4.1) Load positions where variants are present in the ancestors
  {
    #Load data
    Positions <- read_delim("2_data/TetThermRefGenome/Positions.txt", 
                                 "\t", escape_double = FALSE, trim_ws = TRUE)
    
    #Create a variable entry for exact position
    Positions$positions <- paste(Positions$`#CHROM`, Positions$POS)
    
    #Read in scaffold data
    scaffolds <- read_delim("2_data/TetThermRefGenome/GCF_000189635.1_JCVI-TTA1-2.2_ScaffoldLength.txt", 
               "\t", escape_double = FALSE, trim_ws = TRUE)
    scaffolds$name <- scaffolds$`Sequence-Name`
    scaffolds <- scaffolds[order(scaffolds$name),] 
    scaffolds$start <- 1
    for (i in 2:nrow(scaffolds)){
      scaffolds[i, "start"] <- scaffolds[i-1, "start"] + scaffolds[i-1, "Sequence-Length"]
    }
    
  }
  
  #4.2) Create a list for the query files and annotation files
  {
    #Query files
    query.files <- paste(getwd(), list.files(path = "2_data/QueriedVCF/VariantsOnly/", 
                                                           recursive = T, pattern = ".query$", full.names = T), sep="/")[4:19]
    
    #Annotation files
    annot.files  <- paste(getwd(), list.files(path = "2_data/AnnotationsVCF/VariantsOnly/", 
                                              recursive = T, pattern = ".annot$", full.names = T), sep="/")[4:19]
  }
  
  #4.3) Create a function to process 1 pair of query and annotation files
  {
    processFiles <- function(i, query.files, annot.files, scaffolds){
      
      #Get the files names for the query/annotation files
      query.file <- query.files[i]
      annot.file <- annot.files[i]
      
      #Read in and prepare the files
      {
        #Query file
        query <- read_delim(query.file, 
                            "\t", escape_double = FALSE, trim_ws = TRUE)
        query$position <- paste(query$`# [1]CHROM`, query$'[2]POS')
        query <- arrange(query, position)
        
        #Annotation file
        annot <- read_delim(annot.file, 
                            "\t", escape_double = FALSE, col_names = FALSE, 
                            trim_ws = TRUE)
        annot$position <- paste(annot$X1, annot$X2)
        annot <- arrange(annot, position)
        
        #Keep only entries at sites where there are variants in the ancestral populations
        query <- filter(query, position %!in% Positions$positions)
        annot <- filter(annot, position %!in% Positions$positions)
        
        #Remove duplicated entries and ensure positions correspond
        #First order by position and type
        query <- arrange(query, position, `[7]TYPE`)
        query <- query[which(!duplicated(query$position)), ]
        annot <- arrange(annot, position)
        annot <- annot[which(!duplicated(annot$position)), ]
        # query <- filter(query, position %!in% duplicatedEntries$position)
        annot <- filter(annot, position %in% query$position)
        query <- filter(query, position %in% annot$position)
        
        #Remove unneccesary columns of the annotation file
        annot <- select(annot, position, X19)
        
        #Extract the gene id from the annotation file
        geneIDS <- str_split(annot$X19, patter='"')
        annot$geneID <- "."
        annot$geneID_num <- -1
        
        for (i in 1:length(geneIDS)){
          temp <- unlist(geneIDS[i])
          if(length(temp)!=1){
            id <- temp[4]
            annot[i, "geneID"] <- id
            annot[i, "geneID_num"] <- as.numeric(substr(id, start=8, stop=nchar(id)))
          }
        }
        
        #Remove X19 from annot
        annot <- select(annot, -X19)
        
        #Add gene ID to query 
        query$geneID <- annot$geneID
        query$geneID_num <- annot$geneID_num
        
        #Add sample name to the query file
        sample.names <- data.frame(queryfiles <- query.files)
        sample.names$sample <- c("EVO11", "EVO14", "EVO18", "EVO20", "EVO23", "EVO24", "EVO29", "EVO30", "EVO33", "EVO35", "EVO36", "EVO39",
                          "EVO4", "EVO5", "EVO8", "EVO9")
        query$sample <- sample.names[which(sample.names$queryfiles....query.files==query.file), "sample"]
        
        #Add treatment info to the query file
        {
          query$gradient <- ifelse(query$sample %in% c("EVO11", "EVO14", "EVO18", "EVO20",
                                                             "EVO4", "EVO5", "EVO8", "EVO9"), "Gradient", 
                                   ifelse(query$sample %in% c("SB3539", "B2086.2", "CU427.4", "CU428.2"), "Ancestor", "Uniform"))
          
          query$reproduction <- ifelse(query$sample %in% c("EVO11", "EVO14", "EVO24", "EVO23",
                                                       "EVO4", "EVO5", "EVO33", "EVO35"), "Sexual", 
                                   ifelse(query$sample %in% c("SB3539", "B2086.2", "CU427.4", "CU428.2"), "Ancestor", "Asexual"))
          
          query$geneFlow <- ifelse(query$sample %in% c("EVO29", "EVO30", "EVO24", "EVO23",
                                                           "EVO4", "EVO5", "EVO8", "EVO9"), "Present", 
                                       ifelse(query$sample %in% c("SB3539", "B2086.2", "CU427.4", "CU428.2"), "Ancestor", "Absent"))
        }
        
        #Rename column names of query data frame
        colnames(query) <- c("CHROM", "POS", "END", "REF", "ALT", "QUAL", "TYPE", "DP", "AD", "position", "geneID", "geneID_num", "sample", "gradient", "reproduction", "geneFlow")
        
        #Create a loop to calculate all other metrics
        {
          #Define new variables
          query$alt.number <- 0
          query$first.alt.count <- 0
          query$total.alt.count <- 0
          query$ref.count <- 0
          query$total.count <- 0
          query.first.alt.call <- ""
          query$coordinate <- 0
          
          #Perform loop calculating all metrics
          for (i in 1:nrow(query)){
            temp <- query[i, ]
            alt.call <- unlist(str_split(temp$ALT, pattern=","))
            counts <- as.numeric(unlist(str_split(temp$AD, pattern=",")))
            query[i, "alt.number"] <- length(alt.call)
            query[i, "first.alt.count"] <- ifelse(length(counts)==1, 0, counts[2])
            query[i, "ref.count"] <- counts[1]
            query[i, "total.count"] <- sum(counts)
            query[i, "total.alt.count"] <- sum(counts) - counts[1]
            query[i, "first.alt.call"] <- ifelse(length(alt.call)>=1, alt.call[1], "")
            query[i, "coordinate"] <- scaffolds[which(scaffolds$`RefSeq-Accn`==temp$CHROM), "start"] + temp$POS
          }
        }
      }
      
      #Return dataframe
      return(query)
    }
  }
  
  #4.4) Perform parallellized preparation of files, and bind together
  {
    #Parallellized processing
    data.all <- mclapply(1:16, processFiles, query.files, annot.files, scaffolds, mc.cores=6)
    
    #Bind together
    data.all.df <- rbind(data.all[[1]], data.all[[2]], data.all[[3]], data.all[[4]], 
                         data.all[[5]], data.all[[6]], data.all[[7]], data.all[[8]],
                         data.all[[9]], data.all[[10]], data.all[[11]], data.all[[12]],
                         data.all[[13]], data.all[[14]], data.all[[15]], data.all[[16]])
    
    #Calculate Allele frequency
    data.all.df$ref.freq <- data.all.df$ref.count/data.all.df$total.count
    data.all.df$alt.freq <- data.all.df$total.alt.count/data.all.df$total.count
    data.all.df$first.alt.freq <- data.all.df$first.alt.count/data.all.df$total.count
    
    #Make variable for coding/non-coding
    data.all.df$expressed <- ifelse(data.all.df$geneID == ".", "non-expressed", "expressed")
  }
  
  #4.5) Filter variants based on following criteria
  {
    #QUAL>30
    data.all.df <- filter(data.all.df, QUAL>=30)
    
    #Alternative reads >= 1
    data.all.df <- filter(data.all.df, first.alt.count >= 1)
    
    #Total reads >= 30
    data.all.df <- filter(data.all.df, total.count >= 30)
    
    #Make separate files for all versus biallelic sites
    Mutations.all <- data.all.df
    Mutations.biallellic <- filter(data.all.df, alt.number < 2)
  }
  
  #4.6) Save output
  save(Mutations.all, file = "2_data/R_DataFiles/MutationAnalysis/AllMutations.RData")
  save(Mutations.biallellic, file = "2_data/R_DataFiles/MutationAnalysis/BiallellicMutations.RData")
}