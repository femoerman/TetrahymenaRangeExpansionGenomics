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
    query.files <- paste(getwd(), list.files(path = "2_data/QueriedVCF/AllPositions", 
                                             recursive = T, pattern = ".query$", full.names = T), sep="/")
    
    #Annotation files
    annot.files  <- paste(getwd(), list.files(path = "2_data/AnnotationsVCF/AllPositions/", 
                                              recursive = T, pattern = ".annot$", full.names = T), sep="/")
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
        query <- filter(query, position %in% Positions$positions)
        annot <- filter(annot, position %in% Positions$positions)
        
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
        
        #Get the geneID
        fun1 <- function(lst, n){
          sapply(lst, `[`, n)
        }
        annot$geneID <- fun1(geneIDS, 4)
        annot$geneID <- replace_na(annot$geneID, ".")
        
        
        #Get a numeric geneID
        annot <- mutate(annot, geneID_num=substr(geneID, start=8, stop=nchar(geneID)))
        annot$geneID_num <- as.numeric(annot$geneID_num)
        
        #Remove X19 from annot
        annot <- select(annot, -X19)
        
        #Add gene ID to query 
        query$geneID <- annot$geneID
        query$geneID_num <- annot$geneID_num
        
        #Add sample name to the query file
        sample.names <- data.frame(queryfiles <- query.files)
        sample.names$sample <- c("B2086.2", "CU427.4", "CU428.2", "EVO11", "EVO14", "EVO18", "EVO20", "EVO23", "EVO24", "EVO29", "EVO30", "EVO33", "EVO35", "EVO36", "EVO39",
                                 "EVO4", "EVO5", "EVO8", "EVO9", "SB3539")
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
          query$alt.number <- ifelse(query$ALT==".", 0, str_count(query$ALT, ",")+1)
          query$alternatives <-str_split(query$ALT, pattern=",")
          query$alternativecounts <- str_split(query$AD, pattern=",")
          query$total.count <- lapply(query$alternativecounts, as.numeric)
          query$total.count <- unlist(lapply(query$total.count, sum))
          query$ref.count<- as.numeric(unlist(lapply(query$alternativecounts, `[[`, 1)))
          query$total.alt.count <- query$total.count - query$ref.count
          query$alternativecounts <- lapply(query$alternativecounts, c, 0)
          query$first.alt.count <- as.numeric(unlist(lapply(query$alternativecounts, `[[`, 2)))
          query$first.alt.call <- unlist(lapply(query$alternatives, `[[`, 1))
          
          query <- select(query, -alternativecounts, -alternatives)
          
          #Add position (coordinate)
          scaffolds <- as.data.frame(scaffolds)
          rownames(scaffolds) <- scaffolds$`RefSeq-Accn`
          query$coordinate <- scaffolds[query$CHROM, "start"] + query$POS
        }
      }
      
      #Return dataframe
      return(query)
    }
  }
  
  #4.4) Perform parallellized preparation of files, and bind together
  {
    #Parallellized processing
    data.all <- mclapply(1:20, processFiles, query.files, annot.files, scaffolds, mc.cores=6)
    
    #Bind together
    data.all.df <- rbind(data.all[[1]], data.all[[2]], data.all[[3]], data.all[[4]], 
                         data.all[[5]], data.all[[6]], data.all[[7]], data.all[[8]],
                         data.all[[9]], data.all[[10]], data.all[[11]], data.all[[12]],
                         data.all[[13]], data.all[[14]], data.all[[15]], data.all[[16]],
                         data.all[[17]], data.all[[18]], data.all[[19]], data.all[[20]])
    
    #Calculate Allele frequency
    data.all.df$ref.freq <- data.all.df$ref.count/data.all.df$total.count
    data.all.df$alt.freq <- data.all.df$total.alt.count/data.all.df$total.count
    data.all.df$first.alt.freq <- data.all.df$first.alt.count/data.all.df$total.count
    
    #Make variable for coding/non-coding
    data.all.df$expressed <- ifelse(data.all.df$geneID == ".", "non-expressed", "expressed")
    
    #Save temporary file for filtering
    save(data.all.df, file = "2_data/R_DataFiles/AFC_Analysis/DataUnfilteredAFC.RData")
  }
  
  #4.5) Filter variants based on following criteria
  {
    #Reload data
    load("2_data/R_DataFiles/AFC_Analysis/DataUnfilteredAFC.RData")
    
    Positions <- read_delim("2_data/TetThermRefGenome/Positions.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)
    
    #Create a variable entry for exact position
    Positions$positions <- paste(Positions$`#CHROM`, Positions$POS)
    Positions <- filter(Positions, positions %in% data.all.df$position)
    Positions <- as.data.frame(Positions)
    Positions2 <- Positions[which(!duplicated(Positions$positions)), ]
    rownames(Positions2) <- Positions2$positions
    
    #Summarize data frame by positions
    data.all.df.sum <- arrange(data.all.df, position) %>% group_by(position) %>% summarize(minqual = min(QUAL), mindepth=min(total.count), maxdepth=max(total.count),
                                                                                           samples=length(unique(sample)), maxalt=max(first.alt.count),
                                                                                           maxref=max(ref.count), occurrences = n())
    
    #Filter the summarized df based on the quality control criteria
    {
      #1) Minimum quality > 30
      #4) Data available for all 16 samples
      #5) Minimum depth (high quality reads) of 30
      #6) Manixum depth (high quality reads) of 300
      data.all.df.sum.filt <- filter(data.all.df.sum, minqual>=30, samples==20, mindepth>29, maxdepth<301)
    }
    
    #4.5.1) Filter out the positions in the data frame, not included in the filtered summary
    data.all.df.filt <- filter(data.all.df, position %in% data.all.df.sum.filt$position)
    
    #Save intermediate output
    save(data.all.df.filt, file="2_data/R_DataFiles/AFC_Analysis/Filtered_AFC_Data.RData")
  }
  
  #4.6) Save output
  # save(Variants.all, file = "2_data/R_DataFiles/AFC_Analysis/AllMutations.RData")
  # save(Variants.biallellic, file = "2_data/R_DataFiles/AFC_Analysis/BiallellicMutations.RData")
}

#5) Calculate ancestral allele frequency, and than calculate allele frequency change (raw, absolute) and significance for the evolved populations
{
  #Separate data in ancestral and evolved
  data.anc <- filter(data.all.df.filt, gradient=="Ancestor")
  data.evolved <- filter(data.all.df.filt, gradient!="Ancestor")
  
  #Get a list with the ancestral type
  type <- filter(data.anc, TYPE != "REF") %>% filter(!duplicated(position))
  
  #Simplify ancestor dataset, and calculate AFC(Change in ancestral )
  data.anc <- select(data.anc, position, sample, ref.freq)
  data.anc.spread <- spread(data.anc, key=sample, value=ref.freq)
  
  #Load data with ancestor frequencies
  load("2_data/R_DataFiles/ProportionsAncestors_StartExperiment.RData")
  strainprops <- pop_output %>% group_by(strain) %>% summarize(meandens=mean(indiv_per_volume))
  strainprops$prop <- strainprops$meandens/sum(strainprops$meandens)
  strainprops <- as.data.frame(strainprops)
  rownames(strainprops) <- strainprops$strain
  
  #Weigh and sum allele frequencies of ancestors based on the ancestor proportions
  data.anc.spread$AF_anc <- data.anc.spread$B2086.2*strainprops["B2086.2", "prop"] +
    data.anc.spread$CU427.4*strainprops["CU427.4", "prop"] +
    data.anc.spread$CU428.2*strainprops["CU428.2", "prop"] +
    data.anc.spread$SB3539*strainprops["SB3539", "prop"]
  data.anc.spread <- as.data.frame(data.anc.spread)
  rownames(data.anc.spread) <- data.anc.spread$position
  data.anc.spread$type <- type$TYPE
  
  #Calculate allele frequency change for evolved populations
  data.evolved$anc.ref.freq <- data.anc.spread[data.evolved$position, "AF_anc"]
  data.evolved$anc_type <- data.anc.spread[data.evolved$position, "type"]
  data.evolved$AFC_raw <- data.evolved$anc.ref.freq - data.evolved$ref.freq
  data.evolved$AFC_abs <- abs(data.evolved$AFC_raw)
  data.evolved$TYPE <- ifelse(data.evolved$TYPE=="REF", data.evolved$anc_type, data.evolved$TYPE)
  
  #Calculate p-value for binomial test based on the reads and known ancestor frequency
  doBinom <- function(i){
    binom.test(data.evolved[i, "ref.count"], data.evolved[i, "total.count"], data.evolved[i, "anc.ref.freq"])$p.value
  }
  data.evolved <- as.data.frame(data.evolved)
  data.evolved$pval_binom <- unlist(lapply(1:nrow(data.evolved), doBinom))
  
  #Calculate -log10 of the p-value
  data.evolved$p_BinomProb <- -log10(data.evolved$pval_binom)
  
  #Save the output
  save(data.evolved, file="2_data/R_DataFiles/AFC_Analysis/EvolvedData_AFC.RData")
}
