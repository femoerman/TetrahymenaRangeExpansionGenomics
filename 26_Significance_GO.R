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
 setwd("/media/felix/Data/Felix_Data/Science/PhD/12_RangeExpansionGenomics")
}

#4) Load datasets
{
  load(file = "4_output/2_CMH_test_GO_Enrichment/Summary.RData")
  
  #Include 0 values
  data.sum2 <- data.sum %>% select(-Percentage) %>% spread(key = Dataset, value=Count)
  data.sum2[is.na(data.sum2)] <- 0
  data.sum3 <- data.sum2 %>% gather(key=Dataset, value = Count, -Category)

  
}

#7) try 2
{
  bonferroni_gowinda <- read_delim("4_output/2_CMH_test_GO_Enrichment/bonferroni_gowinda.txt", 
                                   "\t", escape_double = FALSE, trim_ws = TRUE)
  
  peaks_gowinda <- read_delim("4_output/2_CMH_test_GO_Enrichment/peaks_gowinda.txt", 
                              "\t", escape_double = FALSE, trim_ws = TRUE)
  
  general_gowinda <- read_delim("4_output/2_CMH_test_GO_Enrichment/general_gowinda.txt", 
                                "\t", escape_double = FALSE, trim_ws = TRUE)
  
  #Read in data
  annot.list <- read_csv("4_output/2_CMH_test_GO_Enrichment/GO_list_groups.csv")
  
  #Assign a treatment to the three lists of GO terms, and combine
  bonferroni_gowinda$Dataset <- "Bonferroni"
  general_gowinda$Dataset <- "General"
  peaks_gowinda$Dataset <- "Peaks"
  
  data.all <- rbind(bonferroni_gowinda, general_gowinda, peaks_gowinda)
  
  annot.list <- as.data.frame(annot.list)
  rownames(annot.list) <- annot.list$GO.list
  data.all$Category <- ""
  for (i in 1:nrow(data.all)){
    go <- as.character(data.all[i, "Description_GO_term"])
    data.all[i, "Category"] <- annot.list[go, "Category"]
  }
  
  #Relevel Dataset
  data.all$Dataset <- relevel(as.factor(data.all$Dataset), ref = "General")
}



# #Test1: Only compare peak and general
# {
#   #Filter dataset
#   data.all2 <- filter(data.all, Dataset != "Bonferroni")
#   
#   #Relevel Dataset
#   data.all2$Dataset <- relevel(as.factor(as.character(data.all2$Dataset)), ref = "General")
#   
#   #Create table and perform test
#   data.all.tab <-  table(data.all2$Category, data.all2$Dataset)
#   chisq <- chisq.test(data.all.tab)
#   chisq
#   chisq$observed
#   round(chisq$expected,2)
#   round(chisq$residuals, 3)
#   library(corrplot)
#   corrplot(chisq$residuals, is.cor = FALSE)
#   contrib <- 100*chisq$residuals^2/chisq$statistic
#   round(contrib, 3)
#   # Visualize the contribution
#   corrplot(contrib, is.cor = FALSE)
#   
#   #Do pairwise comparison
#   library(chisq.posthoc.test)
#   chisq.posthoc.test(data.all.tab, method="BH")
# }

#Test2: Compare bonferroni and general datasets
{
  #Filter dataset
  data.all3 <- filter(data.all, Dataset != "Peaks")

  #Relevel Dataset
  data.all3$Dataset <- relevel(as.factor(as.character(data.all3$Dataset)), ref = "General")

  #Create table and perform test
  data.all.tab2 <-  table(data.all3$Category, data.all3$Dataset)
  chisq2 <- chisq.test(data.all.tab2)
  chisq2
  chisq2$observed
  round(chisq2$expected,2)
  round(chisq2$residuals, 3)
  library(corrplot)
  corrplot(chisq2$residuals, is.cor = FALSE)
  contrib2 <- 100*chisq2$residuals^2/chisq2$statistic
  round(contrib2, 3)
  # Visualize the contribution
  corrplot(contrib2, is.cor = FALSE)

  #Do pairwise comparison
  library(chisq.posthoc.test)
  chisq.posthoc.test(data.all.tab2, method = "BH")
}

#Visualize correlation for both groups
{
  #Relevel Dataset
  data.all$Dataset <- relevel(as.factor(as.character(data.all$Dataset)), ref = "General")
  data.all.taball <-  table(data.all$Category, data.all$Dataset)
  chisq3 <- chisq.test(data.all.taball)
  chisq3
  chisq3$observed
  round(chisq3$expected,2)
  round(chisq3$residuals, 3)
  library(corrplot)
  corrplot(chisq3$residuals, is.cor = FALSE, tl.col = "black")
  
  order <- c("Ion transport and binding", "Mitochondrial functioning and oxidoreductase-reactions", "Membrane functioning, transport and structure",
             "Microtubules and cytoskeleton", "Gene expression (transcription and translation)", "Mitosis, DNA repair and chromosome division", "Metabolism, activity and transport: Nucleic acids",
             "Metabolism, activity and transport: Proteins and amino acids", "Metabolism, activity and transport: Other carbohydrates",
             "Ribosomal structure, activity and functioning", "Signaling", "Vesicle structures", "Other")
  
  data.all$Category <-  factor(data.all$Category, levels = order)
  data.all$Dataset <-  factor(data.all$Dataset, levels = c("General", "Peaks", "Bonferroni"))
  data.all.taball <-  table(data.all$Category, data.all$Dataset)
  chisq3 <- chisq.test(data.all.taball)
  chisq3$residuals
  corrplot(chisq3$residuals, is.cor = FALSE, tl.col = "black")
}

#Visualize plots using ggplot
{
  #Get the order for the categories
  order <- c("Ion transport and binding", "Mitochondrial functioning and oxidoreductase-reactions", "Membrane functioning, transport and structure",
             "Microtubules and cytoskeleton", "Gene expression (transcription and translation)", "Mitosis, DNA repair and chromosome division", "Metabolism, activity and transport: Nucleic acids",
             "Metabolism, activity and transport: Proteins and amino acids", "Metabolism, activity and transport: Other carbohydrates",
             "Ribosomal structure, activity and functioning", "Signaling", "Vesicle structures", "Other")
  
  #Get the residual data from the first comparison (General versus peaks) and plot
  {
    #Get the order for the categories
    order <- c("Ion transport and binding", "Mitochondrial functioning and oxidoreductase-reactions", "Membrane functioning, transport and structure",
               "Microtubules and cytoskeleton", "Gene expression (transcription and translation)", "Mitosis, DNA repair and chromosome division", "Metabolism, activity and transport: Nucleic acids",
               "Metabolism, activity and transport: Proteins and amino acids", "Metabolism, activity and transport: Other carbohydrates",
               "Ribosomal structure, activity and functioning", "Signaling", "Vesicle structures", "Other")
    
    #Get the residuals
    data1 <- as.data.frame(chisq$residuals)
    colnames(data1) <- c("Category", "Dataset", "Residuals")
    
    #Order variables
    data1$Category <-  factor(data1$Category, levels = rev(order))
    data1$Dataset <-  factor(data1$Dataset, levels = c("General", "Peaks"))
    data1$comparison <- "General vs Peaks"
    
    #Plot the data
    x1 <- ggplot(data1, aes(x=Dataset, y=Category, colour=Residuals, size=Residuals)) + geom_point() +
      scale_colour_gradient(low="#0571b0", high="#b30000")
  }
  
  
  #Get the residual data from the second comparison (General versus bonferroni)
  {
    #Get the residuals
    data2 <- as.data.frame(chisq2$residuals)
    colnames(data2) <- c("Category", "Dataset", "Residuals")
    
    #Order variables
    data2$Category <-  factor(data2$Category, levels = rev(order))
    data2$Dataset <-  factor(data2$Dataset, levels = c("General", "Bonferroni"))
    data2$comparison <- "General vs Bonferroni"
    
    #Plot the data
    x2 <- ggplot(data2, aes(x=Dataset, y=Category, colour=Residuals, size=Residuals)) + geom_point() +
      scale_colour_gradient(low="#0571b0", high="#b30000")
    
    data.all <- rbind(data1, data2)
    data.all$comparison <- factor(data.all$comparison, levels=c("General vs Peaks", "General vs Bonferroni"))
    ggplot(data.all, aes(x=Dataset, y=Category, colour=Residuals)) + geom_point(aes(size=Residuals)) +
      scale_colour_gradient(low="#0571b0", high="#b30000") + facet_grid(~comparison, scales = "free_x") + guides(size=F) +
      theme_light() + theme(axis.text=element_text(size=12), legend.text=element_text(size=12),legend.title=element_text(size=12), strip.text.x=element_text(12),
                            axis.title=element_text(size=15), strip.text = element_text(size=12)) + ylab("Gene ontology (GO) category")
    ggsave(file="4_output/2_CMH_test_GO_Enrichment/PostHocCategories.png", device = "png", dpi = 320, units = "mm", width=1034/4, height = 578/4)
  }
  
}
