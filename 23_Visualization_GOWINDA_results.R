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
  setwd("/media/felix/BackupPlus/PhD/12_RangeExpansionGenomics")
}

#4) Load datasets
{
  bonferroni_gowinda <- read_delim("4_output/2_CMH_test_GO_Enrichment/bonferroni_gowinda.txt", 
                                   "\t", escape_double = FALSE, trim_ws = TRUE)
  
  peaks_gowinda <- read_delim("4_output/2_CMH_test_GO_Enrichment/peaks_gowinda.txt", 
                                   "\t", escape_double = FALSE, trim_ws = TRUE)
  
  general_gowinda <- read_delim("4_output/2_CMH_test_GO_Enrichment/general_gowinda.txt", 
                                   "\t", escape_double = FALSE, trim_ws = TRUE)
}

#5) Create a list with the unique gene ontology terms, to assign a broad category
{
  GO.list <- unique(c(bonferroni_gowinda$Description_GO_term, peaks_gowinda$Description_GO_term, general_gowinda$Description_GO_term))
  
  #Save to txt file, in order to add category
  write_tsv(as.data.frame(GO.list), path="4_output/2_CMH_test_GO_Enrichment/GO_list.txt")
}

#6) load the list with assigned groups, and use it to annotate the three existing files
{
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
  
  #Create annot.list without bonferroni dataset
  annot.list2 <- filter(annot.list, GO.list %in% c(general_gowinda$Description_GO_term, peaks_gowinda$Description_GO_term))
  write_csv(annot.list2, "4_output/2_CMH_test_GO_Enrichment/GO_list_groups2.csv")
}

#7) Summarize data and plot output
{
  #Perform analysis without bonferroni dataset
  data.all <- filter(data.all, Dataset != "Peaks")
  
  #Summarize data
  data.sum <- data.all %>% group_by(Dataset, Category) %>% summarize(Count=n())
  
  #Get a relative proportion
  data.sum$Percentage <- 0
  for (i in 1:nrow(data.sum)){
    if (data.sum[i, "Dataset"] == "Bonferroni"){
      data.sum[i, "Percentage"] <- data.sum[i, "Count"]/693 * 100
    } else if (data.sum[i, "Dataset"] == "Peaks"){
      data.sum[i, "Percentage"] <- data.sum[i, "Count"]/74 * 100
    } else if (data.sum[i, "Dataset"] == "General"){
      data.sum[i, "Percentage"] <- data.sum[i, "Count"]/78 * 100
    }
  }
  
  #Relabel categories so it follows a more logical order
  data.sum$Category <- factor(data.sum$Category, levels = c( "Mitosis, DNA repair and chromosome division","Gene expression (transcription and translation)",
                                                             "Ion transport and binding", "Mitochondrial functioning and oxidoreductase-reactions", 
                                                             "Metabolism, activity and transport: Nucleic acids", "Membrane functioning, transport and structure",
                                                             "Microtubules and cytoskeleton","Ribosomal structure, activity and functioning", 
                                                            "Metabolism, activity and transport: Proteins and amino acids", "Metabolism, activity and transport: Other carbohydrates",
                                                            "Vesicle structures", "Signaling", "Other"))
  data.sum$Dataset <- factor(data.sum$Dataset, levels=c("General", "Bonferroni"))
  data.sum$Dataset <- ifelse(data.sum$Dataset=="General", "General adaptation", "Gradient-specific adaptation")
  data.sum$sign <- ifelse(data.sum$Category %in% c( "Mitosis, DNA repair and chromosome division","Gene expression (transcription and translation)",
                                                    "Ion transport and binding", "Mitochondrial functioning and oxidoreductase-reactions", "Metabolism, activity and transport: Nucleic acids"), "1", "0.5")
  #Visualize data
  {
    ggplot(data.sum, aes(x = Category, y = Percentage, colour=Category, fill=Category)) + geom_bar(stat = "identity")+ facet_grid(Dataset~.) + 
      theme_light() + theme(axis.text=element_text(size=12), legend.text=element_text(size=15),legend.title=element_text(size=16), strip.text.x=element_text(24),
                            axis.title=element_text(size=15), strip.text = element_text(size=15), legend.position = "bottom",
                            axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
    
    
    #Create facets on panels
    library(egg)
    tag_facet2 <- function(p, open = "(", close = ")", tag_pool = letters, x = -Inf, y = Inf, 
                           hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
      
      gb <- ggplot_build(p)
      lay <- gb$layout$layout
      tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
      p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                    vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE)
    }
    
    library(RColorBrewer)
    colourCount = length(unique(data.sum$Category))
    getPalette = colorRampPalette(brewer.pal(9, "Set2"))
    colours <- c("#E41A1C", "#1ef8f1", "#3E8E93", "#4DAF4A", "#7E6E85", "#BA5E6C", "#FF7F00", "#FFD421", "#E1C62F", "#A65628", "#DB728C", "#D789B2", "#999999")
    
    x <- ggplot(data.sum, aes(x = "", y = Percentage, colour=Category, fill=Category, alpha=sign)) + geom_bar(stat = "identity", size=1.5) + facet_wrap(~Dataset) + 
      ylab("Percentage of total enriched GO terms") + 
      scale_alpha_discrete(range = c(0, 1), guide = 'none') +
      geom_text(size = 4.5, position = position_stack(vjust = 0.5), mapping = aes(label=paste(Count, " (", round(Percentage, digits = 2), "%)", sep="")), alpha=1, color="black") +
      theme_light() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                            panel.background = element_blank(), axis.text=element_text(size=12), legend.text=element_text(size=12),legend.title=element_text(size=15), strip.text.x=element_text(24),
                            axis.title=element_text(size=15), strip.text = element_text(size=15), axis.title.x = element_blank(),
                            axis.text.x=element_text(angle=90, hjust=1)) + scale_fill_manual(values=colours) +
      scale_colour_manual(values=colours)
  
    
    x <- tag_facet2(x, 
                     x = -Inf, y = 100, 
                     vjust = 0, hjust = -0.3,
                     open = "", close = "",
                     fontface = 4,
                     size = 6,
                     tag_pool = c("A)", "B)"))  
    x
    
    
    ggsave(x, path = "4_output/2_CMH_test_GO_Enrichment/", filename="Percentages.png", device="png", height=180, width=310, dpi=320, units = "mm")
    
    y <- ggplot(data.sum, aes(x = "", y = Count, fill=Category, colour=Category, alpha=sign)) + geom_bar(stat = "identity", size=1.5) + facet_wrap(~Dataset) + 
      ylab("Cpunt of total enriched GO terms") + 
      scale_alpha_discrete(range = c(0, 1), guide = 'none') +
      geom_text(size = 4.5, position = position_stack(vjust = 0.5), mapping = aes(label=paste(Count, " (", round(Count, digits = 2), "%)", sep="")), alpha=1, color="black") +
      theme_light() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                            panel.background = element_blank(), axis.text=element_text(size=12), legend.text=element_text(size=12),legend.title=element_text(size=15), strip.text.x=element_text(24),
                            axis.title=element_text(size=15), strip.text = element_text(size=15), axis.title.x = element_blank(),
                            axis.text.x=element_text(angle=90, hjust=1)) + scale_fill_manual(values=colours) +
      scale_colour_manual(values=colours)
    
    y <- tag_facet2(y, 
                    x = -Inf, y = 700, 
                    vjust = 0, hjust =0.3,
                    open = "", close = "",
                    fontface = 4,
                    size = 5,
                    tag_pool = c("A)", "B)"))  
    y
    ggsave(y, path = "4_output/2_CMH_test_GO_Enrichment/", filename="Counts.png", device="png", height=300, width=310, dpi=320, units = "mm")
    
  }
  
}

#Save output data
save(data.sum, file = "4_output/2_CMH_test_GO_Enrichment/Summary.RData")


#Create output format to make Latex table
{
  data.sum.roud <- data.sum

  data.sum.2 <- select(data.sum.roud, -Count) %>% spread(key=Dataset, value = Percentage)
  data.sum.2$General <- round(replace_na(data.sum.2$General, 0), digits=2)
  data.sum.2$Peaks <- round(replace_na(data.sum.2$Peaks, 0), digits=2)
  data.sum2 <- as.data.frame(data.sum.2)
  write_tsv(data.sum.2, path = "4_output/2_CMH_test_GO_Enrichment/PercentageTable.txt")
  
  data.sum.3 <- select(data.sum.roud, -Percentage) %>% spread(key=Dataset, value = Count)
  data.sum.3$General <- round(replace_na(data.sum.3$General, 0), digits=2)
  data.sum.3$Peaks <- round(replace_na(data.sum.3$Peaks, 0), digits=2)
  data.sum.3 <- as.data.frame(data.sum.3)
  write_tsv(data.sum.3, path = "4_output/2_CMH_test_GO_Enrichment/CountTable.txt")
}
