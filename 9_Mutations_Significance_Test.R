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
  library(car)
  '%!in%' <- function(x,y)!('%in%'(x,y))
}

#3) Set working directory
{
  #setwd("/media/felix/Elements/12_RangeExpansionGenomics")
  setwd("/media/felix/DataDrive2/Documenten/PhD/12_RangeExpansionGenomics")
}

#4) Load data on mutations, and filter to keep significant sites with AFC>0.3 only (exclude clonal selection)
{
  load("2_data/R_DataFiles/MutationAnalysis/AllMutations_Func.RData")
  Mutations.all.func.filt <- filter(Mutations.all.func, alt.freq>0.1)
  
  #Calculate the percentage of snp versus indel
  Mutations.all.func %>% group_by(TYPE) %>% summarize(count=n())
  Mutations.all.func %>% group_by(expressed) %>% summarize(count=n(), perc=n()/27964*100)
  Mutations.all.func %>% group_by(TYPE, expressed) %>% summarize(count=n())
  range(Mutations.all.func$alt.freq)
  nrow(filter(Mutations.all.func.filt, alt.freq>=0.99))
  
  #Do Binomial test on all data
  expected <- c(0.4787*100, (1-0.4787)*100)
  observed <- as.numeric(unlist(Mutations.all.func %>% group_by(expressed) %>% summarize(perc=n()/27964*100) %>% select(perc)))
  chisq.test(expected, observed)
  
  data.sum <-Mutations.all.func %>% group_by(expressed) %>% summarize(count=n())
  x <- xtabs(count ~ expressed, data=data.sum)
  chisq <- chisq.test(x)
  chisq
  #significantly different
  
  
  file_path <- "http://www.sthda.com/sthda/RDoc/data/housetasks.txt"
  housetasks <- read.delim(file_path, row.names = 1)
  
  binom.test(8873, (8873+19091), p=0.4787)
  binom.test(6923, (6923+16835), p=0.4787)
  binom.test(1950, (1950+2256), p=0.4787)
  chisq.test(8873, 19091, p=0.4778)
}

#5) Compare significance between treatments
{
  #5.1) Look at number of variants, depending on different cutoffs for AFC
  {
    #Prepare data
    data.sum.mut <- group_by(Mutations.all.func.filt, gradient, reproduction, geneFlow, sample) %>% summarize(variants=n(), cutoff=0.1, perc.Indel=sum(TYPE == "INDEL")/n())
    cutoffs <- seq(from=0.2, to=0.9, by=0.1)
    for(i in cutoffs){
      temp <- filter(Mutations.all.func.filt, alt.freq>i) %>% group_by(gradient, reproduction, geneFlow, sample) %>% summarize(variants=n(), cutoff=i, perc.Indel=sum(TYPE == "INDEL")/n())
      data.sum.mut <- rbind(data.sum.mut, temp)
    }
    
    #Plot data
    ggplot(data.sum.mut, aes(x=cutoff, y=variants, group=sample, colour=reproduction, shape=geneFlow)) + geom_point() + geom_line() + facet_wrap(~gradient)
    
    #Do a mixed model with the data, separate for gradient/uniform
    
    model.sign <- lme(data=data.sum.mut, variants~reproduction*geneFlow*gradient*cutoff, random= ~1|sample, method="ML")
    dredge(model.sign, rank=BIC)
    best.model <-  lme(data=data.sum.mut, variants~reproduction*cutoff, random= ~1|sample, method="REML")
    summary(best.model)
    Anova(best.model, contrasts=list(topic=contr.sum, sys=contr.sum), type=3)
    
    #Order facets correctly
    data.sum.mut$`Abiotic conditions` <- ifelse(data.sum.mut$gradient=="Gradient", "pH-gradient", "Uniform pH")
    data.sum.mut$`Abiotic conditions` <- factor(data.sum.mut$`Abiotic conditions`, levels=c("Uniform pH", "pH-gradient"))
    
    x <- ggplot(data.sum.mut, aes(x=cutoff, y=(variants), group=sample, colour=reproduction, shape=geneFlow)) + geom_point(size=2) + geom_line() + facet_wrap(~gradient) + 
      ylab("Number of variants") + xlab("Allele frequency cutoff") + 
      scale_colour_manual(breaks = c("Asexual", "Sexual"), values=c("gold3", "blue4")) + 
      theme_light() + theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(24),
                            axis.title=element_text(size=16), strip.text = element_text(size=16)) 
    x
    ggsave(x, path = "4_output/3_significance_Testing/", filename="SignificanceMutations.png", device="png", height=150, width=250, dpi=320, units = "mm")
    
    #Create prediction data
    pred.data <- expand.grid(reproduction=c("Asexual", "Sexual"), cutoff=seq(from=0.1, to=0.9, by=0.1), gradient=c("Gradient", "Uniform"))
    predictions <- predict(best.model, newdata=pred.data, level=0, se.fit=T) 
    pred.data$mean <- predictions$fit
    pred.data$upper <-  predictions$fit + predictions$se.fit*1.96
    pred.data$lower <-  predictions$fit - predictions$se.fit*1.96
    pred.data$Reproduction <- pred.data$reproduction
    data.sum.mut$Reproduction <- data.sum.mut$reproduction
    
    
    pred.data$`Abiotic conditions` <- ifelse(pred.data$gradient=="Gradient", "pH-gradient", "Uniform pH")
    pred.data$`Abiotic conditions` <- factor(pred.data$`Abiotic conditions`, levels=c("Uniform pH", "pH-gradient"))
    
    x <- ggplot(data.sum.mut, aes(x=cutoff, y=(variants), group=sample, colour=Reproduction)) + geom_point(size=3, alpha=0.3) + geom_line(alpha=0.3) + facet_wrap(~`Abiotic conditions`) +
      ylab("Number of variants ") + xlab("Allele frequency of novel variants (cut-off)") +   
      geom_line(inherit.aes = F, data=pred.data, aes(x=cutoff, y=(mean), group=paste(reproduction, gradient), colour=Reproduction), size=1) + 
      geom_ribbon(inherit.aes = F, data=filter(pred.data), aes(x=cutoff, ymax=upper, ymin=lower, group=paste(reproduction, gradient), fill=Reproduction, colour=NA), size=3, alpha=0.3) +
      scale_colour_manual(breaks = c("Asexual", "Sexual"), values=c("gold3", "blue4")) + ylim(-200, 2900) +
      scale_fill_manual(breaks = c("Asexual", "Sexual"), values=c("gold3", "blue4")) +
      theme_light() + theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(24),
                            axis.title=element_text(size=16), strip.text = element_text(size=16), legend.position="bottom") 
    
    x
    tag_facet2 <- function(p, open = "(", close = ")", tag_pool = letters, x = -Inf, y = Inf, 
                           hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
      
      gb <- ggplot_build(p)
      lay <- gb$layout$layout
      tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
      p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                    vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE)
    }
    
    
    x2 <- tag_facet2(x, 
                     x = -Inf, y = 2700, 
                     vjust = 0, hjust = -0,
                     open = "", close = "",
                     fontface = 4,
                     size = 5,
                     tag_pool = c("A)", "B)"))
    x2
    
    ggsave(x2, path = "4_output/3_significance_Testing/", filename="SignificanceMutations_predict.png", device="png", height=150, width=250, dpi=320, units = "mm")
    
    #Plot the proportion of Indels
    x <- ggplot(data.sum.mut, aes(x=cutoff, y=(perc.Indel), group=sample, colour=Reproduction)) + geom_point(size=3) + geom_line() + facet_wrap(~`Abiotic conditions`, labeller = label_both) + 
      ylab("Proportion of indels among novel variants") + xlab("Allele frequency of novel variants (cut-off)") +   
       scale_colour_manual(breaks = c("Asexual", "Sexual"), values=c("gold3", "blue4")) +
      scale_fill_manual(breaks = c("Asexual", "Sexual"), values=c("gold3", "blue4")) +
      theme_light() + theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(24),
                            axis.title=element_text(size=16), strip.text = element_text(size=16), legend.position="bottom") 
    
    x
    ggsave(x, path = "4_output/3_significance_Testing/", filename="Indels.png", device="png", height=150, width=250, dpi=320, units = "mm")
    
    
    
    # #redo for expressed variants only
    # data.sum.AFC.expressed <- filter(Mutations.all.func.filt, expressed=="expressed") %>% group_by(gradient, reproduction, geneFlow, sample) %>% summarize(variants=n(), cutoff=0.3)
    # cutoffs <- seq(from=0.2, to=0.9, by=0.1)
    # for(i in cutoffs){
    #   temp <- filter(Mutations.all.func.filt, alt.freq>i, expressed=="expressed") %>% group_by(gradient, reproduction, geneFlow, sample) %>% summarize(variants=n(), cutoff=i)
    #   data.sum.AFC.expressed <- rbind(data.sum.AFC.expressed, temp)
    # }
    # model.sign.expr <- lme(data=data.sum.AFC.expressed, variants~reproduction*geneFlow*gradient*cutoff, random= ~1|sample, method="ML")
    # dredge(model.sign.expr, rank=BIC)
    # best.model.expr <-  lme(data=data.sum.AFC.expressed, variants~reproduction*cutoff, random= ~1|sample, method="REML")
    # summary(best.model.expr)
    # Anova(best.model.expr, contrasts=list(topic=contr.sum, sys=contr.sum), type=3)
    # 
    # #redo for non-expressed variants only
    # data.sum.AFC.nonexpressed <- filter(Mutations.all.func.filt, expressed!="expressed") %>% group_by(gradient, reproduction, geneFlow, sample) %>% summarize(variants=n(), cutoff=0.3)
    # cutoffs <- seq(from=0.2, to=0.9, by=0.1)
    # for(i in cutoffs){
    #   temp <- filter(Mutations.all.func.filt, alt.freq>i, expressed!="expressed") %>% group_by(gradient, reproduction, geneFlow, sample) %>% summarize(variants=n(), cutoff=i)
    #   data.sum.AFC.nonexpressed <- rbind(data.sum.AFC.nonexpressed, temp)
    # }
    # model.sign.nonexpr <- lme(data=data.sum.AFC.nonexpressed, variants~reproduction*geneFlow*gradient*cutoff, random= ~1|sample, method="ML")
    # dredge(model.sign.nonexpr, rank=BIC)
    # best.model.nonexpr <-  lme(data=data.sum.AFC.nonexpressed, variants ~ reproduction*cutoff, random= ~1|sample, method="REML")
    # summary(best.model.nonexpr)
    # Anova(best.model.nonexpr, contrasts=list(topic=contr.sum, sys=contr.sum), type=3)
  }
}



#7) Plot fitness change to the fitness of the samples
{
  #Load the fitness data
  load("2_data/R_DataFiles/EVO_Fitness.RData")
  
  #Filter to only keep samples included in genomic analysis
  dd.evo <- filter(dd.evo, ID %in% c(18, 20, 8, 9, 11, 14, 4, 5, 36, 39, 29, 30, 33, 35, 23, 24))
  dd.evo$sample <- ifelse(dd.evo$ID==18, "EVO18",
                          ifelse(dd.evo$ID==20, "EVO20",
                                 ifelse(dd.evo$ID==8, "EVO8",
                                        ifelse(dd.evo$ID==9, "EVO9",
                                               ifelse(dd.evo$ID==11, "EVO11",
                                                      ifelse(dd.evo$ID==14, "EVO14",
                                                             ifelse(dd.evo$ID==4, "EVO4",
                                                                    ifelse(dd.evo$ID==5, "EVO5",
                                                                           ifelse(dd.evo$ID==36, "EVO36",
                                                                                  ifelse(dd.evo$ID==39, "EVO39",
                                                                                         ifelse(dd.evo$ID==29, "EVO29",
                                                                                                ifelse(dd.evo$ID==30, "EVO30",
                                                                                                       ifelse(dd.evo$ID==33, "EVO33",
                                                                                                              ifelse(dd.evo$ID==35, "EVO35",
                                                                                                                     ifelse(dd.evo$ID==23, "EVO23",
                                                                                                                            ifelse(dd.evo$ID==24, "EVO24", NA))))))))))))))))
  
  #Add the fitness to the AFC_data
  data.sum.mut$fitness <- 0
  for(i in 1:nrow(data.sum.mut)){
    data.sum.mut[i, "fitness"] <- dd.evo[which(dd.evo$sample==as.character(data.sum.mut[i, "sample"])), "logratio"]
  }

  data.sum.mut$co2 <- paste("Cut-off:", data.sum.mut$cutoff)
  #Plot the fitness values for each of the cutoffs
  y <- ggplot(data.sum.mut, aes(y=fitness, x=variants, group=sample, colour=reproduction, shape=geneFlow)) + geom_point(size=2)  + facet_grid(co2~`Abiotic conditions`) + 
    xlab(expression(paste("Number of ", italic("de novo "), "variants"))) + ylab("Magnitude of fitness change (phenotypic evolution)") + 
    theme_light() + theme(axis.text=element_text(size=10), legend.text=element_text(size=14),legend.title=element_text(size=14), strip.text.x=element_text(24),
                          axis.title=element_text(size=14), strip.text = element_text(size=11), legend.position = "bottom") +
    scale_color_manual(values=c("gold3", "blue4"), breaks=c("Asexual", "Sexual"), name="Reproduction", labels=c("Asexual", "Sexual")) +
    scale_fill_manual(values=c("gold3", "blue4"), breaks=c("Asexual", "Sexual"), name="Reproduction", labels=c("Asexual", "Sexual"))
  y
  
   #Calculate for each of the panels the correlation between change in growth rate and the number of mutations
  {
    temp <- filter(data.sum.mut, cutoff<1) %>% group_by(`Abiotic conditions`, Reproduction, cutoff) %>% summarise(cor= cor.test(fitness, variants)$estimate, pval = cor.test(fitness, variants)$p.value)
    temp$r2 <- temp$cor^2
    
    temp$text1 <- paste("p =", round(temp$pval, digits=4))
    temp$text2 <- paste("R2 =", round(temp$r2, digits=4))
    temp$text3 <- paste("r =", round(temp$cor, digits=4))
  }
  
  
  library(egg)
  tag_facet2 <- function(p, open = "(", close = ")", tag_pool = letters, x = -Inf, y = Inf, 
                         hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
    
    gb <- ggplot_build(p)
    lay <- gb$layout$layout
    tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
    p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                  vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE)
  }
  
  x2 <- tag_facet2(y, 
                   x = -Inf, y = 0.75, 
                   vjust = 0, hjust = -0,
                   open = "", close = "",
                   fontface = 4,
                   size = 3,
                   tag_pool = c("A)", "B)", "C)", "D)", "E)", "F)", "G)", "H)", "I)", "J)", "K)", "L)", "M)", "N)", "O)", "P)", "Q)", "R)"))
  
  #Order the correlations per panel
  t3 <- arrange(temp, Reproduction, `cutoff`) %>% filter(Reproduction=="Sexual")
  t1 <- arrange(temp, Reproduction, `cutoff`) %>% filter(Reproduction=="Asexual")
  
  #Add the correlations to the plot
  x3 <- tag_facet2(x2, 
                   x = c(rep(0, 5), rep(1600, 13)), y = c(rep(0.7, 5), rep(0.15, 13)),
                   vjust = 0, hjust = -0,
                   open = "", close = "",
                   fontface = 4,
                   size = 4,
                   colour="khaki3",
                   tag_pool = c(t1$text1, "", "", "", ""))
  
  x3 <- tag_facet2(x3, 
                   x = c(rep(0, 5), rep(1600, 13)), y = c(rep(0.55, 5), rep(0.00, 13)), 
                   vjust = 0, hjust = -0,
                   open = "", close = "",
                   fontface = 4,
                   size = 4,
                   colour=rep("khaki3",18), 
                   tag_pool = c(t1$text3, "", "", "", ""))
  
  x3 <- tag_facet2(x3, 
                   c(rep(600, 5), rep(2200, 13)), y = c(rep(0.7, 5), rep(0.15, 13)), 
                   vjust = 0, hjust = -0,
                   open = "", close = "",
                   fontface = 4,
                   size = 4,
                   colour= c(rep("cornflowerblue",17), "blue4"),
                   tag_pool = c(t3$text1, "", "", "", ""))
  
  x3 <- tag_facet2(x3, 
                   c(rep(600, 5), rep(2200, 13)), y = c(rep(0.55, 5), rep(0, 13)), 
                   vjust = 0, hjust = -0,
                   open = "", close = "",
                   fontface = 4,
                   size = 4,
                   colour=c(rep("cornflowerblue",17), "blue4"),
                   tag_pool = c(t3$text3, "", "", "", ""))
  
  ggsave(x3, path = "4_output/3_significance_Testing/", filename="Fitness_Mut.png", device = "png",  height=250, width=250, dpi=320, units = "mm")
  
}

# #Repeat with generation scaled data
# {
#   #load the generation data
#   load("2_data/R_DataFiles/Generations.RData")
#   
#   #5.1) Look at number of variants, depending on different cutoffs for AFC
#   {
#     #Prepare data
#     data.sum.mut <- group_by(Mutations.all.func.filt, gradient, reproduction, geneFlow, sample) %>% summarize(variants=n(), cutoff=0.1, perc.Indel=sum(TYPE == "INDEL")/n())
#     cutoffs <- seq(from=0.2, to=0.9, by=0.1)
#     for(i in cutoffs){
#       temp <- filter(Mutations.all.func.filt, alt.freq>i) %>% group_by(gradient, reproduction, geneFlow, sample) %>% summarize(variants=n(), cutoff=i, perc.Indel=sum(TYPE == "INDEL")/n())
#       data.sum.mut <- rbind(data.sum.mut, temp)
#     }
#     
#     #Scale by number of generations
#     sum2$sample <- paste("EVO", sum2$ID, sep="")
#     sum2 <- as.data.frame(sum2)
#     row.names(sum2) <- sum2$sample
#     data.sum.mut$generations <- 0
#     for (i in 1:nrow(data.sum.mut)){
#       ID <- as.character(data.sum.mut[i, "sample"])
#       data.sum.mut[i, "generations"] <- as.numeric(sum2[ID, "totalGenerations"])
#     }
#     data.sum.mut$scaledVariants <- data.sum.mut$variants/data.sum.mut$generations
#     
#     #Plot data
#     ggplot(data.sum.mut, aes(x=cutoff, y=scaledVariants, group=sample, colour=reproduction, shape=geneFlow)) + geom_point() + geom_line() + facet_wrap(~gradient)
#     
#     #Do a mixed model with the data, separate for gradient/uniform
#     
#     model.sign <- lme(data=data.sum.mut, scaledVariants~reproduction*geneFlow*gradient*cutoff, random= ~1|sample, method="ML")
#     dredge(model.sign, rank=BIC)
#     best.model <-  lme(data=data.sum.mut, scaledVariants~cutoff*gradient + reproduction*cutoff, random= ~1|sample, method="REML")
#     summary(best.model)
#     Anova(best.model, contrasts=list(topic=contr.sum, sys=contr.sum), type=3)
#     
#     x <- ggplot(data.sum.mut, aes(x=cutoff, y=(scaledVariants), group=sample, colour=reproduction, shape=geneFlow)) + geom_point(size=2) + geom_line() + facet_wrap(~gradient) + 
#       ylab("Number of variannts") + xlab("Allele frequency cutoff") + 
#       scale_colour_manual(breaks = c("Asexual", "Sexual"), values=c("gold3", "blue4")) + 
#       theme_light() + theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(24),
#                             axis.title=element_text(size=16), strip.text = element_text(size=16)) 
#     x
#     #ggsave(x, path = "4_output/3_significance_Testing/", filename="SignificanceMutations.png", device="png", height=150, width=250, dpi=320, units = "mm")
#     
#     #Create prediction data
#     pred.data <- expand.grid(reproduction=c("Asexual", "Sexual"), cutoff=seq(from=0.1, to=0.9, by=0.1), gradient=c("Gradient", "Uniform"))
#     predictions <- predict(best.model, newdata=pred.data, level=0, se.fit=T) 
#     pred.data$mean <- predictions$fit
#     pred.data$upper <-  predictions$fit + predictions$se.fit*1.96
#     pred.data$lower <-  predictions$fit - predictions$se.fit*1.96
#     pred.data$Reproduction <- pred.data$reproduction
#     data.sum.mut$Reproduction <- data.sum.mut$reproduction
# 
#     
#     x <- ggplot(data.sum.mut, aes(x=cutoff, y=(scaledVariants), group=sample, colour=Reproduction)) + geom_point(size=1.5, alpha=0.3) + geom_line(alpha=0.3) + facet_wrap(~`Abiotic conditions`, labeller = label_both) +
#       ylab("Number of variants ") + xlab("Allele frequency of novel variants (cut-off)") +   
#       geom_line(inherit.aes = F, data=pred.data, aes(x=cutoff, y=(mean), group=paste(reproduction, gradient), colour=Reproduction), size=1) + 
#       geom_ribbon(inherit.aes = F, data=filter(pred.data), aes(x=cutoff, ymax=upper, ymin=lower, group=paste(reproduction, gradient), fill=Reproduction, colour=NA), size=3, alpha=0.3) +
#       scale_colour_manual(breaks = c("Asexual", "Sexual"), values=c("gold3", "blue4")) +
#       scale_fill_manual(breaks = c("Asexual", "Sexual"), values=c("gold3", "blue4")) +
#       theme_light() + theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(24),
#                             axis.title=element_text(size=16), strip.text = element_text(size=16), legend.position="bottom") 
#     
#     x
#     #ggsave(x, path = "4_output/3_significance_Testing/", filename="SignificanceMutations_predict.png", device="png", height=150, width=250, dpi=320, units = "mm")
#     
#     #Plot the proportion of Indels
#     x <- ggplot(data.sum.mut, aes(x=cutoff, y=(perc.Indel), group=sample, colour=Reproduction)) + geom_point(size=1.5) + geom_line() + facet_wrap(~`Abiotic conditions`, labeller = label_both) + 
#       ylab("Proportion of indels among novel variants") + xlab("Allele frequency of novel variants (cut-off)") +   
#       scale_colour_manual(breaks = c("Asexual", "Sexual"), values=c("gold3", "blue4")) +
#       scale_fill_manual(breaks = c("Asexual", "Sexual"), values=c("gold3", "blue4")) +
#       theme_light() + theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(24),
#                             axis.title=element_text(size=16), strip.text = element_text(size=16), legend.position="bottom") 
#     
#     x
#   }
# }
# 
# #8) Plot fitness change to the fitness of the scaledsamples
# {
#   #Load the fitness data
#   load("2_data/R_DataFiles/EVO_Fitness.RData")
#   
#   #Filter to only keep samples included in genomic analysis
#   dd.evo <- filter(dd.evo, ID %in% c(18, 20, 8, 9, 11, 14, 4, 5, 36, 39, 29, 30, 33, 35, 23, 24))
#   dd.evo$sample <- ifelse(dd.evo$ID==18, "EVO18",
#                           ifelse(dd.evo$ID==20, "EVO20",
#                                  ifelse(dd.evo$ID==8, "EVO8",
#                                         ifelse(dd.evo$ID==9, "EVO9",
#                                                ifelse(dd.evo$ID==11, "EVO11",
#                                                       ifelse(dd.evo$ID==14, "EVO14",
#                                                              ifelse(dd.evo$ID==4, "EVO4",
#                                                                     ifelse(dd.evo$ID==5, "EVO5",
#                                                                            ifelse(dd.evo$ID==36, "EVO36",
#                                                                                   ifelse(dd.evo$ID==39, "EVO39",
#                                                                                          ifelse(dd.evo$ID==29, "EVO29",
#                                                                                                 ifelse(dd.evo$ID==30, "EVO30",
#                                                                                                        ifelse(dd.evo$ID==33, "EVO33",
#                                                                                                               ifelse(dd.evo$ID==35, "EVO35",
#                                                                                                                      ifelse(dd.evo$ID==23, "EVO23",
#                                                                                                                             ifelse(dd.evo$ID==24, "EVO24", NA))))))))))))))))
#   
#   #Add the fitness to the AFC_data
#   data.sum.mut$fitness <- 0
#   for(i in 1:nrow(data.sum.mut)){
#     data.sum.mut[i, "fitness"] <- dd.evo[which(dd.evo$sample==as.character(data.sum.mut[i, "sample"])), "logratio"]
#   }
#   
#   #Plot the fitness values for each of the cutoffs
#   x <- ggplot(data.sum.mut, aes(y=fitness, x=variants, group=sample, colour=reproduction, shape=geneFlow)) + geom_point(size=2)  + facet_grid(cutoff~gradient) + 
#     xlab("Number of variants") + ylab("Magnitude of fitness change") + 
#     scale_colour_manual(breaks = c("Asexual", "Sexual"), values=c("gold3", "blue4")) + 
#     theme_light() + theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(24),
#                           axis.title=element_text(size=16), strip.text = element_text(size=16)) 
#   x
#   
#   
#   #ggsave(x, path = "4_output/3_significance_Testing/", filename="SignificanceMutations_Fitness.png", device="png", height=250, width=220, dpi=320, units = "mm")
#   
#   
# }
# 
# #9) Do binomial/chisq test to assess if percentage coding differs from expectation
# {
#   #Determine percentage genes from annotation file
#   {
#     dat.ref <- read_tsv("2_data/TetThermRefGenome/GCF_000189635.1_JCVI-TTA1-2.2_genomic_genesOnly.gtf", col_names = F)
#     lengths <- sum(dat.ref$X5-dat.ref$X4)
#     
#     #Calculated from Eisen et al.: average length of gene (1815.4 bp) * number of genes (27424) divided by genome length (104194423 bp)
#     prop.coding <- 1815.4*27424/104194423
#     
#   }
#   
#   #Do binomial test for all data
#   {
#     coding.all <- nrow(filter(Mutations.all.func, expressed=="expressed"))
#     total.all <- nrow(Mutations.all.func)
#     binom.test(coding.all, total.all, p = prop.coding)
#   }
#   
#   #Do binomial test for indels only
#   {
#     indels <- filter(Mutations.all.func, TYPE=="INDEL")
#     coding.indels <- nrow(filter(indels, expressed=="expressed"))
#     binom.test(coding.indels, nrow(indels), p=prop.coding)
#   }
#   
#   #Do binomial test for SNP only
#   {
#     snp <- filter(Mutations.all.func, TYPE!="INDEL")
#     coding.snp <- nrow(filter(snp, expressed=="expressed"))
#     binom.test(coding.snp, nrow(snp), p=prop.coding)
#   }
#   
#   #Compare prevalence in coding regions between snp and indel
#   {
#     snpcoding <- nrow(filter(snp, expressed=="expressed"))
#     snpnoncoding <- nrow(filter(snp, expressed!="expressed"))
#     indelcoding <- nrow(filter(indels, expressed=="expressed"))
#     indelnoncoding <- nrow(filter(indels, expressed!="expressed"))
#     cont.table <- matrix(data=c(snpcoding, snpnoncoding, indelcoding, indelnoncoding), nrow=2, byrow = T)
#     x <- chisq.test(cont.table)
#     x
#   }
# }