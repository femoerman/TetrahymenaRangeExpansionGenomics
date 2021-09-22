#1) Clear memory
{
  rm(list=ls())
}

#2) Load packages
{
  # library(devtools)
  # install_github("rmcelreath/rethinking")
  library(MCMCglmm)
  library(tidyverse)
  library(MuMIn)
  library(rethinking)
  '%!in%' <- function(x,y)!('%in%'(x,y))
}

#3) Set working directory
{
  setwd("/media/felix/Data/Felix_Data/Science/PhD/12_RangeExpansionGenomics")
  #setwd("/media/felix/Elements/12_RangeExpansionGenomics")
}

#4) Load data on mutations and AFC, and filter to keep onlly sites with change>0.3. Calculate the sum of all changes
{
  #Load mutation data
  load("2_data/R_DataFiles/MutationAnalysis/AllMutations_Func.RData")
  Mutations.all.func.filt <- filter(Mutations.all.func, alt.freq>0.1)
  
  #Prepare data
  data.sum.mut <- group_by(Mutations.all.func.filt, gradient, reproduction, geneFlow, sample) %>% summarize(variants=n(), cutoff=0.1, perc.Indel=sum(TYPE == "INDEL")/n())
  cutoffs <- seq(from=0.2, to=0.9, by=0.1)
  for(i in cutoffs){
    temp <- filter(Mutations.all.func.filt, alt.freq>i) %>% group_by(gradient, reproduction, geneFlow, sample) %>% summarize(variants=n(), cutoff=i, perc.Indel=sum(TYPE == "INDEL")/n())
    data.sum.mut <- rbind(data.sum.mut, temp)
  }
  data.sum.mut <- filter(data.sum.mut, cutoff>=0.3)
  
  #Load data on AFC
  load("2_data/R_DataFiles/AFC_Analysis/EvolvedData_AFC_func.RData")
  data.evolved.func.filt <- filter(data.evolved.func, AFC_abs>0.3, p_BinomProb>10)
  
  #Summarize the dataset for the number of variants
  {
    #Prepare data
    data.sum.AFC <- group_by(data.evolved.func.filt, gradient, reproduction, geneFlow, sample) %>% summarize(variants=n(), cutoff=0.3)
    cutoffs <- seq(from=0.4, to=0.8, by=0.1)
    for(i in cutoffs){
      temp <- filter(data.evolved.func.filt, AFC_abs>i) %>% group_by(gradient, reproduction, geneFlow, sample) %>% summarize(variants=n(), cutoff=i)
      data.sum.AFC <- rbind(data.sum.AFC, temp)
    }
    
    #Plot data
    data.sum.AFC$geneFlow <- factor(data.sum.AFC$geneFlow, levels=c("Absent", "Present"))
    data.sum.AFC$reproduction <- factor(data.sum.AFC$reproduction, levels=c("Asexual", "Sexual"))
    data.sum.AFC$gradient <- factor(data.sum.AFC$gradient, levels=c("Uniform", "Gradient"))
  }
  
 #Create a summed value for all variants
  data.sum.AFC$ID <-paste0(data.sum.AFC$sample, data.sum.AFC$cutoff)
  data.sum.AFC <- data.frame(mutate(data.sum.AFC, AFC=variants)) %>% dplyr::select(-variants, -gradient, -reproduction, -geneFlow, -sample, -cutoff)
  data.sum.mut$ID <-paste0(data.sum.mut$sample, data.sum.mut$cutoff)
  data.sum.mut <- mutate(data.sum.mut, mut=variants) %>% select(-variants, -perc.Indel)
  
  #Merge the two and add cleaned treatment names
  data.sum.both <- merge(data.sum.AFC, data.sum.mut, by="ID", all = T)
  data.sum.both$AFC <- replace_na(data.sum.both$AFC, 0)
  data.sum.both$variants <- data.sum.both$AFC + data.sum.both$mut
  data.sum.both$'Gene flow' <- data.sum.both$geneFlow
  data.sum.both$Reproduction <- data.sum.both$reproduction
}

#5) Make all possible models using the MCMCglmm package, and compare them afterwards based on the BIC criterion
{
  data.stan <- data.frame(select(data.sum.both, gradient, reproduction, geneFlow, sample, variants, cutoff) %>%
                            mutate(gradient=ifelse(gradient=="Uniform", 0, 1), reproduction=ifelse(reproduction=="Asexual", 0, 1),
                                   geneFlow=ifelse(geneFlow=="Absent", 0, 1), sample=sample))
  
  
  m.1 <- map2stan(
    alist(
      variants ~ dpois(mu),
      log(mu) <- a + as[sample]+ a1*cutoff,
      a ~dnorm(6,2.5),
      as[sample] ~dnorm(0, 0.2),
      a1 ~dnorm(0, 2)
    ) ,data=data.stan ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  summary(m.1)
  
  m.all <- map2stan(
    alist(
      variants ~ dpois(mu),
      log(mu) <- a + as[sample] +  a1*cutoff + a2*reproduction + a3*gradient + a4*geneFlow,
      a ~dnorm(6,2.5),
      as[sample] ~dnorm(0, 0.2),
      a1 ~dnorm(0, 2),
      a2 ~dnorm(0, 2),
      a3 ~dnorm(0, 2),
      a4 ~dnorm(0, 2)
    ) ,data=data.stan ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  summary(m.all)
  
  m.int <- map2stan(
    alist(
      variants ~ dpois(mu),
      log(mu) <- a + as[sample] + a1*cutoff + a2*reproduction + a3*gradient + a4*geneFlow +
        a12*cutoff*reproduction + a13*cutoff*gradient + a14*cutoff*geneFlow +
        a23*reproduction*gradient + a24*reproduction*geneFlow + a34*gradient*geneFlow + 
        a123*cutoff*reproduction*gradient + a124*reproduction*gradient*geneFlow + a234*reproduction*gradient*geneFlow + a134*cutoff*gradient*geneFlow + 
        a1234*cutoff*reproduction*gradient*geneFlow,
      a ~dnorm(6,2.5),
      as[sample] ~dnorm(0, 0.2),
      a1 ~dnorm(0, 2),
      a2 ~dnorm(0, 2),
      a3 ~dnorm(0, 2),
      a4 ~dnorm(0, 2),
      a12 ~dnorm(0, 2),
      a13 ~dnorm(0, 2),
      a14 ~dnorm(0, 2),
      a23 ~dnorm(0, 2),
      a24 ~dnorm(0, 2),
      a34 ~dnorm(0, 2),
      a123 ~dnorm(0, 2),
      a124 ~dnorm(0, 2),
      a134 ~dnorm(0, 2),
      a234 ~dnorm(0, 2),
      a1234 ~dnorm(0, 2)
    ) ,data=data.stan ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  summary(m.int)
  plot(m.int)
  
  compare(m.1, m.all, m.int, WAIC = T)
  
  #Save the output from the models
  models.both <- list(simple=m.1, intermediate=m.all, full=m.int)
  save(models.both, file="2_data/R_DataFiles/Posteriors_BothModels.RData")
  
  #Get output from best model
  precis(m.int, digits=4, depth=1, prob = 0.95)
}

#6) Plot raw data on top of posterior prediction (mean and 95 % CI)
{
  #Order facets correctly
  data.sum.both$`Abiotic conditions` <- ifelse(data.sum.both$gradient=="Gradient", "pH-gradient", "Uniform pH")
  data.sum.both$`Abiotic conditions` <- factor(data.sum.both$`Abiotic conditions`, levels=c("Uniform pH", "pH-gradient"))
  
  #Create prediction data
  {
    #Prepare data to calculate predictions
    pred.data <- expand.grid(reproduction=c(0, 1), cutoff=seq(from=0.3, to=0.9, by=0.1), gradient=c(0, 1), geneFlow=c(0, 1))
    
    #Extract samples from the best model
    post <- extract.samples(m.int, n=8e3)
    post.pred=link(m.int, pred.data, n=8e3)
    mu.pred.mean <- apply( post.pred , 2 , mean )
    mu.pred.PI <- apply(post.pred, 2, PI, prob=0.95)
  
    #Create new variables with original character variables for reproduction, gene flow and gradient
    posteriorpred <- cbind(pred.data, data.frame(mean=(mu.pred.mean), upper=(mu.pred.PI[2, ]), lower=(mu.pred.PI[1, ])))
    posteriorpred <- mutate(posteriorpred, reproduction=ifelse(reproduction==1, "Sexual", "Asexual"), geneFlow=ifelse(geneFlow==1, "Present", "Absent"), 
                            gradient=ifelse(gradient==1, "Gradient", "Uniform"))
    posteriorpred$Reproduction <- posteriorpred$reproduction
    posteriorpred$'Gene flow' <- posteriorpred$geneFlow
    posteriorpred$`Abiotic conditions` <- ifelse(posteriorpred$gradient=="Gradient", "pH-gradient", "Uniform pH")
    posteriorpred$`Abiotic conditions` <- factor(posteriorpred$`Abiotic conditions`, levels=c("Uniform pH", "pH-gradient"))
    posteriorpred$gf <- paste("Gene flow:", posteriorpred$`Gene flow`)
    data.sum.both$gf <- paste("Gene flow:", data.sum.both$`Gene flow`)
    
    #Plot raw data and posterior predictions
    x <- ggplot(data.sum.both, aes(x=cutoff, y=(variants), group=sample, colour=Reproduction)) + geom_point(size=3, alpha=0.3) + 
      geom_line(alpha=0.3)+ facet_grid(gf~`Abiotic conditions`) + xlim(0.25, 0.95) +
      ylab("Number of variants") + xlab("Magnitude of genetic change (cut-off)") +  
      geom_line(inherit.aes = F, data=posteriorpred, aes(x=cutoff, y=(mean), group=paste(reproduction, geneFlow, gradient), colour=Reproduction), size=1) + 
      geom_ribbon(inherit.aes = F, data=filter(posteriorpred), aes(x=cutoff, ymax=upper, ymin=lower, group=paste(reproduction, geneFlow), fill=Reproduction, colour=NA), size=3, alpha=0.3) +
      geom_point(inherit.aes = F, data=filter(posteriorpred, cutoff==0.3), aes(x=cutoff, y=(mean), group=paste(reproduction, geneFlow, gradient), colour=Reproduction), size=3) +
      scale_colour_manual(breaks = c("Asexual", "Sexual"), values=c("gold3", "blue4")) + 
      scale_fill_manual(breaks = c("Asexual", "Sexual"), values=c("gold3", "blue4")) +
      theme_light() + theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(24),
                            axis.title=element_text(size=16), strip.text = element_text(size=16), legend.position = "bottom") 
    
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
                     x = -Inf, y = 3500, 
                     vjust = 0, hjust = -0,
                     open = "", close = "",
                     fontface = 4,
                     size = 5,
                     tag_pool = c("A)", "B)", "C)", "D)"))
    x2
  }
  
  #Save the figure 
  ggsave(x2, path = "4_output/3_significance_Testing/", filename="Fitness_Both_Bayesian.png", device = "png", height=150, width=250, dpi=320, units = "mm")
}

#7) Assess correlation with fitness change
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
  data.sum.both$fitness <- 0
  for(i in 1:nrow(data.sum.both)){
    data.sum.both[i, "fitness"] <- dd.evo[which(dd.evo$sample==as.character(data.sum.both[i, "sample"])), "logratio"]
  }
  data.sum.both$`Abiotic conditions` <- data.sum.both$`Abiotic conditions`
  data.sum.both$`Cut-off` <- data.sum.both$cutoff
  
  #Calculate for each of the panels the correlation between change in growth rate and change in allele frequecy
  {
    temp <- filter(data.sum.both, cutoff<0.8) %>% group_by(`Abiotic conditions`, Reproduction, cutoff) %>% summarise(cor= cor.test(fitness, variants)$estimate, pval = cor.test(fitness, variants)$p.value)
    temp$r2 <- temp$cor^2
    
    temp$text1 <- paste("p-value =", round(temp$pval, digits=4))
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
  data.sum.both$co2 <- paste("Cut-off:", data.sum.both$`Cut-off`)
  
  #Plot the fitness values for each of the cutoffs
  x2 <- ggplot(data.sum.both, aes(y=fitness, x=(variants), group=sample, colour=Reproduction, shape=`Gene flow`)) + geom_point(size=3)  + facet_grid(co2~`Abiotic conditions`) + 
    xlab("Number of variants (molecular evolution)") + ylab("Magnitude of fitness change (phenotypic evolution)") +
    scale_colour_manual(breaks = c("Asexual", "Sexual"), values=c("gold3", "blue4")) + 
    theme_light() + theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(24),
                          axis.title=element_text(size=16), strip.text = element_text(size=16), legend.position = "bottom") + ylim(-0.5, 1)
  
  x2
  
  t3 <- arrange(temp, Reproduction, `cutoff`) %>% filter(Reproduction=="Sexual")
  t1 <- arrange(temp, Reproduction, `cutoff`) %>% filter(Reproduction=="Asexual")
  
  x2 <- tag_facet2(x2, 
                   x = -Inf, y = 0.9, 
                   vjust = 0, hjust = -0,
                   open = "", close = "",
                   fontface = 4,
                   size = 5,
                   tag_pool = c("A)", "B)", "C)", "D)", "E)", "F)", "G)", "H)", "I)", "J)", "K)", "L)"))
  
  
  x2 <- tag_facet2(x2, 
                   x = 10, y = -0.25, 
                   vjust = 0, hjust = -0,
                   open = "", close = "",
                   fontface = 4,
                   size = 3,
                   colour="gold3",
                   tag_pool = c(t1$text1, "", "", "", ""))
  
  x2 <- tag_facet2(x2, 
                   x = 10, y = -0.5, 
                   vjust = 0, hjust = -0,
                   open = "", close = "",
                   fontface = 4,
                   size = 3,
                   colour="gold3",
                   tag_pool = c(t1$text3, "", "", "", ""))
  
  x2 <- tag_facet2(x2, 
                   x = 300, y = -0.25, 
                   vjust = 0, hjust = -0,
                   open = "", close = "",
                   fontface = 4,
                   size = 3,
                   colour="blue4",
                   tag_pool = c(t3$text1, "", "", "", ""))
  
  x2 <- tag_facet2(x2, 
                   x = 300, y = -0.5, 
                   vjust = 0, hjust = -0,
                   open = "", close = "",
                   fontface = 4,
                   size = 3,
                   colour="blue4",
                   tag_pool = c(t3$text3, "", "", "", ""))
  
  x2
  }