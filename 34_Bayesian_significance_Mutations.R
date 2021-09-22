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

#4) Load data on mutations, and filter to keep onlly sites with change>0.1
{
  load("2_data/R_DataFiles/MutationAnalysis/AllMutations_Func.RData")
  Mutations.all.func.filt <- filter(Mutations.all.func, alt.freq>0.1)
  
  #Prepare data
  data.sum.mut <- group_by(Mutations.all.func.filt, gradient, reproduction, geneFlow, sample) %>% summarize(variants=n(), cutoff=0.1, perc.Indel=sum(TYPE == "INDEL")/n())
  cutoffs <- seq(from=0.2, to=0.9, by=0.1)
  for(i in cutoffs){
    temp <- filter(Mutations.all.func.filt, alt.freq>i) %>% group_by(gradient, reproduction, geneFlow, sample) %>% summarize(variants=n(), cutoff=i, perc.Indel=sum(TYPE == "INDEL")/n())
    data.sum.mut <- rbind(data.sum.mut, temp)
  }
}

#5) Make all possible models using the MCMCglmm package, and compare them afterwards based on the BIC criterion
{
  data.stan <- data.frame(select(data.sum.mut, gradient, reproduction, geneFlow, sample, variants, cutoff) %>%
                            mutate(gradient=ifelse(gradient=="Uniform", 0, 1), reproduction=ifelse(reproduction=="Asexual", 0, 1),
                                   geneFlow=ifelse(geneFlow=="Absent", 0, 1), sample=sample))
  
  m.1 <- map2stan(
    alist(
      variants ~ dnorm(mu, sigma),
      mu <- a + as[sample]+ a1*cutoff,
      a ~dnorm(1800, 500),
      as[sample] ~dnorm(0, 500),
      a1 ~dnorm(0, 500),
      sigma ~ dcauchy(0, 1)
    ) ,data=data.stan ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  summary(m.1)
  
  m.all <- map2stan(
    alist(
      variants ~ dnorm(mu, sigma),
      mu <- a + as[sample] +  a1*cutoff + a2*reproduction + a3*gradient + a4*geneFlow,
      a ~dnorm(1800, 500),
      as[sample] ~dnorm(0, 500),
      a1 ~dnorm(0, 500),
      a2 ~dnorm(0, 500),
      a3 ~dnorm(0, 500),
      a4 ~dnorm(0, 500),
      sigma ~ dcauchy(0, 1)
    ) ,data=data.stan ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  summary(m.all)
  
  m.int <- map2stan(
    alist(
      variants ~ dnorm(mu, sigma),
      mu <- a + as[sample] + a1*cutoff + a2*reproduction + a3*gradient + a4*geneFlow +
        a12*cutoff*reproduction + a13*cutoff*gradient + a14*cutoff*geneFlow +
        a23*reproduction*gradient + a24*reproduction*geneFlow + a34*gradient*geneFlow + 
        a123*cutoff*reproduction*gradient + a124*reproduction*gradient*geneFlow + a234*reproduction*gradient*geneFlow  + a134*cutoff*gradient*geneFlow + 
        a1234*cutoff*reproduction*gradient*geneFlow,
      a ~dnorm(1800, 500),
      as[sample] ~dnorm(0, 500),
      a1 ~dnorm(0, 500),
      a2 ~dnorm(0, 500),
      a3 ~dnorm(0, 500),
      a4 ~dnorm(0, 500),
      a12 ~dnorm(0, 500),
      a13 ~dnorm(0, 500),
      a14 ~dnorm(0, 500),
      a23 ~dnorm(0, 500),
      a24 ~dnorm(0, 500),
      a34 ~dnorm(0, 500),
      a123 ~dnorm(0, 500),
      a124 ~dnorm(0, 500),
      a134 ~dnorm(0, 500),
      a234 ~dnorm(0, 500),
      a1234 ~dnorm(0, 500),
      sigma ~dcauchy(0, 1)
    ) ,data=data.stan ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  summary(m.int)
  #plot(m.int)
  
  compare(m.1, m.all, m.int, WAIC = T)
  
  #Save the output from the models
  models.mut <- list(simple=m.1, intermediate=m.all, full=m.int)
  save(models.mut, file="2_data/R_DataFiles/Posteriors_MutationModels.RData")
  
  #Get output from best model
  precis(m.int, digits=4, depth=1, prob = 0.95)
  
}

#6) Plot raw data on top of posterior prediction (mean and 95 % CI)
{
  #Order facets correctly
  data.sum.mut$`Abiotic conditions` <- ifelse(data.sum.mut$gradient=="Gradient", "pH-gradient", "Uniform pH")
  data.sum.mut$`Abiotic conditions` <- factor(data.sum.mut$`Abiotic conditions`, levels=c("Uniform pH", "pH-gradient"))
  
  #Create prediction data
  {
    #Prepare data to calculate predictions
    pred.data <- expand.grid(reproduction=c(0, 1), cutoff=seq(from=0.1, to=0.9, by=0.01), gradient=c(0, 1), geneFlow=c(0, 1))
    
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
    data.sum.mut$gf <- paste("Gene flow:", data.sum.mut$geneFlow)
    data.sum.mut$Reproduction <- data.sum.mut$reproduction
    
    #Plot raw data and posterior predictions
    x <- ggplot(data.sum.mut, aes(x=cutoff, y=(variants), group=sample, colour=Reproduction)) + geom_point(size=3, alpha=0.3) + 
      geom_line(alpha=0.3)+ facet_grid(gf~`Abiotic conditions`) + xlim(0.05, 0.95) +
      ylab("Number of variants") + xlab("Allele frequency of novel variants (cut-off)") +  
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
                     x = -Inf, y = 2400, 
                     vjust = 0, hjust = -0,
                     open = "", close = "",
                     fontface = 4,
                     size = 5,
                     tag_pool = c("A)", "B)", "C)", "D)"))
    x2
  }
  
  #Save the figure 
  ggsave(x2, path = "4_output/3_significance_Testing/", filename="Fitness_mut_Bayesian.png", device = "png", height=150, width=250, dpi=320, units = "mm")
}
