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
  library(lme4)
  '%!in%' <- function(x,y)!('%in%'(x,y))
}

#3) Set working directory
{
  setwd("/media/felix/Data/Felix_Data/Science/PhD/12_RangeExpansionGenomics")
  #setwd("/media/felix/Elements/12_RangeExpansionGenomics")
}

#4) Load data on AFC, and filter to keep significant sites with AFC>0.3 only (exclude clonal selection)
{
  load("2_data/R_DataFiles/AFC_Analysis/EvolvedData_AFC_func.RData")
  
  
  # #Compare similarity between 2 filtering methods (sign + effect>0.3 or significant under assumption that expected change=0.27)
  # {
  #   #prepare data 1
  #   dat1 <- filter(data.evolved.func, AFC_abs>0.3, p_BinomProb>10)
  # 
  #   #prepare data 2
  #   {
  #     dat2 <- data.evolved.func
  #     
  #     #Keep only entries that change by more than 0.27
  #     dat2 <- filter(dat2, AFC_abs>0.27)
  #     
  #     dat2$expectedUp <- ifelse(dat2$anc.ref.freq + 0.27>=1, 1, dat2$anc.ref.freq+0.27)
  #     dat2$expecteddown <- ifelse(dat2$anc.ref.freq - 0.27<=0, 0, dat2$anc.ref.freq-0.27)
  # 
  #     #redo binomial test based on expected change
  #     doBinomUpDown <- function(i){
  #       # ifelse(dat2[i, "anc.ref.freq"]>dat2[i, "ref.count"]/dat2[i, "total.count"],
  #       #        binom.test(dat2[i, "ref.count"], dat2[i, "total.count"], dat2[i, "expectedUp"], alternative = "greater")$p.value,
  #       #        binom.test(dat2[i, "ref.count"], dat2[i, "total.count"], dat2[i, "expecteddown"], alternative = "less")$p.value)
  #       up <- binom.test(dat2[i, "ref.count"], dat2[i, "total.count"], dat2[i, "expectedUp"], alternative = "greater")$p.value
  #       down <- binom.test(dat2[i, "ref.count"], dat2[i, "total.count"], dat2[i, "expecteddown"], alternative = "less")$p.value
  #       return(min(up, down))
  #     }
  #     dat2 <- as.data.frame(dat2)
  #     dat2$pval_binom2 <- unlist(lapply(1:nrow(dat2), doBinomUpDown))
  #     dat2$p_BinomProb2 <- -log10(dat2$pval_binom2)
  # 
  #     #Filter dat2 to only keep significant values
  #     bfcorr <- -log10(0.001/nrow(dat2))
  #     dat3 <- filter(dat2, p_BinomProb2 > bfcorr)
  # 
  #   }

    # #Compare similarity between dat1 and dat3
    # {
    #   common <- inner_join(dat1, dat3)
    #   rows.shared <- nrow(common)
    #   rows.unique1 <- nrow(dat1)-nrow(common)
    #   rows.unique3 <- nrow(dat3) - nrow(common)
    # }
    # 
    # #Rename p-value variables and save this output
    # {
    #   dat3$p_BinomProb <- dat3$p_BinomProb2
    #   dat3$pval_binom <- dat3$pval_binom2
    #   save(dat3, file = "2_data/R_DataFiles/AFC_Analysis/EvolvedData_AFC_func2.RData")
    # }
  # }

  #Reload data with more restrictive filtering
  
   data.evolved.func.filt <- filter(data.evolved.func, AFC_abs>0.3, p_BinomProb>10)
   data.evolved.func.sign <-  filter(data.evolved.func, p_BinomProb>10)
  
  #Plot histogram of unfiltered data
  y <- ggplot(data.evolved.func, aes(x=AFC_abs)) + geom_histogram() + facet_wrap(~sample) + xlab("Absolute change in allele frequency") + ylab("Count")+ xlim(-0.05, 1.05) +
    theme_light() + theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(24),
                          axis.title=element_text(size=16), strip.text = element_text(size=16), axis.text.x = element_text(angle = 90)) 
  
  ggsave(y, path = "4_output/3_significance_Testing/", filename="Hist_unfiltered.png", device="png", height=250, width=350, dpi=320, units = "mm")
  
  y <- ggplot(data.evolved.func.filt, aes(x=AFC_abs)) + geom_histogram() + facet_wrap(~sample) + xlab("Absolute change in allele frequency") + ylab("Count")+ xlim(-0.05, 1.05) +
    theme_light() + theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(24),
                          axis.title=element_text(size=16), strip.text = element_text(size=16), axis.text.x = element_text(angle = 90)) 
  
  ggsave(y, path = "4_output/3_significance_Testing/", filename="Hist_filtered.png", device="png", height=250, width=350, dpi=320, units = "mm")
}

#5) Compare significance between treatments
{
  load("2_data/R_DataFiles/AFC_Analysis/EvolvedData_AFC_func2.RData")
  # data.evolved.func.filt <- dat3
  # data.evolved.func.filt <- dat1
  
  #5.1) Look at number of variants, depending on different cutoffs for AFC
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
    data.sum.AFC$'Gene flow' <- data.sum.AFC$geneFlow
    data.sum.AFC$Reproduction <- data.sum.AFC$reproduction
    data.sum.AFC$`Gene flow`
    x <- ggplot(data.sum.AFC, aes(x=cutoff, y=(variants), group=sample, colour=Reproduction, shape=`Gene flow`)) + geom_point(size=2) + geom_line() + facet_wrap(~gradient) + 
      ylab("Number of variants (square root transformed)") + xlab("Magnitude of allele frequency change (cutoff)") + 
      scale_colour_manual(breaks = c("Asexual", "Sexual"), values=c("gold3", "blue4")) + 
      theme_light() + theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(24),
                            axis.title=element_text(size=16), strip.text = element_text(size=16)) 
    x
    ggsave(x, path = "4_output/3_significance_Testing/", filename="SignificanceAFC.png", device="png", height=150, width=250, dpi=320, units = "mm")
    
    #Plot data for all significant datapoints
    ggplot(filter(data.evolved.func), aes(x=AFC_raw, group=sample)) + geom_histogram() + facet_wrap(~sample)
    ggplot(filter(data.evolved.func), aes(x=AFC_abs, group=sample)) + geom_histogram() + facet_wrap(~sample)
    ggplot(filter(data.evolved.func.sign, AFC_abs>=0.3), aes(x=AFC_raw, group=sample)) + geom_histogram()
    ggplot(filter(data.evolved.func.sign, (AFC_raw <0.3 |  AFC_raw > 0)), aes(x=AFC_raw, group=sample)) + geom_histogram()
    
    #Do a mixed model with the data, separate for gradient/uniform
    model.sign <- lme(data=data.sum.AFC, (variants)~reproduction*geneFlow*gradient*(cutoff), random= ~1|sample, method="ML")
    mod_glmer1<-glmer(data=data.sum.AFC, variants~reproduction*geneFlow*gradient*cutoff+(1|sample) ,family="poisson", na.action="na.fail")
    #   plot(mod_glmer1)
     dredge(mod_glmer1, rank=BIC)
    best.model <-  glmer(data=data.sum.AFC, variants~reproduction*geneFlow*gradient*cutoff+(1|sample) ,family="poisson", na.action="na.fail")
    #best.model <-   lme(data=data.sum.AFC, (variants)~reproduction*(cutoff), random= ~1|sample, method="REML")
    plot(best.model)
    summary(best.model)
    Anova(best.model, contrasts=list(topic=contr.sum, sys=contr.sum), type=3)
    
    #Create prediction data
    pred.data <- expand.grid(geneFlow=c("Absent", "Present"), reproduction=c("Asexual", "Sexual"), gradient=c("Gradient", "Uniform"), cutoff=seq(from=0.3, to=0.8, by=0.01))
    predictions <- predict(best.model,re.form=NA,newdata=pred.data)
    #predictions <- predict(best.model, newdata=pred.data, level=0, se.fit=T) 
    # pred.data$mean <- predictions$fit
    # pred.data$upper <-  predictions$fit + predictions$se.fit*1.96
    # pred.data$lower <-  predictions$fit - predictions$se.fit*1.96
    
    
    library("sjPlot")
    pred.data$mean=exp(predictions)

    ??sjp.glmer
    se <- sqrt(diag(vcov(best.model)))
    tab <- cbind(Est = fixef(best.model), LL = fixef(best.model) - 1.96 * se, UL = fixef(best.model) + 1.96 * se)

    devtools::install_github("remkoduursma/bootpredictlme4")
    library(bootpredictlme4)
    pred2 <- predict(best.model, newdata=pred.data, re.form=NA, se.fit=TRUE, nsim=100)

    pred.data$mean <- exp(pred2$fit)
    pred.data$upper <-  exp(pred2$fit + 1.96*pred2$se.fit)
    pred.data$lower <-  exp(pred2$fit - 1.96*pred2$se.fit)
    pred.data$'Gene flow' <- pred.data$geneFlow
    pred.data$Reproduction <- pred.data$reproduction
   
    #Order facets correctly
    pred.data$`Abiotic conditions` <- ifelse(pred.data$gradient=="Gradient", "pH-gradient", "Uniform pH")
    pred.data$`Abiotic conditions` <- factor(pred.data$`Abiotic conditions` , levels=c("Uniform pH", "pH-gradient"))
    data.sum.AFC$`Abiotic conditions` <- ifelse(data.sum.AFC$gradient=="Gradient", "pH-gradient", "Uniform pH")
    data.sum.AFC$`Abiotic conditions` <- factor(data.sum.AFC$`Abiotic conditions`, levels=c("Uniform pH", "pH-gradient"))
    pred.data$gf <- paste("Gene flow:", pred.data$`Gene flow`)
    data.sum.AFC$gf <- paste("Gene flow:", data.sum.AFC$`Gene flow`)
    
    
    #Plot data + predictions
    x <- ggplot(data.sum.AFC, aes(x=cutoff, y=(variants), group=sample, colour=Reproduction)) + geom_point(size=3, alpha=0.3) + 
      geom_line(alpha=0.3)+ facet_grid(gf~`Abiotic conditions`) + xlim(0.25, 0.85) +
      ylab("Number of variants") + xlab("Magnitude of allele frequency change (cut-off)") +  
      geom_line(inherit.aes = F, data=pred.data, aes(x=cutoff, y=(mean), group=paste(reproduction, geneFlow, gradient), colour=Reproduction), size=1) + 
      geom_ribbon(inherit.aes = F, data=filter(pred.data), aes(x=cutoff, ymax=upper, ymin=lower, group=paste(reproduction, geneFlow), fill=Reproduction, colour=NA), size=3, alpha=0.3) +
      geom_point(inherit.aes = F, data=filter(pred.data, cutoff==0.3), aes(x=cutoff, y=(mean), group=paste(reproduction, geneFlow, gradient), colour=Reproduction), size=3) +
      scale_colour_manual(breaks = c("Asexual", "Sexual"), values=c("gold3", "blue4")) + 
      scale_fill_manual(breaks = c("Asexual", "Sexual"), values=c("gold3", "blue4")) +
      theme_light() + theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(24),
                            axis.title=element_text(size=16), strip.text = element_text(size=16), legend.position = "bottom") 
    x
    
    x2 <- tag_facet2(x, 
                     x = -Inf, y = 1200, 
                     vjust = 0, hjust = -0,
                     open = "", close = "",
                     fontface = 4,
                     size = 5,
                     tag_pool = c("A)", "B)", "C)", "D)"))
    x2
    
    #temp <- group_by(dat3, sample) %>%  summarise(count=n(), meandep=mean(total.count), reproduction=unique(reproduction), geneFlow=unique(geneFlow), gradient=unique(gradient))
    #plot(temp$count~temp$meandep)
    #cor.test(temp$count, temp$meandep)
    ggsave(x2, path = "4_output/3_significance_Testing/", filename="SignificanceAFC_predict.png", device="png", height=150, width=250, dpi=320, units = "mm")
    
  #   #redo for expressed variants only
  #   data.sum.AFC.expressed <- filter(data.evolved.func.filt, expressed=="expressed") %>% group_by(gradient, reproduction, geneFlow, sample) %>% summarize(variants=n(), cutoff=0.3)
  #   cutoffs <- seq(from=0.4, to=0.9, by=0.1)
  #   for(i in cutoffs){
  #     temp <- filter(data.evolved.func.filt, AFC_abs>i, expressed=="expressed") %>% group_by(gradient, reproduction, geneFlow, sample) %>% summarize(variants=n(), cutoff=i)
  #     data.sum.AFC.expressed <- rbind(data.sum.AFC.expressed, temp)
  #   }
  #   model.sign.expr <- lme(data=data.sum.AFC.expressed, variants~reproduction*geneFlow*gradient*cutoff, random= ~1|sample, method="ML")
  #   dredge(model.sign.expr, rank=BIC)
  #   best.model.expr <-  lme(data=data.sum.AFC.expressed, variants~geneFlow*gradient*reproduction + cutoff*gradient*reproduction, random= ~1|sample, method="REML")
  #   summary(best.model.expr)
  #   Anova(best.model.expr, contrasts=list(topic=contr.sum, sys=contr.sum), type=3)
  #   
  #   #redo for non-expressed variants only
  #   data.sum.AFC.nonexpressed <- filter(data.evolved.func.filt, expressed!="expressed") %>% group_by(gradient, reproduction, geneFlow, sample) %>% summarize(variants=n(), cutoff=0.3)
  #   cutoffs <- seq(from=0.4, to=0.9, by=0.1)
  #   for(i in cutoffs){
  #     temp <- filter(data.evolved.func.filt, AFC_abs>i, expressed!="expressed") %>% group_by(gradient, reproduction, geneFlow, sample) %>% summarize(variants=n(), cutoff=i)
  #     data.sum.AFC.nonexpressed <- rbind(data.sum.AFC.nonexpressed, temp)
  #   }
  #   model.sign.nonexpr <- lme(data=data.sum.AFC.nonexpressed, variants~reproduction*geneFlow*gradient*cutoff, random= ~1|sample, method="ML")
  #   dredge(model.sign.nonexpr, rank=BIC)
  #   best.model.nonexpr <-  lme(data=data.sum.AFC.nonexpressed, variants~gradient*reproduction*cutoff, random= ~1|sample, method="REML")
  #   summary(best.model.nonexpr)
  #   Anova(best.model.nonexpr, contrasts=list(topic=contr.sum, sys=contr.sum), type=3)
  #   
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
  data.sum.AFC$fitness <- 0
  for(i in 1:nrow(data.sum.AFC)){
    data.sum.AFC[i, "fitness"] <- dd.evo[which(dd.evo$sample==as.character(data.sum.AFC[i, "sample"])), "logratio"]
  }
  data.sum.AFC$`Abiotic conditions` <- data.sum.AFC$`Abiotic conditions`
  data.sum.AFC$`Cut-off` <- data.sum.AFC$cutoff
  
  #Calculate for each of the panels the correlation between change in growth rate and change in allele frequecy
  {
    temp <- filter(data.sum.AFC, cutoff<0.8) %>% group_by(`Abiotic conditions`, Reproduction, cutoff) %>% summarise(cor= cor.test(fitness, variants)$estimate, pval = cor.test(fitness, variants)$p.value)
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
  data.sum.AFC$co2 <- paste("Cut-off:", data.sum.AFC$`Cut-off`)
  
  
  #Plot the fitness values for each of the cutoffs
  x2 <- ggplot(data.sum.AFC, aes(y=fitness, x=(variants), group=sample, colour=Reproduction, shape=`Gene flow`)) + geom_point(size=3)  + facet_grid(co2~`Abiotic conditions`) + 
    xlab("Number of variants (genotypic evolution)") + ylab("Magnitude of fitness change (phenotypic evolution)") +
    scale_colour_manual(breaks = c("Asexual", "Sexual"), values=c("gold3", "blue4")) + 
    theme_light() + theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(24),
                          axis.title=element_text(size=16), strip.text = element_text(size=16), legend.position = "bottom") + ylim(0, 1) + xlim(-50, 1200)
  
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
            x = 700, y = c(0.9, rep(0.25, 11)), 
            vjust = 0, hjust = -0,
            open = "", close = "",
            fontface = 4,
            size = 4,
            colour="khaki3",
            tag_pool = c(t1$text1, "", "", "", ""))
  
  x2 <- tag_facet2(x2, 
                  x = 700, y = c(0.75, rep(0.10, 11)), 
                  vjust = 0, hjust = -0,
                  open = "", close = "",
                  fontface = 4,
                  size = 4,
                  colour="khaki3",
                  tag_pool = c(t1$text3, "", "", "", ""))
  
  x2 <- tag_facet2(x2, 
                   x = 1000, y = c(0.9, rep(0.25, 11)), 
                  vjust = 0, hjust = -0,
                  open = "", close = "",
                  fontface = 4,
                  size = 4,
                  colour= c("cornflowerblue", "cornflowerblue", rep(c("blue4", "cornflowerblue"), 5)),
                  tag_pool = c(t3$text1, "", "", "", ""))
  
  x2 <- tag_facet2(x2, 
                   x = 1000, y = c(0.75, rep(0.10, 11)), 
                  vjust = 0, hjust = -0,
                  open = "", close = "",
                  fontface = 4,
                  size = 4,
                  colour=c("cornflowerblue", "cornflowerblue", rep(c("blue4", "cornflowerblue"), 5)),
                  tag_pool = c(t3$text3, "", "", "", ""))
  
  x2
  ggsave(x2, path = "4_output/3_significance_Testing/", filename="SignificanceAFC_Fitness.png", device="png", height=250, width=250, dpi=320, units = "mm")
 
  #Save the output data
  
  
  # 
  # dd <- filter(data.sum.AFC, cutoff==0.3)
  # test <- lm(data=dd, fitness ~ variants*reproduction*geneFlow*gradient, na.action="na.fail")
  # test2 <- lm(data=dd, fitness ~ variants, na.action="na.fail")
  # summary(test)
  # Anova(test, contrasts=list(topic=contr.sum, sys=contr.sum), type=3)
  # 
  # x2 <- predict(test)
  # plot(x2~dd$fitness)
  # test2 <- lm(data=dd, fitness ~ variants*reproduction, na.action="na.fail")
  # library(MuMIn)
  # dredge(test)

  #Compare fitness among treatments
  {
    #Calculate mean growth rate per pH value
    dd.evo
    
    #Do model selection
    fullmodel2 <- lm(data = dd.evo, logratio ~ Sex*Gene.Flow*Gradient, na.action = "na.fail")
    t2 <- dredge(fullmodel2)
    model.avg(t2)
    best.model2 <- lm(data = dd.evo, logratio ~ Sex + Gene.Flow + Gradient + Sex*Gene.Flow, na.action = "na.fail")
    anova(best.model2)
    summary(best.model2)
    Anova(best.model2, contrasts=list(topic=contr.sum, sys=contr.sum), type=3)
    
    #Plot the data
    #Create prediction data
    prd.data2 <- expand.grid(Sex = unique(dd.evo$Sex), Gene.Flow = unique(dd.evo$Gene.Flow), Gradient = unique(dd.evo$Gradient))
    predictions2 <- predict(best.model2, newdata = prd.data2, se.fit = T)
    prd.data2$meanr0 <-predictions2$fit
    prd.data2$upperr0 <- predictions2$fit + 1.96*predictions2$se.fit
    prd.data2$lowerr0 <- predictions2$fit - 1.96*predictions2$se.fit
    
    dd.evo$`Gene flow` <- dd.evo$Gene.Flow
    prd.data2$`Gene flow` <- prd.data2$Gene.Flow
    
    dd.evo$`Abiotic conditions` <- dd.evo$Gradient
    dd.evo$`Abiotic conditions` <- ifelse(dd.evo$`Abiotic conditions`=="Gradient", "pH-gradient", "Uniform pH")
    dd.evo$`Abiotic conditions` <- factor(dd.evo$`Abiotic conditions` , levels=c("Uniform pH", "pH-gradient"))
    prd.data2$`Abiotic conditions` <- prd.data2$Gradient
    prd.data2$`Abiotic conditions` <- ifelse(prd.data2$`Abiotic conditions`=="Gradient", "pH-gradient", "Uniform pH")
    prd.data2$`Abiotic conditions` <- factor(prd.data2$`Abiotic conditions` , levels=c("Uniform pH", "pH-gradient"))
    dd.evo$gf <- paste("Gene flow:", dd.evo$`Gene flow`)
    prd.data2$gf <- paste("Gene flow:", prd.data2$`Gene flow`)
    
    
    #Plot data + predictions
    y <- ggplot(dd.evo, aes(y=logratio, x=Sex, color = Sex)) + geom_point(size=2) + facet_grid(gf~`Abiotic conditions`) +
      geom_boxplot(inherit.aes = F, data = prd.data2, mapping = aes(fill = Sex, colour = NA, middle = meanr0, ymax = upperr0, ymin = lowerr0, upper = upperr0, lower = lowerr0, x = Sex), stat = "identity", alpha = 0.3) + 
      geom_boxplot(inherit.aes = F, data = prd.data2, mapping = aes(fill = NA, middle = meanr0, ymax = meanr0, ymin = meanr0, upper = meanr0, lower = meanr0, x = Sex), stat = "identity", alpha = 0.3) + 
      ylab("Magnitude of fitness change (genotypic evolution)") + 
      theme_light() + theme(axis.text=element_text(size=10), legend.text=element_text(size=14),legend.title=element_text(size=14), strip.text.x=element_text(24),
                            axis.title=element_text(size=14), strip.text = element_text(size=11), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                            axis.title.x = element_blank(), legend.position = "bottom") +
      scale_color_manual(values=c("gold3", "blue4"), breaks=c("Asexual", "Sexual"), name="Reproduction", labels=c("Asexual", "Sexual")) +
      scale_fill_manual(values=c("gold3", "blue4"), breaks=c("Asexual", "Sexual"), name="Reproduction", labels=c("Asexual", "Sexual"))
    
    x2 <- tag_facet2(y, 
                     x = -Inf, y = 0.9, 
                     vjust = 0, hjust = -0,
                     open = "", close = "",
                     fontface = 4,
                     size = 5,
                     tag_pool = c("A)", "B)", "C)", "D)"))
    x2
    
    ggsave(x2, path = "4_output/3_significance_Testing/", filename="Fitness.png", device = "png", width = 7, height = 6, dpi = 300)
    
  }
  
  #get some basic stats
  nrow(data.evolved.func.filt)
  data.evolved.func.filt %>% group_by(expressed) %>% summarize(count=n())
  3902/13278
  data.evolved.func.filt %>% group_by(TYPE) %>% summarize(count=n()/13278)
  
  sum <- dat1 %>% group_by(expressed, TYPE) %>% summarize(n())
  sum2 <- dat1 %>% group_by(expressed) %>% summarize(n())
  
  #Do binomial/chisq test to assess if percentage coding differs from expectation
  {
    #Determine percentage genes from annotation file
    {
      dat.ref <- read_tsv("2_data/TetThermRefGenome/GCF_000189635.1_JCVI-TTA1-2.2_genomic_genesOnly.gtf", col_names = F)
      lengths <- sum(dat.ref$X5-dat.ref$X4)
      
      #Calculated from Eisen et al.: average length of gene (1815.4 bp) * number of genes (27424) divided by genome length (104194423 bp)
      prop.coding <- 1815.4*27424/104194423
        
    }
    
    #Do binomial test for all data
    {
      coding.all <- nrow(filter(data.evolved.func.filt, expressed=="expressed"))
      total.all <- nrow(data.evolved.func.filt)
      binom.test(coding.all, total.all, p = prop.coding)
    }
    
    #Do binomial test for indels only
    {
      indels <- filter(data.evolved.func.filt, TYPE=="INDEL")
      coding.indels <- nrow(filter(indels, expressed=="expressed"))
      binom.test(coding.indels, nrow(indels), p=prop.coding)
    }
    
    #Do binomial test for SNP only
    {
      snp <- filter(data.evolved.func.filt, TYPE!="INDEL")
      coding.snp <- nrow(filter(snp, expressed=="expressed"))
      binom.test(coding.snp, nrow(snp), p=prop.coding)
    }
    
    #Compare prevalence in coding regions between snp and indel
    {
      snpcoding <- nrow(filter(snp, expressed=="expressed"))
      snpnoncoding <- nrow(filter(snp, expressed!="expressed"))
      indelcoding <- nrow(filter(indels, expressed=="expressed"))
      indelnoncoding <- nrow(filter(indels, expressed!="expressed"))
      cont.table <- matrix(data=c(snpcoding, snpnoncoding, indelcoding, indelnoncoding), nrow=2, byrow = T)
      x <- chisq.test(cont.table)
      x$method
    }
    
  }
  

}