#1) Clear memory
{
  rm(list=ls())
}

#2) Load packages
{
  library(tidyverse)
  library(ggpubr)
  '%!in%' <- function(x,y)!('%in%'(x,y))
}

#3) Set working directory
{
  setwd("/media/felix/BackupPlus/PhD/12_RangeExpansionGenomics")
  #setwd("/media/felix/Elements/12_RangeExpansionGenomics")
}

#4) read in data
{
  timepoints <- c( "t1", "t2", "t4", "t5", "t6", "t7", "t8", "t9", "t10", "t11", "t12")
  load("2_data/Old/5_EvolvedStrains/t0/Population_Data.RData")
  dd <- rbind(pop_output, pop_output, pop_output, pop_output, pop_output, pop_output, pop_output, pop_output)
  dd$hours <- 0
  dd$timepoint <- 't0'
  dd$testpH <- c(rep(6.5, 49), rep(6.0, 49), rep(5.5, 49), rep(5, 49), rep(4.5, 49), rep(4, 49), rep(3.5, 49), rep(3, 49))
  for (p in timepoints){
    path <- paste("2_data/Old/5_EvolvedStrains/", p, "/Population_Data.RData", sep="")
    load(path)
    if (p=="t5"){
      pop_output$time <- pop_output$hours
    }
    if (p=="t1"){
      pop_output$time <- c(rep(8, 36), rep(9, 46), 10, rep(9, 40), rep(10, 73), rep(11, 36),rep(12, 13), rep(11, 49), rep(9, 12), rep(10, 25), rep(11, 12), rep(8, 25), rep(9, 24))
    }
    pop_output$timepoint <- p
    dd <- rbind(dd, pop_output)
  }
  
  #Read in summarized data
  load("2_data/Old/6_EvolvedStrainsPosteriors/SummarisedPosteriorsEvolved.RData")
  sumdata <-sumoutput
}

#5) Data managing
{
  #5.1) Define hours variable
  dd$curveID <- paste(dd$ID, dd$testpH)
  dd$days<-  as.integer(as.Date(dd$date, "%d.%m.%y"))-17868
  dd$hours <- dd$days*24+dd$time - ifelse(dd$ID<23, 17, 18)
  
  #5.2) create indiv_per_ml variable
  dd$indiv_per_ml <- dd$indiv_per_volume*dd$dilution*1000
  
  #5.3) Remove missing datapoints
  dd <- dd[complete.cases(dd), ]
  dd$pHfact <- as.character(dd$testpH)
  
  # #Filter data to only have pH 6.5 and 4
  # dd <- filter(dd, testpH %in% c(6.5, 4))
  
  
  # #Work only with evolution block
  # dd <- filter(dd, strain=="mix")
  
  #Create identifier
  dd$ident <- paste(dd$ID, dd$testpH)
  dd$Sex <- ifelse(dd$Sex=="y", "Sex", "No Sex")
  dd$Gradient <- ifelse(dd$Gradient=="y", "Gradient", "No Gradient")
  dd$Gene.Flow <- ifelse(dd$Gene.Flow=="y", "Gene flow", "No Gene Flow")
  
  
  #calculate t0, inflection point and K point
  #Make dataset for t1
  dd.t0 <- filter(dd, timepoint=="t1")
  
  #Calculate mid log phase and K for all populations
  {
    sumdata <- drop_na(sumdata)
    dd.K <- data.frame()
    dd.inf <- data.frame()
    
    
    row.names(sumdata) <- sumdata$ident
    for (i in sumdata$ident){
      tempdata <- filter(dd, ident==i)
      timeK <- min(filter(tempdata, indiv_per_ml>=0.99*exp(sumdata[i, ]$logK.mean))$hours)
      timeK <- ifelse(timeK==Inf, tempdata[which.max(tempdata$indiv_per_ml), ]$hours, timeK)
      if(timeK> 0){
        timeInflect <- tempdata[which.min((tempdata[which(tempdata$hours<timeK), ]$indiv_per_ml - exp(sumdata[i, ]$logK.mean)/2)^2), ]$hours
        dd.inf <- rbind(dd.inf, filter(tempdata, hours==timeInflect))
      }
      dd.K <- rbind(dd.K, filter(tempdata, hours==timeK))
    }
  }
  
  #3.6) remove datapoints where no good datapoint was available (i.e. where population crashed and never reached K or mid log phase)
  dd.K <- na.omit(dd.K)
  dd.inf <- na.omit(dd.inf)
  
  #Filter for only evolution data
  dd.t0 <- filter(dd.t0, strain=="mix")
  dd.K <- filter(dd.K, strain=="mix")
  dd.inf <- filter(dd.inf, strain=="mix")
  
  #Filter for test pH
  dd.t0 <- filter(dd.t0, testpH %in% c(4, 6.5))
  dd.K <- filter(dd.K, testpH %in% c(4, 6.5))
  dd.inf <- filter(dd.inf, testpH %in% c(4, 6.5))
  
  #Collect in one nice dataset
  dd.t0$tp <- "t0"
  dd.K$tp <- "K"
  dd.inf$tp <- "inflection"
  dd.all <- rbind(dd.t0, dd.K, dd.inf)
  
  #Add info on r0 to the dataset
  dd.all$r0 <- 0
  for (i in unique(dd.all$ident)){
    dd.all[which(dd.all$ident==i), "r0"] <- sumoutput[which(sumoutput$ident==i), "logr0.mean"]
  }
}

#7) Calculate correlation and plot link in a) different environments and b) at different points in the growth curve
{
  #Rename and order facet names
  dd.all$pH <- ifelse(dd.all$testpH==4, "pH 4.0", "pH 6.5")
  dd.all$curvepoint <- ifelse(dd.all$tp=="t0", "Low density", ifelse(dd.all$tp=="inflection", "Medium density", "High density"))
  dd.all$curvepoint <- factor(dd.all$curvepoint, levels=c("Low density", "Medium density", "High density"), labels = c("Low density", "Medium density", "High density"))
  
  #First quick plot
  ggplot(dd.all, aes(x=r0, y=gross_speed_mean)) + geom_point() + facet_grid(curvepoint~pH)
  
  #Plot with correlation plotted for movement speed
  x <- ggscatter(dd.all, x="r0", y="gross_speed_mean", add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "Intrinsic rate of increase r0 (1/h; log-scale)", ylab = "Gross movement speed (μm/s)") + facet_grid(curvepoint~pH) + 
    theme_light()+ theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(24),
                         axis.title=element_text(size=16), strip.text = element_text(size=16), legend.position = "bottom") 
  x
  ggsave(x, path = "4_output/3_significance_Testing/", filename="r0_movement_corr.png", device = "png", height=150, width=200, dpi=320, units = "mm")
  
  
  
  #Plot with correlation plotted for turning speed
  y <- ggscatter(dd.all, x="r0", y="sd_turning_mean", add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "Intrinsic rate of increase r0 (1/h; log-scale)", ylab = "Turning speed (rad/s)") + facet_grid(curvepoint~pH) + 
    theme_light()+ theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(24),
                         axis.title=element_text(size=16), strip.text = element_text(size=16), legend.position = "bottom")
  y
  ggsave(y, path = "4_output/3_significance_Testing/", filename="r0_turning_corr.png", device = "png", height=150, width=200, dpi=320, units = "mm")
  
}