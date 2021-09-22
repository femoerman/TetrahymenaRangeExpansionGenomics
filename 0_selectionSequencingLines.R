#1) Clear memory
rm(list=ls())

#2) Load packages
{
  library(tidyverse)
}

#3) Set working directory and load in data
{
  #setwd("/media/felix/DataDrive2/Documenten/PhD/12_RangeExpansionGenomics")
  setwd("F:/Documenten/PhD/12_RangeExpansionGenomics")
  load("2_data/5_EvolvedStrains/1_RawDensityData.RData")
  load("2_data/6_EvolvedStrainsPosteriors/SummarisedPosteriorsEvolved.RData")
  evolved <- sumoutput
  load("2_data/summarizedAncestorData.RData")
  ancestor <- sumdata
}

#4) Combine posterior data with raw data and un-log transform posteriors
{
  dd <- dd %>% filter(hours==0) %>% arrange(curveID) %>% select(curveID, strain, Treatment, testpH, Gradient, Sex, Gene.Flow, ID)
  plotdata <- cbind(dd, sumoutput)
  
  plotdata$r0.mean <- ifelse(is.na(plotdata$logr0.mean), 0, exp(plotdata$logr0.mean))
  plotdata$K.mean <- ifelse(is.na(plotdata$logK.mean), NA, exp(plotdata$logK.mean))
  plotdata$alpha.mean <- ifelse(is.na(plotdata$logalpha.mean), NA, exp(plotdata$logalpha.mean))
  plotdata$d.mean <- ifelse(is.na(plotdata$logd.mean), NA, exp(plotdata$logd.mean))
  plotdata$Sex <- ifelse(plotdata$Sex=="y", "Sex", "No Sex")
  plotdata$Gradient <- ifelse(plotdata$Gradient=="y", "Gradient", "No Gradient")
  plotdata$Gene.Flow <- ifelse(plotdata$Gene.Flow=="y", "Gene flow", "No Gene Flow")

  sumdata$testpH <- paste(sumdata$pH, "j")
}

#Set bad fits to 0?
sumdata$r0mean <- ifelse(is.na(sumdata$r0mean), 0, sumdata$r0mean)
sumdata <- mutate(sumdata, r0mean = ifelse(pH==3.5 & strain != "B2086.2", 0, r0mean))
sumdata <- mutate(sumdata, r0mean = ifelse(pH==4 & strain == "CU427.4", 0, r0mean))


#5) Isolate highest r0 value lines for gradient treatment at pH 4
{
  #Filter data to get evolution lines with gradient
  gradient <- filter(plotdata, Gradient == "Gradient" & Treatment == "EVO")
  gradient <- filter(gradient, testpH==4)
  gr.sex.geneflow <- filter(gradient, Sex=="Sex" & Gene.Flow == "Gene flow") %>% arrange(-logr0.mean) %>% slice(1:2)
  gr.nosex.geneflow <- filter(gradient, Sex!="Sex" & Gene.Flow == "Gene flow") %>% arrange(-logr0.mean) %>% slice(1:2)
  gr.sex.nogeneflow <- filter(gradient, Sex=="Sex" & Gene.Flow != "Gene flow") %>% arrange(-logr0.mean) %>% slice(1:2)
  gr.nosex.nogeneflow <- filter(gradient, Sex!="Sex" & Gene.Flow != "Gene flow") %>% arrange(-logr0.mean) %>% slice(1:2)
}

#5) Isolate highest r0 value lines for no gradient treatment at pH 4
{
  #Filter data to get evolution lines with gradient
  nogradient <- filter(plotdata, Gradient != "Gradient" & Treatment == "EVO")
  nogradient <- filter(nogradient, testpH==6.5)
  nogr.sex.geneflow <- filter(nogradient, Sex=="Sex" & Gene.Flow == "Gene flow") %>% arrange(-logr0.mean) %>% slice(1:2)
  nogr.nosex.geneflow <- filter(nogradient, Sex!="Sex" & Gene.Flow == "Gene flow") %>% arrange(-logr0.mean) %>% slice(1:2)
  nogr.sex.nogeneflow <- filter(nogradient, Sex=="Sex" & Gene.Flow != "Gene flow") %>% arrange(-logr0.mean) %>% slice(1:2)
  nogr.nosex.nogeneflow <- filter(nogradient, Sex!="Sex" & Gene.Flow != "Gene flow") %>% arrange(-logr0.mean) %>% slice(1:2)
}