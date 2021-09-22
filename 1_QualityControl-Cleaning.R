#1) Clear memory and set number of cores
{
  rm(list = ls())
  library(parallel)
  cores <- detectCores() - 1
  library(tidyverse)
}

#2) Set working directory
{
  setwd("/media/felix/Elements/12_RangeExpansionGenomics")
}

#3) Get list with all fastq files
{
  block1 <- list.files(path = "2_data/Project_3191_20190919_230834_b", recursive = T, pattern = "\\.fastq$", full.names = T)
  block2 <- list.files(path = "2_data/Project_3191_20191030_113041_b", recursive = T, pattern = "\\.fastq$", full.names = T)
}

all.fastq <- c(block1, block2)
all.fastq <- paste(getwd(), all.fastq, sep = "/")
#4) Perform fastqc quality controll for all files (parallel, over all -1 available cores)
{
  block1 <- paste(getwd(), block1, sep = "/")
  qc.fastq <- function(i){system(paste("fastqc", i))}
  # mclapply(qc.fastq, qc.fastq, mc.cores = cores)
}

#5) Perform multiqc for each folder
{
  #Needs to be run maually from shell, doesn't work from r --> try to find out why
}

#6) Exclude the two files without reads
{
  bad1 <- paste(getwd(), "/2_data/Project_3191_20190919_230834_b/resources/workunit_202888/resource_1325223/20190830.B-MIX_EVO_35_P3_R1.fastq", sep = "")
  bad2 <- paste(getwd(), "/2_data/Project_3191_20190919_230834_b/resources/workunit_202888/resource_1325224/20190830.B-MIX_EVO_35_P3_R2.fastq", sep = "")
  
  block1 <- paste(getwd(), block1, sep="/")
  all.fastq.filt <- all.fastq[-which(all.fastq %in% c(bad1, bad2))]
  block1.filt <- block1[-which(block1 %in% c(bad1, bad2))]
}

#7) Filter Fastqc files using trimmomatic
{
  #Filtering is just done by filtering for illumina adapters using illuminaclip + additional quality check (leading, trailing, minlen)
  trimming <- function(i, fastq.list){
    raw1 <- fastq.list[i]
    raw2 <- fastq.list[i+1]
    
    #Create a subfolder for trimmed files
    folder1 <- unlist(strsplit(fastq.list[i], "/"))
    folder1 <- paste(c(folder1[1:length(folder1)-1], "Trimmed"), collapse = "/")
    system(paste("mkdir -p", folder1))
    
    folder2 <- unlist(strsplit(fastq.list[i+1], "/"))
    folder2 <- paste(c(folder2[1:length(folder2)-1], "Trimmed"), collapse = "/")
    system(paste("mkdir -p", folder2))
    
    #Create the names of the output files
    o.R1.paired <- paste(folder1, "R1_paired.fastq", sep="/")
    o.R1.unpaired <- paste(folder1, "R1_unpaired.fastq", sep="/")
    o.R2.paired <- paste(folder2, "R2_paired.fastq", sep="/")
    o.R2.unpaired <- paste(folder2, "R2_unpaired.fastq", sep="/")
    
    #Run trimmomatic command for the paired read files
    command <- paste("java -jar /home/felix/anaconda3/share/trimmomatic/trimmomatic.jar PE -phred33", raw1, raw2, o.R1.paired, o.R1.unpaired, o.R2.paired, o.R2.unpaired, "ILLUMINACLIP:/home/felix/anaconda3/share/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36")
    system(command)
    
    #Run fastqc for the paired read files
    system(paste("fastqc", o.R1.paired))
    system(paste("fastqc", o.R2.paired))
  }
  block2.b <- block2[!grepl("Trimmed", block2)]
  block2.b <- block2.b[!grepl("trimmed", block2.b)]
  block2.b <- paste(getwd(), block2.b, sep="/")
  # trimming(1, all.fastqc.filt)
  # block1.b <- block1[c(81, 84)]
  #mclapply(seq(1, length(block2.b), by=2), trimming, block2.b, mc.cores = 1)
  # trimming(1, block1.b)
}

#8) Map reads to reference genome using BWA
{
  #Define location of reference genome to map to
  ref.genome.fna <- paste(getwd(), "2_data/TetThermRefGenome/GCF_000189635.1_JCVI-TTA1-2.2_genomic.fna", sep="/")
  ref.genome.gbff <- paste(getwd(), "2_data/TetThermRefGenome/GCF_000189635.1_JCVI-TTA1-2.2_genomic.gbff", sep="/")
  ref.genome.fa <- paste(getwd(), "2_data/TetThermRefGenome/GCF_000189635.1_JCVI-TTA1-2.2_genomic.fa", sep="/")
  
  #Create Samtools, PICARD and BWA index
  {
    #BWA
    #system(paste("bwa index -a is", ref.genome.fna))
    
    #Samtools
    #system(paste("samtools faidx", ref.genome.fna))
    
    #PICARD
    picard.jar.location <- paste(getwd(), "5_Packages/picard.jar", sep = "/")
    #system(paste("java -jar ", picard.jar.location, " CreateSequenceDictionary R=", ref.genome.fna, " O=", getwd(), "/2_data/TetTherm_RefGenome/", "GCF_000189635.1_JCVI-TTA1-2.2_genomic.dict", sep=""))
  }
  
  #Perform read mapping using BWA
  
  mapping <- function(i, fastqc.list, ref.genome.fna){
    foldermain1 <- unlist(str_split(fastqc.list[i], "/"))
    foldermain2 <- unlist(str_split(fastqc.list[i+1], "/"))
    folder1 <- paste(paste(foldermain1[1:length(foldermain1)-1], collapse = "/"), "/Bamfiles", sep = "")
    paired1 <- paste(paste(foldermain1[1:length(foldermain1)-1], collapse = "/"), "/Trimmed/R1_paired.fastq", sep = "")
    paired2 <- paste(paste(foldermain2[1:length(foldermain2)-1], collapse = "/"), "/Trimmed/R2_paired.fastq", sep = "")
    system(paste("mkdir -p", folder1))
    output <- paste(folder1, "/MappedReads.sam", sep="")
    # system(paste("bwa mem ", ref.genome.fna,  paired1, paired2,  ">",  output, sep=" "))
    
    #Convert bam to sam
    system(paste("samtools view -Sb", output, ">", paste(folder1, "/Mapped_Reads.bam", sep=""))) #depending on version of samtools
    
    #Sort bam file
    #system(paste("samtools sort", paste(folder1, "/MappedReads.bam", sep=""), "-o", paste(folder1, "/SortedMappedReads.bam", sep="")))
    
    #Index bam file
    #system(paste("samtools index", paste(folder1, "/SortedMappedReads.bam", sep="")))
    
    #Remove sam file to save space
    #system(paste("rm", output))
  }
  
  #Do parallelised mapping
  #mapping(1, all.fastq.filt, ref.genome.fna)
  all.fastq.filt <- all.fastq.filt[!grepl("trimmed", all.fastq.filt)]
  block1.b <- block1.filt[!grepl("Trimmed", block1.filt)]
  # block2.b <- block2[!grepl("Trimmed", block2)]
  # block2.b <- block2.b[!grepl("trimmed", block2.b)]
  # mclapply(seq(25, 38, by=2), mapping, block1.b, ref.genome.fna, mc.cores = 7)
  # mapping(1, block1.b, ref.genome.fna)
}

#9) Merge bam files of replicate populations
{
  #Function to merge bam files
  
  # mergeBam <- function(foldername, bam1, bam2, ident){
  #   if (foldername == "MIX_EVO_35_P3_AACGAGGCCG-CACGGAACAA"){
  #     inputbam <- grep(ident[foldername, 3], bam2, value=TRUE)[1]
  #     output <- paste(getwd(), "/2_data/MergedBam/", foldername, "_merged.bam", sep="")
  #     #system(paste("mv ", inputbam, output))
  #   }
  #   else {
  #     input1 <- grep(ident[foldername, 2], bam1, value=TRUE)[1]
  #     input2 <- grep(ident[foldername, 3], bam2, value=TRUE)[1]
  #     output <- paste(getwd(), "/2_data/MergedBam/", foldername, "_merged.bam", sep="")
  # 
  #     #system(paste("samtools merge", output, input1, input2))
  #   }
  #   output2 <- paste(getwd(), "/2_data/MergedBam/", foldername, "_ordered.bam", sep="")
  #   output3 <- paste(getwd(), "/2_data/MergedBam/", foldername, "_ordered_Mdup.bam", sep="")
  #   duplicates <- paste(getwd(), "/2_data/MergedBam/", foldername, "_duplicates.bam", sep="")
  #   system(paste("samtools sort ", output, " -o ", output2, sep=""))
  #   system(paste("samtools index", output2))
  #   system(paste("java -jar ", picard.jar.location, " MarkDuplicates ", " I=", output2, " O=", output3 , " M=", duplicates, sep=""))
  # 
  # }
  # 
  # #Make a list with all filenames of filtered bam files
  # # bam1 <- paste(getwd(), list.files(path = "2_data/Project_3191_20190919_230834_b", recursive = T, pattern = "\\.bam$", full.names = T), sep="/")
  # # bam2 <- paste(getwd(), list.files(path = "2_data/Project_3191_20191030_113041_b", recursive = T, pattern = "\\.bam$", full.names = T), sep="/")
  # 
  # all.fastq.filt <- all.fastq.filt[!grepl("Trimmed", all.fastq.filt)]
  # all.fastq.filt <- all.fastq.filt[!grepl("trimmed", all.fastq.filt)]
  # 
  # listnames <- paste(getwd(), list.files(path = "2_data", recursive = T, pattern = "\\R1.fastq$", full.names = T), sep="/")
  # ident <- data.frame()
  # ident <- mutate(ident, name = NA, location1 = NA, location2 = NA)
  # for (i in strsplit(listnames, split = "/")){
  #   temp1 <- as.character(i[13])
  #   temp1 <- substr(temp1, 12, nchar(temp1)-9)
  # 
  #   temp2 <- as.character(i[12])
  #   if (temp1 != "MIX_EVO_35_P3"){
  #     if(temp1 %in% row.names(ident)){
  #       ident[temp1, ] <- c(temp1, ident[temp1, 2], temp2)
  #     }
  #     else {
  #       ident[temp1, ] <- c(temp1, temp2, NA)
  #     }
  #   }
  # }
  # ident_names <- row.names(ident)
  # ident["MIX_EVO_35_P3_AACGAGGCCG-CACGGAACAA", ] <- c(ident["MIX_EVO_35_P3", 1], ident["MIX_EVO_35_P3", 3], ident["MIX_EVO_35_P3", 2])

  #Make a folder to save merged bam folders
  # system(paste("mkdir -p", paste(getwd(), "/2_data/MergedBam", sep="")))

  
  #Make new function to do just ordering and marking duplicates for bam files
  
  bamfiles <- paste(getwd(), list.files(path = "2_data/MergedBam/", recursive = T, pattern = "\\.bam$", full.names = T), sep="/")
  {
    MarkDupl <- function(bamfile){
      filename <- substr(bamfile, 1, nchar(bamfile)-4)
      output2 <- paste(bamfile, "_ordered.bam", sep="")
      output3 <- paste(bamfile, "_ordered_Mdup.bam", sep="")
      duplicates <- paste(bamfile, "_duplicates.bam", sep="")
      # system(paste("samtools sort ", bamfile, " -o ", output2, sep=""))
      # system(paste("samtools index", output2))
      system(paste("java -jar ", picard.jar.location, " MarkDuplicates ", " I=", output2, " O=", output3 , " M=", duplicates, sep=""))
    }
  }
  #perform merging step for all bam files
  # mclapply(bamfiles, MarkDupl, mc.cores = 1)
} 
Listfiles <-  paste(getwd(), list.files(path = "2_data/MergedBam/", recursive = T, pattern = "Mdup\\.bam$", full.names = T), sep="/")

#10) Do variant calling without -no--BAQ flag
{
  #10.1) Create output for all sites
  {
    doSamtools <- function(i, outputfolder, ref.genome.fa){
      filename <- unlist(strsplit(i, split = "/"))[length(unlist(strsplit(i, split = "/")))]
      outputfile <- paste(substr(filename, start = 1, stop = nchar(filename)-4), ".vcf", sep="")
      output <- paste(outputfolder, outputfile, sep="/")
      system(paste("bcftools mpileup -a AD,DP,SP,FORMAT/AD,FORMAT/DV,FORMAT/DPR -Ou -f", ref.genome.fa, i, "| bcftools call -m -f GQ,GP -A -o", output))
    }
    
    #Do Samtools variant calling for ancestors 
    outputfolder <- paste(getwd(), "2_data/VariantCallingComplete", sep="/")
    Listfiles <-  paste(getwd(), list.files(path = "2_data/MergedBam/", recursive = T, pattern = "Mdup\\.bam$", full.names = T), sep="/")
    mclapply(Listfiles, doSamtools, outputfolder, ref.genome.fa, mc.cores=6)
  }
  
  #10.2) Do the same, but only output variant sites
  {
    doSamtools <- function(i, outputfolder, ref.genome.fa){
      filename <- unlist(strsplit(i, split = "/"))[length(unlist(strsplit(i, split = "/")))]
      outputfile <- paste(substr(filename, start = 1, stop = nchar(filename)-4), ".vcf", sep="")
      output <- paste(outputfolder, outputfile, sep="/")
      system(paste("bcftools mpileup -a AD,DP,SP,FORMAT/AD,FORMAT/DV,FORMAT/DPR -Ou -f", ref.genome.fa, i, "| bcftools call -m -f GQ,GP -A -v -o", output))
    }
    
    #Do Samtools variant calling for ancestors 
    outputfolder <- paste(getwd(), "2_data/VariantCalling_OnlyDeviating", sep="/")
    Listfiles <-  paste(getwd(), list.files(path = "2_data/MergedBam/", recursive = T, pattern = "Mdup\\.bam$", full.names = T), sep="/")
    # mclapply(Listfiles, doSamtools, outputfolder, ref.genome.fa, mc.cores=6)
  }
  
}

#11) Do variant calling with -no--BAQ flag
# {
#   #11.1) Create output for all sites
#   {
#     doSamtools <- function(i, outputfolder, ref.genome.fa){
#       filename <- unlist(strsplit(i, split = "/"))[length(unlist(strsplit(i, split = "/")))]
#       outputfile <- paste(substr(filename, start = 1, stop = nchar(filename)-4), ".vcf", sep="")
#       output <- paste(outputfolder, outputfile, sep="")
#       system(paste("bcftools mpileup -a AD,DP,SP,FORMAT/AD,FORMAT/DV,FORMAT/DPR -Ou -f", ref.genome.fa, i, "| bcftools call -m -f GQ,GP -A -o", output))
#     }
#     
#     #Do Samtools variant calling for ancestors 
#     outputfolder <- paste(getwd(), "2_data/7_Allel_frequencies_Ancestors/noBAQ_VariantCalling/Samtools_Complete_noBAQ", sep="/")
#     Listfiles <-  paste(getwd(), list.files(path = "2_data/MergedBam/", recursive = T, pattern = "Mdup\\.bam$", full.names = T), sep="/")
#     #mclapply(Listfiles, doSamtools, outputfolder, ref.genome.fa, mc.cores=6)
#   }
#   
#   #11.2) Do the same, but only output variant sites
#   {
#     doSamtools <- function(i, outputfolder, ref.genome.fa){
#       filename <- unlist(strsplit(i, split = "/"))[length(unlist(strsplit(i, split = "/")))]
#       outputfile <- paste(substr(filename, start = 1, stop = nchar(filename)-4), ".vcf", sep="")
#       output <- paste(outputfolder, outputfile, sep="")
#       system(paste("bcftools mpileup -a AD,DP,SP,FORMAT/AD,FORMAT/DV,FORMAT/DPR -Ou -f", ref.genome.fa, i, "| bcftools call -m -f GQ,GP -A -v -o", output))
#     }
#     
#     #Do Samtools variant calling for ancestors 
#     outputfolder <- paste(getwd(), "2_data/7_Allel_frequencies_Ancestors/noBAQ_VariantCalling/Samtools_reduced_noBAQ", sep="/")
#     Listfiles <-  paste(getwd(), list.files(path = "2_data/MergedBam/", recursive = T, pattern = "Mdup\\.bam$", full.names = T), sep="/")
#     #mclapply(Listfiles, doSamtools, outputfolder, ref.genome.fa, mc.cores=6)
#   }
#   
# }