setwd("/Users/fazal2/Desktop/MethylationData_Updated/Methylation_CorelationAnalysis_Selected/TCGA/")
library(tidyverse)
library(plyr)
library(readr)
library(stringr)
library(reshape2)
suppressPackageStartupMessages(library(dplyr))
library(data.table)
library(vroom)

#read multiplefiles
files=fs::dir_ls(glob = "*.txt")
files
data <- vroom(files,col_select = list(c(1:3),probe=`Composite Element REF`), id="SampleID")

#group by each sample
bysample=data %>% 
  group_by(SampleID) 
head(bysample)

#split each group
fit=group_split(bysample, .keep = TRUE)
length(fit)

#read DMPs to compare
DMPs<-read_tsv("/Users/fazal2/Desktop/tidy_practice/B3_Case1_Beta0.5_0.05FDR.txt")
head(DMPs)
test=list()
#compare each sample with DMPs and write in output
for(i in 1:length(fit)){
  test[[i]]=fit[[i]] %>% 
  filter(probe %in% DMPs$Probe) %>% 
  group_by(SampleID) %>% 
  pivot_wider(names_from = SampleID, values_from = Beta_value)
}

full_join(test,by="probe")
full=Reduce(merge,test)
write_tsv(full, path = "test.tsv")