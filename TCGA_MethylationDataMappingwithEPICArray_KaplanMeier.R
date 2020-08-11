setwd("/Users/fazal2/Desktop/MethylationData_Updated/Methylation_CorelationAnalysis_Selected/TCGA/")
library(tidyverse)
require(data.table)
library("survival")
library("survminer")

DMPs<-read_tsv("/Users/fazal2/Desktop/tidy_practice/B3_Case1_Beta0.4_0.05FDR.txt")
clinical<-read_tsv("/Users/fazal2/Desktop/tidy_practice/tgct_tcga_clinical_data.tsv")

colnames(clinical)<-c("ID","PID","Sample_ID","CancerType","DSF_Time","DFS_Status")
#reform_clin=melt(clinical)
#colnames(reform_clin)<-c("ID","PID","Sample_ID","CancerType","DFS_Status","variable","Time")
head(DMPs)
#head(reform_clin)

#create a list of the files from your target directory
file_list <- list.files()
test=list()
for (i in 1:length(file_list)){
  temp_data <- fread(file_list[i], stringsAsFactors = F, select = c(1,2))
  colnames(temp_data)<-c("Probe","beta")
  temp_data$SampleID<-file_list[i]
  test[[i]]=temp_data %>% 
  filter(Probe %in% DMPs$Probe) %>% 
  group_by(SampleID) %>% 
  pivot_wider(names_from = SampleID, values_from = beta) %>% 
  drop_na()
}
full=Reduce(merge,test)
dim(full)

#calculate mean of each sample
Samples_Mean=full %>% 
  summarise(across(where(is.numeric), mean))

#extract sample Names from file names
names(Samples_Mean)<-substr(names(Samples_Mean), start=1, stop=15)
head(Samples_Mean)
reform_means=melt(Samples_Mean)
colnames(reform_means)<-c("Sample_ID","Beta")

#mapp clinical data
Methy_Clinical=clinical %>% 
  filter(Sample_ID%in%reform_means$Sample_ID) %>% 
  full_join(reform_means,by="Sample_ID") %>% 
  drop_na() %>% 
  select(c("Sample_ID","CancerType","DFS_Status","DSF_Time","Beta"))

#spliting data by median Beta value
Median_Beta=median(Methy_Clinical$Beta)
Methy_Clinical$level = NA # creates a new variable filled with NAs
high = Methy_Clinical$Beta>=Median_Beta
low =  Methy_Clinical$Beta<Median_Beta
Methy_Clinical$level[high]="High"
Methy_Clinical$level[low]="Low"

Methy_Clinical$event = NA # creates a new variable filled with NAs
event_occured = Methy_Clinical$DFS_Status=="1:Recurred/Progressed"
no_event =  Methy_Clinical$DFS_Status=="0:DiseaseFree"
Methy_Clinical$event[event_occured]="1"
Methy_Clinical$event[no_event]="0"


#Methy_Clinical$CancerType <- as.character(Methy_Clinical$CancerType)
Methy_Clinical$CancerType[Methy_Clinical$CancerType == "Non-Seminomatous Germ Cell Tumor" | Methy_Clinical$CancerType == "Non-Seminomatous Germ Cell Tumo" | Methy_Clinical$CancerType == "Embryonal Carcinoma" ] <- "Non-Seminoma"

write_tsv(Methy_Clinical, path = "/Users/fazal2/Desktop/tidy_practice/B3_Case1_Beta0.4_0.05FDR_Methylation_Clinical.tsv")

#run KM
fit <- survfit(Surv(as.numeric(Methy_Clinical$DSF_Time), as.numeric(Methy_Clinical$event)) ~ Methy_Clinical$level, data = Methy_Clinical)
print(fit)
summary(fit)

ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


surv_diff <- survdiff(Surv(as.numeric(Methy_Clinical$DSF_Time), as.numeric(Methy_Clinical$event)) ~ Methy_Clinical$level, data = Methy_Clinical)
surv_diff

