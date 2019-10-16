#Main code
#####################################################################################################
#####################################################################################################
#Figure 1A was drawn manually based on the source data.
#####################################################################################################
#####################################################################################################
#Figure 1B: Radar plot for TCGA cancer data
rm(list=ls())
rm(list=ls())
rm(list=ls())
library(ggplot2)
path_raw<-"E:\\PKM_data_code\\Source_data\\Figure_1B\\"
setwd(path_raw)
matrix<-as.matrix(read.csv("Figure_1B_TCGA.txt",header=T,row.names = 1,sep="\t"))
colnames(matrix)[1]="PKM1"
colnames(matrix)[2]="PKM2"
matrix<-log2(matrix+1)
require(tidyverse)
a<-matrix %>%
  as_tibble(rownames = "Tissue") %>%
  gather(Gene, value,  -Tissue) %>%
  mutate(order = case_when(Gene == "PKM1" ~ 1,
                           Gene == "PKM2" ~ 2,
                           Gene == "ENST00000389093" ~ 3,
                           Gene == "ENST00000561609" ~ 4,
                           Gene == "ENST00000568883" ~ 5,
                           Gene == "ENST00000568459" ~ 6,
                           Gene == "ENST00000562997" ~ 7,
                           Gene == "ENST00000449901" ~ 8,
                           Gene == "ENST00000564178" ~ 9,
                           Gene == "ENST00000565154" ~ 10,
                           Gene == "ENST00000565184" ~ 11,
                           Gene == "ENST00000566809" ~ 12,
                           Gene == "ENST00000567087" ~ 13,
                           Gene == "ENST00000569050" ~ 14
  )) %>%
  mutate(order = as.factor(order)) 

ggplot(data=a,aes(x = Tissue, group = Gene, y = value, colour = Gene, fill = Gene)) +
  geom_polygon(data = a %>%
                 filter(order == "14"), aes(color=Gene),size=1.7, alpha = 0) +
  geom_polygon(data = a %>%
                 filter(order == "13"), aes(color=Gene),size=1.7, alpha = 0) +
  geom_polygon(data = a %>%
                 filter(order == "12"), aes(color=Gene),size=1.7, alpha = 0) +
  geom_polygon(data = a %>%
                 filter(order == "11"), aes(color=Gene),size=1.7, alpha = 0) +
  geom_polygon(data = a %>%
                 filter(order == "10"), aes(color=Gene),size=1.7, alpha = 0) +
  geom_polygon(data = a %>%
                 filter(order == "9"), aes(color=Gene),size=1.7, alpha = 0) +
  geom_polygon(data = a %>%
                 filter(order == "8"), aes(color=Gene),size=1.7, alpha = 0) +
  geom_polygon(data = a %>%
                 filter(order == "7"), aes(color=Gene),size=1.7, alpha = 0) +
  geom_polygon(data = a %>%
                 filter(order == "6"), aes(color=Gene),size=1.7, alpha = 0) +
  geom_polygon(data = a %>%
                 filter(order == "5"), aes(color=Gene),size=1.7, alpha = 0) +
  geom_polygon(data = a %>%
                 filter(order == "4"), aes(color=Gene),size=1.7, alpha = 0) +
  geom_polygon(data = a %>%
                 filter(order == "3"), aes(color=Gene),size=1.7, alpha = 0) +
  geom_polygon(data = a %>%
                 filter(order == "2"), aes(color=Gene),size=1.7, alpha = 0) +
  geom_polygon(data = a %>%
                 filter(order == "1"), aes(color=Gene),size=1.7, alpha = 0) +
  labs(x = "", y = "") +
  coord_polar() +
  scale_colour_manual(values = c("PKM1" = "yellow green", "PKM2" = "red","ENST00000389093"= "blue","ENST00000561609"="magenta","ENST00000568883"="gold","ENST00000568459"="purple","ENST00000562997"="cyan", "ENST00000449901"="gray","ENST00000564178"="gray","ENST00000565154"="gray","ENST00000565184"="gray","ENST00000566809"="gray","ENST00000567087"="gray","ENST00000569050"="gray")) +
  #scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10,11)) +
  guides(col = guide_legend(ncol = 7)) +
  theme(
    legend.position = "top",
    axis.text = element_text(face = "bold",size=20),
    panel.background = element_blank(),
    panel.grid.major.y = element_line(color = alpha("gray", alpha = 0.9)),
    panel.grid.major.x = element_line(color = alpha("gray", alpha = 0.9)),
    axis.text.y = element_text(),
    axis.ticks.y = element_line()
  )
path_out<-paste0(path_raw,"Output\\")
dir.create(path_out)
setwd(path_out)
ggsave("./Output_Figure_1B_TCGA.pdf", width = 12, height = 12, dpi=800)
#####################################################################################################
#Figure_1B: Radar plot for GTEX normal tissues
rm(list=ls())
rm(list=ls())
rm(list=ls())
library(ggplot2)
path_raw<-"E:\\PKM_data_code\\Source_data\\Figure_1B\\"
setwd(path_raw)
matrix<-as.matrix(read.csv("Figure_1B_GTEx.txt",header=T,row.names = 1,sep="\t"))
colnames(matrix)[1]="PKM1"
colnames(matrix)[2]="PKM2"
matrix<-log2(matrix+1)
require(tidyverse)
a<-matrix %>%
  as_tibble(rownames = "Tissue") %>%
  gather(Gene, value,  -Tissue) %>%
  mutate(order = case_when(Gene == "PKM1" ~ 1,
                           Gene == "PKM2" ~ 2,
                           Gene == "ENST00000389093" ~ 3,
                           Gene == "ENST00000561609" ~ 4,
                           Gene == "ENST00000568883" ~ 5,
                           Gene == "ENST00000568459" ~ 6,
                           Gene == "ENST00000562997" ~ 7,
                           Gene == "ENST00000449901" ~ 8,
                           Gene == "ENST00000564178" ~ 9,
                           Gene == "ENST00000565154" ~ 10,
                           Gene == "ENST00000565184" ~ 11,
                           Gene == "ENST00000566809" ~ 12,
                           Gene == "ENST00000567087" ~ 13,
                           Gene == "ENST00000569050" ~ 14
  )) %>%
  mutate(order = as.factor(order)) 
ggplot(data=a,aes(x = Tissue, group = Gene, y = value, colour = Gene, fill = Gene)) +
  #geom_point()+
  geom_polygon(data = a %>%
                 filter(order == "14"), aes(color=Gene),size=1.7, alpha = 0) +
  geom_polygon(data = a %>%
                 filter(order == "13"), aes(color=Gene),size=1.7, alpha = 0) +
  geom_polygon(data = a %>%
                 filter(order == "12"), aes(color=Gene),size=1.7, alpha = 0) +
  geom_polygon(data = a %>%
                 filter(order == "11"), aes(color=Gene),size=1.7, alpha = 0) +
  geom_polygon(data = a %>%
                 filter(order == "10"), aes(color=Gene),size=1.7, alpha = 0) +
  geom_polygon(data = a %>%
                 filter(order == "9"), aes(color=Gene),size=1.7, alpha = 0) +
  geom_polygon(data = a %>%
                 filter(order == "8"), aes(color=Gene),size=1.7, alpha = 0) +
  geom_polygon(data = a %>%
                 filter(order == "7"), aes(color=Gene),size=1.7, alpha = 0) +
  geom_polygon(data = a %>%
                 filter(order == "6"), aes(color=Gene),size=1.7, alpha = 0) +
  geom_polygon(data = a %>%
                 filter(order == "5"), aes(color=Gene),size=1.7, alpha = 0) +
  geom_polygon(data = a %>%
                 filter(order == "4"), aes(color=Gene),size=1.7, alpha = 0) +
  geom_polygon(data = a %>%
                 filter(order == "3"), aes(color=Gene),size=1.7, alpha = 0) +
  geom_polygon(data = a %>%
                 filter(order == "2"), aes(color=Gene),size=1.7, alpha = 0) +
  geom_polygon(data = a %>%
                 filter(order == "1"), aes(color=Gene),size=1.7, alpha = 0) +
  labs(x = "", y = "") +
  
  coord_polar() +
  scale_colour_manual(values = c("PKM1" = "yellow green", "PKM2" = "red","ENST00000389093"= "blue","ENST00000561609"="magenta","ENST00000568883"="gold","ENST00000568459"="purple","ENST00000562997"="cyan", "ENST00000449901"="gray","ENST00000564178"="gray","ENST00000565154"="gray","ENST00000565184"="gray","ENST00000566809"="gray","ENST00000567087"="gray","ENST00000569050"="gray")) +
  #scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10)) +
  guides(col = guide_legend(ncol = 7)) +
  theme(
    legend.position = "top",
    axis.text = element_text(face = "bold",size=20),
    panel.background = element_blank(),
    panel.grid.major.y = element_line(color = alpha("gray", alpha = 0.9)),
    panel.grid.major.x = element_line(color = alpha("gray", alpha = 0.9)),
    axis.text.y = element_text(),
    axis.ticks.y = element_line()
  )
path_out<-paste0(path_raw,"Output\\")
setwd(path_out)
ggsave("./Output_Figure_1B_GTEx.pdf", width = 12, height = 12,dpi=800)
##############################################################################################
##############################################################################################
#figure 2A: clustering plot based on KM P values
rm(list=ls())
rm(list=ls())
rm(list=ls())

path_raw<-"E:\\PKM_data_code\\Source_data\\Figure_2A\\"
setwd(path_raw)
mean_int_trans<-as.matrix(read.csv("Figure_2A_cluster_mean_trans_exp.txt",header=T,row.names = 1,sep="\t"))
coef_end<-as.matrix(read.csv("Figure_2A_cluster_KM_coef.txt",header=T,row.names = 1,sep="\t"))
p_end<-as.matrix(read.csv("Figure_2A_cluster_KM_p.txt",header=T,row.names = 1,sep="\t"))

index_1<-which(mean_int_trans<5)
p_end[index_1]=1
log_p_value=-log10(p_end)

index<-which(coef_end<0)
log_p_value[index]=-log_p_value[index]
#use the P value of KM do clustering analysis
library(corrplot)
library(Hmisc)
library(gplots)
path_out<-paste0(path_raw,"Output\\")
dir.create(path_out)
setwd(path_out)

tiff(filename = "Figure_2A_cluster.tiff",width=1500,height=1200)
hr <- hclust(as.dist(1-cor(t(log_p_value), method="spearman")), method="ward.D2")
my_palete <- colorRampPalette(c("blue","white","red3"))(n=299)
col_breaks = c(seq(-1.3,-1,length=100),
               + seq(-0.999,1,length=100),
               + seq(1.009,1.3,length=100))  
heatmap.2(log_p_value,Rowv=FALSE, Colv=FALSE,scale="none",density.info="none",trace="none",col=my_palete,breaks=col_breaks,margins =c(15,15),cexRow=1.5,cexCol=1.5,sepwidth=c(0.5,0.5),sepcolor="black")
dev.off()
##################################################################################
#Figure 2A: Kaplan-Meier plot for TCGA KIRC
#survival data (days)
rm(list=ls())
rm(list=ls())
rm(list=ls())

library("grid")
library("xlsx")
library("XLConnect")
library("ggplot2")
library("ggvis")
library("rgl")
library("dplyr")
library("tidyr")
library("stringr")
library("lubridate")
require(gplots)
library("survival")
require("xlsx")
library("DESeq2")
library("biomaRt")
library("DESeq2")
library("piano")
library("Biobase")
setwd("E:\\PKM_data_code\\Code\\")
source('Cheng_toolbox_beta.R')
gene_list<-as.matrix(c("ENST00000335181","ENST00000561609","ENST00000389093","ENST00000568883"))
dataDir<-"E:\\PKM_data_code\\Source_data\\Figure_2A\\"
setwd(dataDir)
TXT_name<-"Figure_2A_KM_plot_TCGA_KIRC.txt"
path_out<-paste0(dataDir,"Output\\")
setwd(path_out)

  for (j in 1:length(gene_list)){
    output<-Cheng_generateSurvInputfromTCGA(gene_list[j,1],TXT_name,dataDir)
    result<-Cheng_generateKMplot(output,outFile=paste0("TCGA_KIRC","_",gene_list[j,1]))
  }
##################################################################################
##################################################################################
#Figure 2B: bubble plot
rm(list=ls())
rm(list=ls())
rm(list=ls())

trans_list<-as.matrix(c("ENST00000335181","ENST00000561609","ENST00000389093","ENST00000568883"))
path_raw<-"E:\\PKM_data_code\\Source_data\\Figure_2B\\"
setwd(path_raw)
path_out<-paste0(path_raw,"Output\\")
dir.create(path_out)
go_count_end<-as.matrix(read.csv("Figure_2B_GoTerm_count_inform.txt",header=F,sep="\t"))

for(i in 1:length(trans_list)){
setwd(path_raw)
int_trans<-trans_list[i]
file_name<-paste0("Figure_2B_GoTerm_inform_",int_trans,".txt")
matrix<-as.matrix(read.csv(file_name,header=T,sep="\t"))
goid_all<-as.matrix(matrix[,"goid_all"])
Direction<-as.matrix(matrix[,"Direction"])
Generality<-as.matrix(matrix[,"Generality"])
rownames(Direction)<-matrix[,1]
rownames(Generality)<-matrix[,1]
mode(Direction)<-"numeric"
mode(Generality)<-"numeric"

index<-match(goid_all,go_count_end[,1])
goid_count<-as.matrix(go_count_end[index,2])
mode(goid_count)<-"numeric"
size=log10(goid_count)
bubbleRes<-cbind(Direction,Generality,size)
colnames(bubbleRes)<-c("Direction","Generality","size")
bubbleRes<-as.data.frame(bubbleRes)

tes=gsub(' ','.',rownames(read.csv("Figure_2B_common_GO_terms.txt",row.names = 1,sep='\t')))
tes=gsub(',','.',tes)
tes=gsub('-','.',tes)
bubbleRes2=na.omit(bubbleRes[tes,])
bubbleRes=bubbleRes[!rownames(bubbleRes)%in%tes,]
if(int_trans=="ENST00000335181"|int_trans=="ENST00000561609"){
  color_com<-"dark green"
}else{
  color_com<-"red"
}
require(ggplot2)
ggOut = ggplot(bubbleRes, aes(x=Direction, y=Generality, size=size))+
  geom_jitter(aes(size =size,colour=Generality,alpha=.02))+
  scale_colour_gradient(guide = FALSE)+
  xlim(-25,25)+ylim(0,25)+theme(axis.ticks.length = unit(.2, "cm"), axis.line = element_line(colour = "black"), text = element_text(size=17), legend.position="none", panel.background=element_rect(fill="white"))+
  geom_vline(xintercept = 0,linetype="dashed", color = "black", size=1)+
  geom_jitter(data = bubbleRes2, colour=color_com, aes(size =size,colour=Generality,alpha=.02))
#print(ggOut)
setwd(path_out)
ggsave(paste0("./Output_Figure_2B_BubblePlot_",int_trans,".pdf"), width = 9, height = 9)

}

##################################################################################
##################################################################################
#Figure 3A: Kaplan-Meier plot in Japanese cohort
rm(list=ls())
rm(list=ls())
rm(list=ls())

library("grid")
library("xlsx")
library("XLConnect")
library("ggplot2")
library("ggvis")
library("rgl")
library("dplyr")
library("tidyr")
library("stringr")
library("lubridate")
require(gplots)
library("survival")
require("xlsx")
library("DESeq2")
library("biomaRt")
library("DESeq2")
library("piano")
library("Biobase")
setwd("E:\\PKM_data_code\\Code\\")
source('Cheng_toolbox_beta.R')
gene_list<-as.matrix(c("ENST00000335181","ENST00000561609","ENST00000389093","ENST00000568883"))
dataDir<-"E:\\PKM_data_code\\Source_data\\Figure_3A\\"
setwd(dataDir)
TXT_name<-"Figure_3A_KM_plot_Japanese_KIRC.txt"
path_out<-paste0(dataDir,"Output\\")
dir.create(path_out)
setwd(path_out)

for (j in 1:length(gene_list)){
  output<-Cheng_generateSurvInputfromTCGA(gene_list[j,1],TXT_name,dataDir)
  result<-Cheng_generateKMplot(output,outFile=paste0("Japanese_KIRC","_",gene_list[j,1]))
}

##################################################################################
##################################################################################
#Figure 3B: KM plot for validating biomarker
rm(list=ls())
rm(list=ls())
rm(list=ls())

library(survival)
th<-2# at least 2 of 4 prognostic transcripts vote high risk
int_trans<-as.matrix(c("ENST00000335181","ENST00000561609","ENST00000389093","ENST00000568883"))
path_raw<-"E:\\PKM_data_code\\Source_data\\Figure_3B\\"
setwd(path_raw)
cutoff<-as.matrix(read.csv("Fgiure_3B_cutoff_TCGA_KIRC.txt",header=F,row.names = 1,sep="\t"))
#cutoff<-as.matrix(read.csv("Fgiure_3B_cutoff_Japanese.txt",header=F,row.names = 1,sep="\t"))
exp<-as.matrix(read.csv("Fgiure_3B_TCGA_KIRC_exp_survival.txt",header=T,row.names = 1,sep="\t"))
#exp<-as.matrix(read.csv("Fgiure_3B_Japanese_KIRC_exp_survival.txt",header=T,row.names = 1,sep="\t"))
survival<-exp[,c("OS","status")]
exp<-exp[,int_trans]

loc_1=exp[,int_trans[1]]<cutoff[1] #high risk is 1
loc_2=exp[,int_trans[2]]<cutoff[2] #high risk is 1
loc_3=exp[,int_trans[3]]>cutoff[3] #high risk is 1
loc_4=exp[,int_trans[4]]>cutoff[4] #high risk is 1
label_1<-as.matrix(round(loc_1))
label_2<-as.matrix(round(loc_2))
label_3<-as.matrix(round(loc_3))
label_4<-as.matrix(round(loc_4))
label_raw<-cbind(label_1,label_2,label_3,label_4)
colnames(label_raw)<-int_trans
rownames(label_raw)<-rownames(exp)
label_raw<-as.matrix(rowSums(label_raw))

label<-matrix(NA,length(label_raw),1)
index<-which(label_raw>=th)
label[index]=1
label[which(is.na(label))]=0
rownames(label)<-rownames(label_raw)

survival_time<-as.matrix(survival[,1])
status<-as.matrix(survival[,2])

survival_1<-as.matrix(survival_time[which(label==1),1])
survival_0<-as.matrix(survival_time[which(label==0),1])
num_1<-length(survival_1)
num_0<-length(survival_0)

summary(coxph(Surv(survival_time,status) ~ label))
surv_obj<-survfit(Surv(survival_time,status) ~ label)
surv_obj

path_out<-paste0(path_raw,"Output\\")
dir.create(path_out)
setwd(path_out)
pdf(file = "Output_Figure_3B_TCGA_KITC.pdf",width=5.8,height=6)
#pdf(file = "Output_Figure_3B_Japanese_KITC.pdf",width=5.8,height=6)
plot(surv_obj,col=c("black","red"), mark.time=T, cex=1.4,xlab="Time (year)",lty =1, ylab = "Survival Probability",las=1,cex.lab=1.4)
legend("bottomleft",legend =c(paste0("Low risk (n=",num_0,")"), paste0("High risk (n=",num_1,")")),col=c("black","red"),text.font=2,lty=c(1,1),lwd=2,bty="n",cex=1.4)
survtest <- survdiff(Surv(survival_time, status) ~ label)
log_rank_p<-1 - pchisq(survtest$chisq, 1) 
log_rank_p<-formatC(log_rank_p,format="e",digits=2)
mode(log_rank_p)<-"character"
legend("topright", legend =paste0("P=",log_rank_p),text.font=2,bty="n",cex=1.4)
dev.off()

##################################################################################
##################################################################################
#Figure 3C
#Figure 3C was drawn by Cytoscape software (Version 3.6.1).
#In file "Figure_3C_GoTerm_inform.txt", We took the column "Term_name" as node and "Label" as node attribute. The size of each node (go term) depends on the number of genes enriched in this go term (column "size").
#In the column "Label",the overlapped go terms identified from TCGA and Japanese cohort were denoted as 1. Among these overlapped go terms, the go terms which were associated with the four prognostic transcripts were denoted as 4. The go terms exclusively identified from TCGA cohort were denoted as 2.The go terms exclusively identified form Japanese cohort were denoted as 3.
##################################################################################
##################################################################################
#Fgiure 4A
#The Figure 4A was drawn manually based on the source data.
##################################################################################
##################################################################################
#Fgiure 4B
#The homology models were built using Schrodinger Suite software (Schrödinger/2019-3, LLC, New York, NY). 
##################################################################################
##################################################################################
#Figure 4C and 4D: 
#The Figure 4C and 4D were drawn manually based on the source data.
##################################################################################
##################################################################################
#Fgiure 4E: showing peptide intensity
rm(list=ls())
rm(list=ls())
rm(list=ls())

path_raw<-"E:\\PKM_data_code\\Source_data\\Figure_4E\\"
setwd(path_raw)

library(tidyverse) 
library(magrittr)
library(ggsci)
library(gridExtra)
library(gsubfn)
library(reshape2)
library(ggsci)
library(ggthemes)
library(stringr)

file <- "Figure_3E_peptides_whole_proteome.csv"
df <- read.csv(file, stringsAsFactors = FALSE) %>% as.tibble %>% na_if(., 0)

ENST00000389093 <- "MSKPHSEAGTAFIQTQQLHAAMADTFLEHMCRLDIDSPPITARNTGIICTIGPASRSVELKKGATLKITLDNAYMEKCDENILWLDYKNICKVVEVGSKIYVDDGLISLQVKQKGADFLVTEVENGGSLGSKKGVNLPGAAVDLPAVSEKDIQDLKFGVEQDVDMVFASFIRKASDVHEVRKVLGEKGKNIKIISKIENHEGVRRFDEILEASDGIMVARGDLGIEIPAEKVFLAQKMMIGRCNRAGKPVICATQMLESMIKKPRPTRAEGSDVANAVLDGADCIMLSGETAKGDYPLEAVRMQHLIAREAEAAIYHLQLFEELRRLAPITSDPTEATAVGAVEASFKCCSGAIIVLTKSGRSAHQVARYRPRAPIIAVTRNPQTARQAHLYRGIFPVLCKDPVQEAWAEDVDLRVNFAMNVGKARGFFKKGDVVIVLTGWRPGSGFTNTMRVVPVP"
ENST00000568883 <- "MSKPHSEAGTAFIQTQQLHAAMADTFLEHMCRLDIDSPPIKKGVNLPGAAVDLPAVSEKDIQDLKFGVEQDVDMVFASFIRKASDVHEVRKVLGEKGKNIKIISKIENHEGVRRFDEILEASDGIMVARGDLGIEIPAEKVFLAQKMMIGRCNRAGKPVICATQMLESMIKKPRPTRAEGSDVANAVLDGADCIMLSGETAKGDYPLEAVRMQHLIAREAEAAMFHRKLFEELVRASSHSTDLMEAMAMGSVEASYKCLAAALIVLTESGRSAHQVARYRPRAPIIAVTRNPQTARQAHLYRGIFPVLCKDPVQEAWAEDVDLRVNFAMNVGKARGFFKKGDVVIVLTGWRPGSGFTNTMRVVPVP"
length_49kDA <- nchar(ENST00000389093)
length_40kDA <- nchar(ENST00000568883)

df_883 <- df %>% mutate(ENST00000568883 = str_detect(ENST00000568883, regex(df$Sequence))) %>%
  filter(as.logical(.$ENST00000568883)) %>%
  mutate(no_proteins = lengths(strsplit(.$Proteins, ";")))

df_093 <- df %>% mutate(ENST00000389093 = str_detect(ENST00000389093, regex(df$Sequence))) %>%
  filter(as.logical(.$ENST00000389093)) %>%
  mutate(no_proteins = lengths(strsplit(.$Proteins, ";")))


df_883$start_883 <- apply(df_883, 1, function(x) length(strsplit(strapply(ENST00000568883, paste("(.*)", x[1], sep = ""), simplify = TRUE), "")[[1]]) + 1)
df_883$end_883 <- df_883$start_883 + lengths(strsplit(df_883$Sequence, "")) - 1

df_093$start_093 <- apply(df_093, 1, function(x) length(strsplit(strapply(ENST00000389093, paste("(.*)", x[1], sep = ""), simplify = TRUE), "")[[1]]) + 1)
df_093$end_093 <- df_093$start_093 + lengths(strsplit(df_093$Sequence, "")) - 1


### Top band - ENST00000389093 ###
df_t <-  df_093 %>% select(c(1:6, grep("_t", names(df_093)), 19:ncol(df_093))) %>% rowwise() %>% 
  mutate(median = median(c(Intensity.250ug_t, Intensity.rep1_t, Intensity.rep2_t, Intensity.rep3_t), na.rm = TRUE),
         mean = mean(c(Intensity.250ug_t, Intensity.rep1_t, Intensity.rep2_t, Intensity.rep3_t), na.rm = TRUE))

df_t$No_rep <-  df_t %>% 
  select(7:10) %>% 
  is.na %>% 
  `!` %>% rowSums
df_t <- df_t[df_t$No_rep != 0,]


### Bottom band - ENST00000568883 ###
df_b <-  df_883 %>% select(c(1:6, grep("_b", names(df_883)), 19:ncol(df_883))) %>% rowwise() %>% 
  mutate(median = median(c(Intensity.250ug_b, Intensity.rep1_b, Intensity.rep2_b, Intensity.rep3_b), na.rm = TRUE),
         mean = mean(c(Intensity.250ug_b, Intensity.rep1_b, Intensity.rep2_b, Intensity.rep3_b), na.rm = TRUE))

df_b$No_rep <-  df_b %>% 
  select(7:10) %>% 
  is.na %>% 
  `!` %>% rowSums
df_b <- df_b[df_b$No_rep != 0,]

##### Make plots #####
## Top ##
df_t_melted <- melt(df_t, id.vars = c("Sequence", "median"), measure.vars = c("start_093", "end_093"))
df_t_test <- df_t
df_t_test$median <- 500000
df_t_test_melted <- melt(df_t_test, id.vars = c("Sequence", "median"), measure.vars = c("start_093", "end_093"))

df_t_test_melted$col <- "b"
df_t_melted$col <- "a" 
df_t_melted_bound <- rbind(df_t_melted, df_t_test_melted)
df_t_melted_bound$name <- paste(df_t_melted_bound$Sequence, df_t_melted_bound$col)

df_t_line <- df_t %>% arrange(Start.position) %>% select(c(1, median, start_093, end_093))

df_empty <- matrix(NA, nrow = nrow(df_t_line), ncol = length_49kDA) %>% data.frame() %>% `colnames<-`(1:length_49kDA)

df_t_line <- cbind(df_t_line, df_empty)

i=1
for(i in 1:nrow(df_t_line)){
  df_t_line[i, c((df_t_line[i,"start_093"] + 4) : (df_t_line[i,"end_093"] + 4))] <- df_t_line$median[i]
  
}

df_t_line_melted <- colSums(df_t_line[,c(5:(length_49kDA+4))], na.rm = TRUE) %>% melt() %>% mutate(pos = c(1:length_49kDA)) %>% filter(value != 0) %>%
  `colnames<-`(c("median", "value")) %>% mutate(Sequence = "c", name = "c", col = "c", variable = "c", )

df_t_melted_bound_2 <- rbind(df_t_melted_bound, df_t_line_melted) 
df_t_melted_bound_3 <- df_t_melted_bound_2 %>% filter(col != "a")

line_t_III <- ggplot(df_t_melted_bound_3, aes(value, median, group = name, colour = col)) +
  geom_line(size = 1, linetype = 1, alpha = 1)  + geom_rangeframe(aes(colour = "black")) + theme_tufte() +
  ylab("log10 Intensity [IU]") + xlab("") + theme(plot.margin = unit(c(0,-2,0,0.5), "lines")) + scale_y_log10(breaks=c(1e6,1e8)) +
  theme(legend.position = "none") + scale_colour_manual(values=c("#3C5488B2", "black", "#00A087B2")) 
path_out<-paste0(path_raw,"Output\\")
dir.create(path_out)
setwd(path_out)
#pdf("Top_band_20190327_III.pdf",width=15,height=6)
pdf("Ouput_Figure_4E_093.pdf",width=15,height=6)
grid.arrange(line_t_III, ncol=2, nrow=2, widths=c(12, 2), heights=c(1, 0.1))
dev.off()

## Bottom ##
df_b_melted <- melt(df_b, id.vars = c("Sequence", "median"), measure.vars = c("start_883", "end_883"))
df_b_test <- df_b
df_b_test$median <- 50000
df_b_test_melted <- melt(df_b_test, id.vars = c("Sequence", "median"), measure.vars = c("start_883", "end_883"))

df_b_test_melted$col <- "b"
df_b_melted$col <- "a" 
df_b_melted_bound <- rbind(df_b_melted, df_b_test_melted)
df_b_melted_bound$name <- paste(df_b_melted_bound$Sequence, df_b_melted_bound$col)

df_b_line <- df_b %>% arrange(Start.position) %>% select(c(1, median, start_883, end_883))
df_empty_b <- matrix(NA, nrow = nrow(df_b_line), ncol = length_40kDA) %>% data.frame() %>% `colnames<-`(1:length_40kDA)
df_b_line <- cbind(df_b_line, df_empty_b)

for(i in 1:nrow(df_b_line)){
  df_b_line[i, c((df_b_line[i,"start_883"] + 4) : (df_b_line[i,"end_883"] + 4))] <- df_b_line$median[i]
}

df_b_line_melted <- colSums(df_b_line[,c(5:(length_40kDA+4))], na.rm = TRUE) %>% melt() %>% mutate(pos = c(1:length_40kDA)) %>% filter(value != 0) %>%
  `colnames<-`(c("median", "value")) %>% mutate(Sequence = "c", name = "c", col = "c", variable = "c", )

df_b_melted_bound_2 <- rbind(df_b_melted_bound, df_b_line_melted)
df_b_melted_bound_3 <- df_b_melted_bound_2 %>% filter(col != "a")

line_b_III <- ggplot(df_b_melted_bound_3, aes(value, median, group = name, colour = col)) +
  geom_line(size = 1, linetype = 1, alpha = 1)  + geom_rangeframe(aes(colour = "black")) + theme_tufte() +
  ylab("log10 Intensity [IU]") + xlab("") + theme(plot.margin = unit(c(0,-2,0,0.5), "lines")) + scale_y_log10(breaks=c(1e6,1e8)) +
  theme(legend.position = "none") + scale_colour_manual(values=c("#3C5488B2", "black", "#00A087B2")) 
setwd(path_out)
#pdf("Bottom_band_20190401_III.pdf",width=15,height=6)
pdf("Output_Figure_4E_883.pdf",width=15,height=6)
grid.arrange(line_b_III, ncol=2, nrow=2, widths=c(12, 2), heights=c(1, 0.1))
dev.off()
##################################################################################
##################################################################################
#Figure S1: Overlapped DEGs 
rm(list=ls())
rm(list=ls())
rm(list=ls())

library(ggplot2)

path_raw<-"E:\\PKM_data_code\\Source_data\\Fgiure_S1\\"
setwd(path_raw)
matrix<-as.matrix(read.csv("Figure_S1_TCGA_KIRC_DEGs_overlap.txt",header=T,row.names = 1,sep="\t"))
matrix<-matrix[,c("overlap","consis","not_consis")]
ratio_matrix<-cbind(matrix[,"consis"]/matrix[,"overlap"],matrix[,"not_consis"]/matrix[,"overlap"])
colnames(ratio_matrix)<-c("consis","not_consis")
#ratio_matrix
path_out<-paste0(path_raw,"Output\\")
dir.create(path_out)
for (i in 1:dim(ratio_matrix)[1]){
  df = data.frame("Classification" = c("Consistent ratio","Not consistent ratio"),
                  "ratio" = c(ratio_matrix[i,1],ratio_matrix[i,2]))
  pie = ggplot(df, aes(x="", y=ratio, fill=Classification)) + geom_bar(stat="identity", width=1)
  
  # Convert to pie (polar coordinates) and add labels
  pie = pie + coord_polar("y", start=0)
  
  # Add color scale (hex colors)
  pie = pie + scale_fill_manual(values=c("lightcoral", "lightskyblue")) 
  
  # Tidy up the theme
  pie = pie + theme_classic() + theme(axis.line = element_blank(),
                                      axis.text = element_blank(),
                                      axis.ticks = element_blank())
  setwd(path_out)
  ggsave(paste0("./",rownames(ratio_matrix)[i],".pdf"), width = 5, height = 5)
  
}

##################################################################################
##################################################################################
#Figure S2
rm(list=ls())
rm(list=ls())
rm(list=ls())

path_raw<-"E:\\PKM_data_code\\Source_data\\Figure_S2\\"
setwd(path_raw)
matrix_raw<-as.matrix(read.csv("Figure_S2_GoTerm_log10P.txt",header=T,sep="\t"))
#Figure_S1_GoTerm_log10P.txt: negtive log10 transformation of P values, then negtive transformation for the pathways enriched with down-regulated genes
matrix<-matrix_raw[,c(2,3,4,5)]
mode(matrix)<-"numeric"
rownames(matrix)<-matrix_raw[,1]

library(corrplot)
library(Hmisc)
library(gplots)
path_out<-paste0(path_raw,"Output\\")
dir.create(path_out)
setwd(path_out)
tiff(filename = "Output_Figure_S2_pathway_heatmap.tiff",width=1500,height=1200)
hr <- hclust(as.dist(1-cor(t(matrix), method="spearman")), method="ward.D2")
#hc <- hclust(as.dist(1-cor(matrix, method="spearman")), method="ward.D2") 
my_palete <- colorRampPalette(c("blue","white","red3"))(n=299)
col_breaks = c(seq(-1.3,-1,length=100),
               + seq(-0.999,1,length=100),
               + seq(1.009,1.3,length=100))  
#heatmap.2(log_p_value,Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),scale="none",density.info="density",trace="none",col=my_palete,breaks=col_breaks,margins =c(15,15),cexRow=1.5,cexCol=1.5,sepwidth=c(0.5,0.5),sepcolor="black")
heatmap.2(matrix,Rowv=as.dendrogram(hr), Colv=FALSE,scale="none",density.info="none",trace="none",col=my_palete,breaks=col_breaks,margins =c(15,60),cexRow=1.5,cexCol=1.5,sepwidth=c(0.5,0.5),sepcolor="black")
dev.off()

##################################################################################
##################################################################################
##################################################################################
#Supplementary Tables
##################################################################################
##################################################################################
#Table S1
#Table S1 shows the exton-intron structure, which was downloaded from the Ensembl website (Version 83, GENCODE Version 24).
##################################################################################
##################################################################################
#Table S2: mean expression value of PKM transcripts in TCGA cancer-types
rm(list=ls())
rm(list=ls())
rm(list=ls())

path_raw<-"E:\\PKM_data_code\\Source_data\\Table_S2\\"
setwd(path_raw)
folder_list<-as.matrix(read.csv("Cancer_list.txt",header=T,sep="\t"))
num<-NULL
mean_matrix<-NULL
for (i in 1:length(folder_list)){
  print(i)
  int_exp<-as.matrix(read.csv(paste0(folder_list[i],"_TPM_exp",".txt"),header=T,row.names = 1,sep="\t"))
  num_value<-dim(int_exp)[1]
  mean_value<-as.matrix(colMeans(int_exp))
  colnames(mean_value)<-folder_list[i]
  num<-rbind(num,num_value)
  mean_matrix<-cbind(mean_matrix,mean_value)
}
mean_matrix<-t(mean_matrix)
rownames(num)<-folder_list
colnames(num)="Size"
mean_matrix<-cbind(folder_list,num,mean_matrix)

path_out<-paste0(path_raw,"Output\\")
dir.create(path_out)
setwd(path_out)

write.table(mean_matrix,file="Output_Table_S2.txt",sep="\t",row.names=F,col.names=T,quote=F)

##################################################################################
##################################################################################
#Table S3: mean expression values of PKM transcripts in GTEx normal tissues
rm(list=ls())
rm(list=ls())
rm(list=ls())

path_raw<-"E:\\PKM_data_code\\Source_data\\Table_S3\\"
setwd(path_raw)
tissue_list<-as.matrix(read.csv("Tissue_list.txt",header=T,sep="\t"))
num<-NULL
mean_matrix<-NULL
for (i in 1:length(tissue_list)){
  print(i)
  int_exp<-as.matrix(read.csv(paste0("GTEx_",tissue_list[i],"_TPM_exp",".txt"),header=T,row.names = 1,sep="\t"))
  num_value<-dim(int_exp)[1]
  mean_value<-as.matrix(colMeans(int_exp))
  colnames(mean_value)<-tissue_list[i]
  num<-rbind(num,num_value)
  mean_matrix<-cbind(mean_matrix,mean_value)
}

mean_matrix<-t(mean_matrix)
rownames(num)<-tissue_list
colnames(num)="Size"
mean_matrix<-cbind(tissue_list,num,mean_matrix)

path_out<-paste0(path_raw,"Output\\")
dir.create(path_out)
setwd(path_out)
write.table(mean_matrix,file="Output_Table_S3.txt",sep="\t",row.names=F,col.names=T,quote=F)
##################################################################################
##################################################################################
#Table S4: Log-rank p value of PKM and its 14 transcripts in 25 cancer types
rm(list=ls())
rm(list=ls())
rm(list=ls())

library("grid")
library("xlsx")
library("XLConnect")
library("ggplot2")
library("ggvis")
library("rgl")
library("dplyr")
library("tidyr")
library("stringr")
library("lubridate")
require(gplots)
library("survival")
require("xlsx")
library("DESeq2")
library("biomaRt")
library("DESeq2")
library("piano")
library("Biobase")

setwd("E:\\PKM_data_code\\Code\\")
source('Cheng_toolbox_beta.R')

gene_list<-as.matrix(c("PKM","ENST00000335181","ENST00000561609","ENST00000389093","ENST00000568883","ENST00000319622","ENST00000562997","ENST00000568459"))
path_raw<-"E:\\PKM_data_code\\Source_data\\Table_S4\\"
setwd(path_raw)
folder_list<-as.matrix(read.csv("Cancer_list.txt",header=T,sep="\t"))

path_out<-paste0(path_raw,"Output\\")
dir.create(path_out)

p_end<-NULL
#cutoff_end<-NULL
#coef_end<-NULL
for (i in 1:length(folder_list)){
  print(i)
  cancerType<-as.matrix(folder_list[i,1])
  dataDir<-path_raw
  setwd(dataDir)
  TXT_name<-paste0(folder_list[i],"_TPM_exp_SurvivalData.txt")
  setwd(path_out)
  cancer_p<-NULL
  #cancer_cutoff<-NULL
  #cancer_coef<-NULL
  for (j in 1:length(gene_list)){
    output<-Cheng_generateSurvInputfromTCGA(gene_list[j,1],TXT_name,dataDir)
    result<-Cheng_generateKMplot(output,outFile=paste0(cancerType,"_",gene_list[j,1]))
    log_rank_p<-as.matrix(result$logRankP)
    #cut_off<-as.matrix(result$EXPcut)
    #coef<-as.matrix(result$coef)
    colnames(log_rank_p)<-gene_list[j,1]
    #colnames(cut_off)<-gene_list[j,1]
    #colnames(coef)<-gene_list[j,1]
    cancer_p<-cbind(cancer_p,log_rank_p)
    #cancer_cutoff<-cbind(cancer_cutoff,cut_off)
    #cancer_coef<-cbind(cancer_coef,coef)
  }

  p_end<-rbind(p_end,cancer_p)
  #cutoff_end<-rbind(cutoff_end,cancer_cutoff)
  #coef_end<-rbind(coef_end,cancer_coef)
}

rownames(p_end)<-folder_list
#rownames(cutoff_end)<-folder_list
#rownames(coef_end)<-folder_list
colnames(p_end)<-gene_list
#colnames(cutoff_end)<-gene_list
#colnames(coef_end)<-gene_list
p_end<-cbind(folder_list,p_end)
setwd(path_out)
write.table(p_end,file="Output_Table_S4.txt",sep="\t",row.names=F,col.names=T,quote=F)

##################################################################################
##################################################################################
#Table S5: The enriched of GO terms with DEGs for all transcripts in 25 cancer types
rm(list=ls())
rm(list=ls())
rm(list=ls())

int_trans<-"ENST00000568883"#"ENST00000335181","ENST00000561609","ENST00000389093","ENST00000568883")
path_raw<-"E:\\PKM_data_code\\Source_data\\Table_S5\\"
path_int_trans<-paste0(path_raw,int_trans,"\\")
setwd(path_int_trans)
goid_all<-as.matrix(read.csv("GoTerm_id.txt",header=F,sep="\t"))
p_matrix_up<-as.matrix(read.csv("Up_GoTerm_p_matrix.txt",header=T,sep="\t"))
rownames(p_matrix_up)<-p_matrix_up[,1]
p_matrix_up<-p_matrix_up[,-1]
mode(p_matrix_up)<-"numeric"

p_matrix_down<-as.matrix(read.csv("Down_GoTerm_p_matrix.txt",header=T,sep="\t"))
rownames(p_matrix_down)<-p_matrix_down[,1]
p_matrix_down<-p_matrix_down[,-1]
mode(p_matrix_down)<-"numeric"


binary_up=p_matrix_up
binary_up[which(!is.na(p_matrix_up))]=1
binary_up[which(is.na(p_matrix_up))]=0
row_sum_up<-as.matrix(rowSums(binary_up))#the frequency of each path in cancers

binary_down=p_matrix_down
binary_down[which(!is.na(p_matrix_down))]=1
binary_down[which(is.na(p_matrix_down))]=0
row_sum_down<-as.matrix(rowSums(binary_down))

binary_all=binary_up+binary_down
binary_all[which(binary_all==2)]=1#the pathwy is simultaneously up- and down-regulated 
Generality<-as.matrix(rowSums(binary_all))
Direction<-as.matrix(row_sum_up-row_sum_down)

drop_index<-which(Generality==0)
Generality<-as.matrix(Generality[-drop_index,])
Direction<-as.matrix(Direction[-drop_index,])
goid_all<-as.matrix(goid_all[-drop_index,])
go_inform<-cbind(goid_all,rownames(Direction),Direction,Generality)
colnames(go_inform)<-c("GO_ID","term_name","Direction","Generality")

path_out<-"E:\\PKM_data_code\\Source_data\\Table_S5\\Output\\"
setwd(path_out)
#save(file=paste0(int_trans,"_GoTerm_Generality_direction.Rdata"),Generality,Direction,goid_all)
write.table(go_inform,file=paste0(int_trans,"_GoTerm_Generality_direction.txt"),sep="\t",row.names=F,col.names=T,quote=F)

##################################################################################
##################################################################################
#Table S6-9: DEGs overlapping
rm(list=ls())
rm(list=ls())
rm(list=ls())

setwd("E:\\PKM_data_code\\Code\\")
source('DEGs_DEseq_table_overlap_2.R')

#int_transcript<-as.matrix(c("ENST00000335181","ENST00000561609","ENST00000389093","ENST00000568883"))#TCGA-KIRC
#path_raw<-"E:\\PKM_data_code\\Source_data\\Table_S6\\TCGA_KIRC_DEGs\\"#TCGA-KIRC
#int_transcript<-as.matrix(c("ENST00000335181","ENST00000389093","ENST00000568883"))#TCGA-CESC
#path_raw<-"E:\\PKM_data_code\\Source_data\\Table_S7\\TCGA_CESC_DEGs\\"#TCGA-CESC
#int_transcript<-as.matrix(c("ENST00000335181","ENST00000389093","ENST00000568883"))#TCGA-PAAD
#path_raw<-"E:\\PKM_data_code\\Source_data\\Table_S8\\TCGA_PAAD_DEGs\\"#TCGA-PAAD
int_transcript<-as.matrix(c("ENST00000335181","ENST00000389093","ENST00000568883"))#TCGA-BRCA
path_raw<-"E:\\PKM_data_code\\Source_data\\Table_S9\\TCGA_BRCA_DEGs\\"#TCGA-BRCA

fdr_th<-0.00001 #for extracting DEGs
index_p<-t(combn(dim(int_transcript)[1],2))
int_pair<-cbind(int_transcript[index_p[,1],1],int_transcript[index_p[,2],1])
path_out<-paste0(path_raw,"Output\\")
dir.create(path_out)

over_result_common_end<-NULL
for (i in 1:dim(int_pair)[1]){
  project_1<-as.matrix(int_pair[i,1])
  project_2<-as.matrix(int_pair[i,2])
  
  path_1<-paste0(path_raw,project_1,"\\") 
  setwd(path_1)
  deg_DEseq_1<-as.matrix(read.csv("DEseq_DEGs_result.txt",header=T,sep="\t"))
  
  path_2<-paste0(path_raw,project_2,"\\") 
  setwd(path_2)
  deg_DEseq_2<-as.matrix(read.csv("DEseq_DEGs_result.txt",header=T,sep="\t"))
  
  output<-DEGs_DEseq_table_overlap_2(deg_DEseq_1,deg_DEseq_2,fdr_th)
  over_result=output$over_result
  file_over_result<-paste0(project_1," VS ",project_2)
  over_result_common<-t(as.matrix(over_result[2,]))
  rownames(over_result_common)<-file_over_result
  over_result_common_end<-rbind(over_result_common_end,over_result_common)
}

comparison<-as.matrix(rownames(over_result_common_end))
over_result_common_end<-cbind(comparison,over_result_common_end)
over_result_common_end<-over_result_common_end[,-c(4,8)]
colnames(over_result_common_end)<-c("Comparison",	"Number of DEGs for transcript 1"	,"Number of DEGs for transcript 2","Overlaps",	"Consistent DEGs",	"Concordance ratio",	"Consistent up-regulated DEGs",	"Consistent down-regulated DEGs",	"P value")
setwd(path_out)
#write.table(over_result_common_end,file="Output_Table_S6.txt",sep="\t",row.names=F,col.names=T,quote=F) 
#write.table(over_result_common_end,file="Output_Table_S7.txt",sep="\t",row.names=F,col.names=T,quote=F) 
#write.table(over_result_common_end,file="Output_Table_S8.txt",sep="\t",row.names=F,col.names=T,quote=F) 
#write.table(over_result_common_end,file="Output_Table_S8.txt",sep="\t",row.names=F,col.names=T,quote=F) 
write.table(over_result_common_end,file="Output_Table_S9.txt",sep="\t",row.names=F,col.names=T,quote=F) 

##################################################################################
##################################################################################
#Table S10: DEGs overlapping
rm(list=ls())
rm(list=ls())
rm(list=ls())

setwd("E:\\PKM_data_code\\Code\\")
source('DEGs_DEseq_table_overlap_2.R')

int_transcript<-as.matrix(c("ENST00000335181","ENST00000389093"))#TCGA-COAD
path_raw<-"E:\\PKM_data_code\\Source_data\\Table_S10\\TCGA_COAD_DEGs\\"#TCGA-COAD

fdr_th<-0.00001 #for extracting DEGs
index_p<-t(combn(dim(int_transcript)[1],2))
int_pair<-cbind(int_transcript[index_p[,1],1],int_transcript[index_p[,2],1])
path_out<-paste0(path_raw,"Output\\")
dir.create(path_out)

over_result_common_end<-NULL
for (i in 1:dim(int_pair)[1]){
  project_1<-as.matrix(int_pair[i,1])
  project_2<-as.matrix(int_pair[i,2])
  
  path_1<-paste0(path_raw,project_1,"\\") 
  setwd(path_1)
  deg_DEseq_1<-as.matrix(read.csv("DEseq_DEGs_result.txt",header=T,sep="\t"))
  
  path_2<-paste0(path_raw,project_2,"\\") 
  setwd(path_2)
  deg_DEseq_2<-as.matrix(read.csv("DEseq_DEGs_result.txt",header=T,sep="\t"))
  
  output<-DEGs_DEseq_table_overlap_2(deg_DEseq_1,deg_DEseq_2,fdr_th)
  over_result=output$over_result
  file_over_result<-paste0(project_1," VS ",project_2)
  over_result_common<-t(as.matrix(over_result[2,]))
  rownames(over_result_common)<-file_over_result
  over_result_common_end<-rbind(over_result_common_end,over_result_common)
}

comparison<-as.matrix(rownames(over_result_common_end))
over_result_common_end<-cbind(comparison,over_result_common_end)
over_result_common_end<-t(as.matrix(over_result_common_end[,-c(4,8)]))

colnames(over_result_common_end)<-c("Comparison",	"Number of DEGs for transcript 1"	,"Number of DEGs for transcript 2","Overlaps",	"Consistent DEGs",	"Concordance ratio",	"Consistent up-regulated DEGs",	"Consistent down-regulated DEGs",	"P value")
setwd(path_out)
write.table(over_result_common_end,file="Output_Table_S10.txt",sep="\t",row.names=F,col.names=T,quote=F)
##################################################################################
##################################################################################
#Table S11: cox analysis of four transcripts in TCGA KIRC
rm(list=ls())
rm(list=ls())
rm(list=ls())
library(survival)
path_raw<-"E:\\PKM_data_code\\Source_data\\Table_S11\\"
setwd(path_raw)
exp_trans<-as.matrix(read.csv("Table_S11_TCGA_KIRC_trans_TPM.txt",header=T,row.names = 1,sep="\t"))
survival<-as.matrix(read.csv("Table_S11_TCGA_KIRC_SurvivalData.txt",header=T,row.names = 1,sep="\t"))
est_622<-as.matrix(exp_trans[,"ENST00000319622"])
est_181<-as.matrix(exp_trans[,"ENST00000335181"])
est_609<-as.matrix(exp_trans[,"ENST00000561609"])
est_093<-as.matrix(exp_trans[,"ENST00000389093"])
est_883<-as.matrix(exp_trans[,"ENST00000568883"])
summary(coxph(Surv(as.matrix(survival[,1]), as.matrix(survival[,2])) ~ est_181))
summary(coxph(Surv(as.matrix(survival[,1]), as.matrix(survival[,2])) ~ est_609))
summary(coxph(Surv(as.matrix(survival[,1]), as.matrix(survival[,2])) ~ est_093))
summary(coxph(Surv(as.matrix(survival[,1]), as.matrix(survival[,2])) ~ est_883))

summary(coxph(Surv(as.matrix(survival[,1]), as.matrix(survival[,2])) ~ est_609+est_181))
summary(coxph(Surv(as.matrix(survival[,1]), as.matrix(survival[,2])) ~ est_093+est_181))
summary(coxph(Surv(as.matrix(survival[,1]), as.matrix(survival[,2])) ~ est_883+est_181))

summary(coxph(Surv(as.matrix(survival[,1]), as.matrix(survival[,2])) ~ est_181+est_622))
summary(coxph(Surv(as.matrix(survival[,1]), as.matrix(survival[,2])) ~ est_609+est_622))
summary(coxph(Surv(as.matrix(survival[,1]), as.matrix(survival[,2])) ~ est_093+est_622))
summary(coxph(Surv(as.matrix(survival[,1]), as.matrix(survival[,2])) ~ est_883+est_622))

##################################################################################
##################################################################################
#Table S12: DEGs overlapping
rm(list=ls())
rm(list=ls())
rm(list=ls())

quantile_th<-0.2
setwd("E:\\PKM_data_code\\Code\\")
source('DEGs_DEseq_table_overlap_quantile_th.R')

int_trans<-as.matrix(c("ENST00000335181","ENST00000561609","ENST00000389093","ENST00000568883"))
path_raw<-"E:\\PKM_data_code\\Source_data\\Table_S12\\"
path_out<-paste0(path_raw,"Output\\")
dir.create(path_out)

over_result_common_end<-NULL
for (i in 1:length(int_trans)){
  path_1<-paste0(path_raw,"TCGA_KIRC_DEGs\\",int_trans[i],"\\")
  setwd(path_1)
  deg_DEseq_1<-as.matrix(read.csv("DEseq_DEGs_result.txt",header=T,sep="\t"))
  
  path_2<-paste0(path_raw,"Japanese_KIRC_DEGs\\",int_trans[i],"\\")
  setwd(path_2)
  deg_DEseq_2<-as.matrix(read.csv("DEseq_DEGs_result.txt",header=T,sep="\t"))
  
  output_end<-DEGs_DEseq_table_overlap_quantile_th(deg_DEseq_1,deg_DEseq_2,quantile_th)
  over_result=output_end$over_result
  file_over_result<-int_trans[i]
  over_result_common<-t(as.matrix(over_result[2,]))
  rownames(over_result_common)<-file_over_result
  over_result_common_end<-rbind(over_result_common_end,over_result_common)
}
comparison<-as.matrix(rownames(over_result_common_end))
over_result_common_end<-cbind(comparison,over_result_common_end)
over_result_common_end<-over_result_common_end[,-c(4,8)]

colnames(over_result_common_end)<-c("Comparison",	"Number of DEGs for transcript 1"	,"Number of DEGs for transcript 2","Overlaps",	"Consistent DEGs",	"Concordance ratio",	"Consistent up-regulated DEGs",	"Consistent down-regulated DEGs",	"P value")
setwd(path_out)
write.table(over_result_common_end,file="Output_Table_S12.txt",sep="\t",row.names=F,col.names=T,quote=F) 

##################################################################################
##################################################################################
#Table S13: Go enrichment
rm(list=ls())
rm(list=ls())
rm(list=ls())

library(clusterProfiler)
library(dplyr)
library(tidyr)
library(DOSE)
library(GO.db)
library(org.Hs.eg.db)
library(GSEABase)

fdr_goterm<-0.00001 #go enrichment cutoff
setwd("E:\\PKM_data_code\\Code\\")
source('deg_GoTerm_clusterProfiler.R')
path_raw<-"E:\\PKM_data_code\\Source_data\\Table_S13\\"
setwd(path_raw)
background<-as.matrix(read.csv("Background_gene.txt",header=F,sep="\t"))

#deg_DEGseq<-as.matrix(read.csv("TCGA_DEGs_quantile_th_0.2.txt",header=F,sep="\t"))
deg_DEGseq<-as.matrix(read.csv("Japanese_DEGs_quantile_th_0.2.txt",header=F,sep="\t"))
deg_up<-as.matrix(deg_DEGseq[which(as.numeric(deg_DEGseq[,3])=="1"),2])#extract the gene symbol
pathway_up=enrichGO(gene=deg_up,OrgDb='org.Hs.eg.db',keyType = "SYMBOL",ont="BP",universe = background ,pAdjustMethod = "BH",qvalueCutoff=0.05)
path_table_up<-deg_GoTerm_clusterProfiler(pathway_up)

p.adjust.up<-as.matrix(path_table_up[,"p.adjust"])
mode(p.adjust.up)="numeric"
sig_path_up<-as.matrix(path_table_up[which(p.adjust.up<fdr_goterm),])

path_out<-paste0(path_raw,"Output\\")
dir.create(path_out)
setwd(path_out)

#write.table(sig_path_up,file="TCGA_GoTerms_up.txt",sep="\t",row.names=F,col.names=T,quote=F)
write.table(sig_path_up,file="Janpanese_GoTerms_up.txt",sep="\t",row.names=F,col.names=T,quote=F)
##################################################################################
##################################################################################
#Table S14: Alignment of amino acid sequences of PKM transcripts
#The amino acid sequences of different transcripts were aligned using Uniprot website.
##################################################################################
##################################################################################
#Table S15: Relative protein level
#The quantification of protein is based on the western blots image.







