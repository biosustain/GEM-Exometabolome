## RESULTS ANALYSIS - LOOPLESS+COMPACTED+ROOM (ALL REACTIONS)|| JUNE 28  2022
#ERs & SL for test 5 - There is no SLs in this test.  

# Libraries ####
library(readr);library(tidyverse)
################################## SINGLE LETHAL REACTIONS ############################
newERS <- ERs_ROOMALL_28Jun <- read_delim("Desktop/ERs_ROOMALL_28Jun.csv", delim = ";", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)


colnames(newERS)<-c('RXNS','RXNS','RXNS','RXNS','RXNS','RXNS','RXNS','RXNS','RXNS') #Rename to combine
A2<-newERS[,1] #split into vector
B2<-newERS[,2]
C2<-newERS[,3]
D2<-newERS[,4]
E2<-newERS[,5]
G2<-newERS[,6]
H2<-newERS[,7]
I2<-newERS[,8]
J2<-newERS[,9]
newERS<-rbind(A2,B2,C2,D2,E2,G2,H2,I2) #join into single one
newERS<- newERS %>% drop_na() #:) 150 Total
#
## Values for contingency table

TPlist<-newERS%>%filter(newERS$RXNS %in% LethalRxns_Exp$RXNS)
TP<-44

FPlist<-newERS%>%filter(newERS$RXNS %!in% LethalRxns_Exp$RXNS)
FP<-106

TNlist<- AllRxns %>% filter(AllRxns$X1 %!in% newERS$RXNS) #all model reaction not in list of predictions and not in exp list
TNlist <- TNlist%>% filter(TNlist$X1 %!in%  LethalRxns_Exp$RXNS) #2210
TN<- 2216 

FP<-2583-44-106-2216

perc=(2216+44)/2583*100 # 87.49516 % <- Prediction ratio
##
## Statistics ## 
#Order or the contingency table / matrix -> (TP/FP/FN/TN)
ContMatrix_Compact<- as.data.frame(matrix(c(44,106,217,2216),nrow=2)) #make matrix

fisherCompact<-fisher.test(ContMatrix_Compact) #Fisher's Exact test
fisherCompact$p.value #7.413246e-12<- P-value


##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##################### RESULTS WITH THE CORRECTION FOR SOME REACTIONS (ERs) ########

Gene_correctedLethalRxns <- read_csv("Desktop/ERs_SLs/ExperimentalData/Gene_correctedLethalRxns.csv")
TPlist<-newERS%>%filter(newERS$RXNS %in% Gene_correctedLethalRxns$Rxns) #43
TP<-43

FPlist<-newERS%>%filter(newERS$RXNS %!in% Gene_correctedLethalRxns$Rxns)
FP<-107

TNlist<- AllRxns %>% filter(AllRxns$X1 %!in% newERS$RXNS) #all model reaction not in list of predictions and not in exp list
TNlist <- TNlist%>% filter(TNlist$X1 %!in%  Gene_correctedLethalRxns$Rxns) #2210
TN<-2221

FN<-2583-2221-43-107 #212

perc=(2221+43)/2583*100 # 87.65002 % <- Prediction ratio
ContMatrix_Compact<- as.data.frame(matrix(c(43,107,212,2221),nrow=2)) #make matrix

fisherCompact<-fisher.test(ContMatrix_Compact) #Fisher's Exact test
fisherCompact$p.value #1.382327e-11<- P-value








