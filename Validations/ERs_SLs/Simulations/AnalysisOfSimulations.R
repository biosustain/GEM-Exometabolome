## RESULTS ANALYSIS - LOOPLESS+COMPACTED+PARALLEL SL+ROOM & ORIGINAL MODEL + FASTSL|| JUNE 28  2022
#ERs & SL for tests 1 and 4 - original and gene corrected versions. 


# Libraries ####
library(readr);library(tidyverse)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~COMPACTED+LOOPLESS + ROOM  #~~~~~~
################################## SINGLE LETHAL REACTIONS ############################
#
##Load data
LethalRxns_Pipe<- read_delim("/Users/isabelmontejano/DTU-final/ERs_SLs/ROOM+Compact+Loopless/ERs_Compacted_Loopless_28jun.csv", 
                                           delim = ";", escape_double = FALSE, col_names = FALSE, 
                                           trim_ws = TRUE)

LethalRxns_Exp_7Jun <- read_csv("/Users/isabelmontejano/DTU-final/ERs_SLs/ExperimentalData/LethalRxns_Exp_7Jun.csv") #"Experimental" lethal reaction - from lethal genes

## ORDER RXNS from MATLAB table into a single row with all the reactions ##
#LethalRxns_Pipe<-Jsl_Compacted_Loopless_28jun 
colnames(LethalRxns_Pipe)<-c('RXNS','RXNS','RXNS','RXNS','RXNS','RXNS','RXNS','RXNS','RXNS') #Rename to combine
A<-LethalRxns_Pipe[,1] #split into vector
B<-LethalRxns_Pipe[,2]
C<-LethalRxns_Pipe[,3]
D<-LethalRxns_Pipe[,4]
E<-LethalRxns_Pipe[,5]
G<-LethalRxns_Pipe[,6]
H<-LethalRxns_Pipe[,7]
I<-LethalRxns_Pipe[,8]
J<-LethalRxns_Pipe[,9]
LethalRxns_Pipe<-rbind(A,B,C,D,E,G,H,I,J) #join into single one
LethalRxns_Pipe_Final %>% drop_na() #:) s

LethalRxns_Pipe<-unique(LethalRxns_Pipe) #Predicted from pipeline w/ compact model - check for duplicates
#LethalRxns_Pipe$RXNS<-as.character(LethalRxns_Pipe$RXNS) #as character to compare

LethalRxns_SL<-LethalRxns_SL_9Jun #Normal method Rxns
LethalRxns_SL$X1<- as.character(LethalRxns_SL_9Jun$X1)#as character to compare
colnames(LethalRxns_SL)<-'RXNS'

LethalRxns_Exp<-LethalRxns_Exp_7Jun #Experimental rxns
LethalRxns_Exp$Rxns<-as.character(LethalRxns_Exp$Rxns)#as character to compare
colnames(LethalRxns_Exp)<-'RXNS'

##
## Values for contingency matrix ## -- 160 total predictions
TPList<- LethalRxns_Pipe %>% filter(LethalRxns_Pipe$RXNS %in% LethalRxns_Exp$RXNS)
 TP<-47  #47 ,- true positive
FPlist<- LethalRxns_Pipe %>% filter(LethalRxns_Pipe$RXNS %!in% LethalRxns_Exp$RXNS)#113 <- false positive
FP<-113
TNlist<- AllRxns %>% filter(AllRxns$X1 %!in% LethalRxns_Pipe$RXNS) #all model reaction not in list of predictions and not in exp list
TNlist <- TNlist%>% filter(TNlist$X1 %!in%  LethalRxns_Exp$RXNS) #2210
  TN<- 2209 #<- True negative
  FN<- 214 #<- False negative ---- All - (TN,+TP+FP)
perc=(2209+47)/2583*100 # 87.3403 % <- Prediction ratio
##
## Statistics ## 
#Order or the contingency table / matrix -> (TP/FP/FN/TN)
  ContMatrix_Compact<- as.data.frame(matrix(c(47,113,214,2209),nrow=2)) #make matrix

  
  fisherCompact<-fisher.test(ContMatrix_Compact) #Fisher's Exact test
  fisherCompact$p.value #1.18232e-12<- P-value
  ################################## SYNTHETIC LETHAL REACTIONS ############################
  #~~~~~~~~~~~~~~
  # Load data ##
Synthetic_lethalRxns_Exp <- read_csv("/Users/isabelmontejano/DTU-final/ERs_SLs/ExperimentalData/ExperimentalData/Synthetic lethalRxns_Exp_7Jun.csv")
Doubles1 <- read_delim("DTU-final/ERs_SLs/ROOM+Compact+Loopless/SL1_CompactedLoopless_28Jun.csv", delim = ";", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
Doubles2 <- read_delim("DTU-final/ERs_SLs/ROOM+Compact+Loopless/SL2_CompactedLoopless_28Jun.csv", delim = ";", escape_double = FALSE, col_names = FALSE,  trim_ws = TRUE) 
#
##Values for contingency matrix 
TPSLlist<- AllCombo %>% filter((AllCombo$Var1 %in% Synthetic_lethalRxns_Exp$X1 & AllCombo$Var2 %in% Synthetic_lethalRxns_Exp$X2) |
                                 (AllCombo$Var1 %in% Gene_correctedLethalRxns$Rxns & AllCombo$Var2 %in% Gene_correctedLethalRxns$Rxns)) 
TP<- 12 

FPSLlist<- AllCombo %>% filter((AllCombo$Var1 %!in% Synthetic_lethalRxns_Exp$X1 | AllCombo$Var2 %!in% Synthetic_lethalRxns_Exp$X2) &
                                 (AllCombo$Var1 %!in% Gene_correctedLethalRxns$Rxns | AllCombo$Var2 %!in% Gene_correctedLethalRxns$Rxns)) 
FP<- 49 


FNlist<-Synthetic_lethalRxns_Exp %>% filter(Synthetic_lethalRxns_Exp$X1 %!in% AllCombo$Var1 & Synthetic_lethalRxns_Exp$X2 %!in% AllCombo$Var2)
FN<-202

TN<-possibleCombinations-TP-FP-FN #16894
ContMatrix_Normal<- as.data.frame(matrix(c(TP,FN,FP,TN),nrow=2)) #make matrix
ContMatrix_Normal
percNorm<- ((TN+TP)/possibleCombinations)*100 
percNorm
##
## Statistics ## 
#Order or the contingency table / matrix -> (TP/FP/FN/TN)

fisherNormal<-fisher.test(ContMatrix_Normal) #Fisher's Exact test
fisherNormal$p.value #5.358634e-13 <- P-value


##################### RESULTS WITH THE CORRECTION FOR SOME REACTIONS (ERs & SLs) ########
#~~~~~~~~~~~~~~~~~~ Essential reactions

TPNormlist<- LethalRxns_Pipe %>% filter(LethalRxns_Pipe$RXNS %in% Gene_correctedLethalRxns$Rxns)  
TP<-dim(TPNormlist)[1]
FPNorm<- LethalRxns_Pipe %>% filter(LethalRxns_Pipe$RXNS %!in% Gene_correctedLethalRxns$Rxns)  
FP<-dim(FPNorm)[1]

TNNormlist<- AllRxns %>% filter(AllRxns$X1 %!in% LethalRxns_Pipe$RXNS) 
TNNormlist <- TNNormlist %>% filter(TNNormlist$X1 %!in%  Gene_correctedLethalRxns$Rxns) 



TN<- dim(TNNormlist)[1]
FN<- 2583-TP-FP-TN 
ContMatrix_Normal<- as.data.frame(matrix(c(TP,FN,FP,TN),nrow=2))
ContMatrix_Normal
percNorm<- ((TN+TP)/2583)*100 #87.49516 <- Prediction ratio
##
## Statistics ## 
fisherNormal<-fisher.test(ContMatrix_Normal) #Fisher's Exact test
fisherNormal$p.value #2.091731e-12 <- P-value

#~~~~~~~~~~~~~~~~ Synthetic lethals
colnames(AllCombo)<-c('Var1','Var2')
#
##Values for contingency matrix 
TPSLlist<- AllCombo %>% filter((AllCombo$Var1 %in% Gene_CorrectedSynthetic_lethal$X1 & AllCombo$Var2 %in% Gene_CorrectedSynthetic_lethal$X2) |
                                 (AllCombo$Var1 %in% Gene_correctedLethalRxns$Rxns & AllCombo$Var2 %in% Gene_correctedLethalRxns$Rxns)) 
TP<- 13 

FPSLlist<- AllCombo %>% filter((AllCombo$Var1 %!in% Gene_CorrectedSynthetic_lethal$X1 | AllCombo$Var2 %!in% Gene_CorrectedSynthetic_lethal$X2) &
                                 (AllCombo$Var1 %!in% Gene_correctedLethalRxns$Rxns | AllCombo$Var2 %!in% Gene_correctedLethalRxns$Rxns)) 
FP<- 48 
FNlist<-Gene_CorrectedSynthetic_lethal %>% filter(Gene_CorrectedSynthetic_lethal$X1 %!in% AllCombo$Var1 & Gene_CorrectedSynthetic_lethal$X2 %!in% AllCombo$Var2)
FN<-203

TN<-possibleCombinations-TP-FP-FN #16894
ContMatrix_Normal<- as.data.frame(matrix(c(TP,FN,FP,TN),nrow=2)) #make matrix
ContMatrix_Normal
percNorm<- ((TN+TP)/possibleCombinations)*100 #98.537 <- Prediction ratio
percNorm
##
## Statistics ## 
#Order or the contingency table / matrix -> (TP/FP/FN/TN)

fisherNormal<-fisher.test(ContMatrix_Normal) #Fisher's Exact test
fisherNormal$p.value #5.358634e-13 <- P-value
  
  
  
  
  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ORIGINAL MODEL + FASTSL #~~~~~~
################################## SINGLE LETHAL REACTIONS ############################
#
## load data 
LethalRxns_SL_9Jun - read_csv("/Users/isabelmontejano/DTU-final/ERs_SLs/ROOM+Compact+Loopless/ERs_Original_28Jun.csv", col_names = FALSE)#Predicted lethal reactions  - normal model 

  LethalRxns_SL <- read_csv("Desktop/ERs_Original_28Jun.csv", 
                                 col_names = FALSE)
  colnames(LethalRxns_SL)<-'RXNS'
## Values for contingency matrix  ##
TPNormlist_O<- LethalRxns_SL %>% filter(LethalRxns_SL$RXNS %in% LethalRxns_Exp$RXNS)  #94 <-true positive
  TPNormlist<-94
  FPNorm_O<- LethalRxns_SL %>% filter(LethalRxns_SL$RXNS %!in% LethalRxns_Exp$RXNS)  #175 <- false positive
TNNormlist_O<- AllRxns %>% filter(AllRxns$X1 %!in% LethalRxns_SL$RXNS) 
TNNormlist_O <- TNNormlist_O %>% filter(TNNormlist_O$X1 %!in%  LethalRxns_Exp$RXNS) #2147 <- true negative
  TNNorm<- 2147
FN<- 167 #<- False negative
FNlist_O<-AllRxns %>% filter(AllRxns$X1 %!in% TPNormlist_O$RXNS) 
FNlist_O<-FNlist_O %>% filter(FNlist_O$X1 %!in% TNNormlist_O$X1) 
FNlist_O<-FNlist_O %>% filter(FNlist_O$X1 %!in% FPNorm_O$RXNS)


percNorm<- ((2147+94)/2583)*100 #86.759 <- Prediction ratio
##
## Statistics ## 
#Order or the contingency table / matrix -> (TP/FP/FN/TN)
ContMatrix_Normal<- as.data.frame(matrix(c(94,175,167,2147),nrow=2)) #make matrix
    row.names(ContMatrix_Normal) = c("Positive", "Negative")
    colnames(ContMatrix_Normal) = c("Positive", "Negative")
    mosaicplot(ContMatrix_Normal,main='Contincengy table - Non-compact model')

fisherNormal<-fisher.test(ContMatrix_Normal) #Fisher's Exact test
fisherNormal$p.value #5.241556e-33 <- P-value

################################## SYNTHETIC LETHALS ############################
#SL_normal <- read_delim("/Users/isabelmontejano/DTU-final/ERs_SLs/ROOM+Compact+Loopless/SL_Original+FASTSL.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)

SL_normal<- read_delim("/Users/isabelmontejano/DTU-final/ERs_SLs/ROOM+Compact+Loopless/SLs_Original_28Jun.csv", delim = ";", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
colnames(SL_normal)<- c('RXN1','RXN2')


## Values for contingency matrix ##
TPSLlist<- SL_normal %>% filter((SL_normal$Var1 %in% Synthetic_lethalRxns_Exp$X1 & SL_normal$Var2 %in% Synthetic_lethalRxns_Exp$X2) |
                                  (SL_normal$Var1 %in% Gene_correctedLethalRxns$Rxns & SL_normal$Var2 %in% Gene_correctedLethalRxns$Rxns)) 
TP<- 25 

FPSLlist<- SL_normal %>% filter((SL_normal$Var1 %!in% Synthetic_lethalRxns_Exp$X1 | SL_normal$Var2 %!in% Synthetic_lethalRxns_Exp$X2) &
                                  (SL_normal$Var1 %!in% Gene_correctedLethalRxns$Rxns | SL_normal$Var2 %!in% Gene_correctedLethalRxns$Rxns)) 
FP<- 236 

FNlist<-Synthetic_lethalRxns_Exp %>% filter(Synthetic_lethalRxns_Exp$X1 %!in% SL_normal$Var1 & Synthetic_lethalRxns_Exp$X2 %!in% SL_normal$Var2)
FN<-136

TN<-possibleCombinations-25-136-236 #
ContMatrix_Normal<- as.data.frame(matrix(c(TP,FN,FP,TN),nrow=2)) #make matrix
ContMatrix_Normal
percNorm<- ((TN+TP)/possibleCombinations)*100 #98.41464 <- Prediction ratio

##
## Statistics ## 
#Order or the contingency table / matrix -> (TP/FP/FN/TN)

fisherNormal<-fisher.test(ContMatrix_Normal) #Fisher's Exact test
fisherNormal$p.value #2.114852e-18<- P-value

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##################### RESULTS WITH THE CORRECTION FOR SOME REACTIONS (ERs & SLs) ########
#~~~~~~~~~~~~~~~~~~ Essential reactions
ContMatrix_Normal_Corrected<- as.data.frame(matrix(c(94,175,161,2153),nrow=2)) #make matrix
perc_normal_Corrected<-((94+2153)/2583) *100 #86.99187


fisherNormal_Corrected<-fisher.test(ContMatrix_Normal_Corrected) #Fisher's Exact test
fisherNormal_Corrected$p.value #5.529644e-34 <- P-value

TPNormlist<- LethalRxns_SL %>% filter(LethalRxns_SL$RXNS %in% Gene_correctedLethalRxns$Rxns)  #94 <-true positive
TP<-94
FPNorm<- LethalRxns_SL %>% filter(LethalRxns_SL$RXNS %!in% Gene_correctedLethalRxns$Rxns)  #175 <- false positive
TNNormlist<- AllRxns %>% filter(AllRxns$X1 %!in% LethalRxns_SL$RXNS) 
TNNormlist <- TNNormlist %>% filter(TNNormlist$X1 %!in%  Gene_correctedLethalRxns$Rxns) #2147 <- true negative
TN<- 2153
FN<- 2583-94-2153-175 #161<- False negative
percNorm<- ((2153+94)/2583)*100 #86.99187 <- Prediction ratio
##
## Statistics ## 
#Order or the contingency table / matrix -> (TP/FP/FN/TN)
ContMatrix_Normal<- as.data.frame(matrix(c(94,161,175,2153),nrow=2)) #make matrix

fisherNormal<-fisher.test(ContMatrix_Normal) #Fisher's Exact test
fisherNormal$p.value #5.529644e-34 <- P-value

#~~~~~~~~~~~~~~~~ Synthetic lethals
colnames(SL_normal)<-c('Var1','Var2')
colnames(Gene_CorrectedSynthetic_lethal)<-c('X1','X2')
#
##Values for contingency matrix 
TPSLlist_O<- SL_normal %>% filter((SL_normal$Var1 %in% Gene_CorrectedSynthetic_lethal$X1 & SL_normal$Var2 %in% Gene_CorrectedSynthetic_lethal$X2) |
                                    (SL_normal$Var1 %in% Gene_correctedLethalRxns$Rxns & SL_normal$Var2 %in% Gene_correctedLethalRxns$Rxns)) 
TP<- 28 

FPSLlist_O<- SL_normal %>% filter((SL_normal$Var1 %!in% Gene_CorrectedSynthetic_lethal$X1 | SL_normal$Var2 %!in% Gene_CorrectedSynthetic_lethal$X2) &
                                    (SL_normal$Var1 %!in% Gene_correctedLethalRxns$Rxns | SL_normal$Var2 %!in% Gene_correctedLethalRxns$Rxns)) 
FP<- 233 

FNSLlist_O2<-Gene_CorrectedSynthetic_lethal %>% 
  filter(Gene_CorrectedSynthetic_lethal$X1 %!in% SL_normal$Var1 &
           Gene_CorrectedSynthetic_lethal$X1 %!in% SL_normal$Var2 &
           Gene_CorrectedSynthetic_lethal$X2 %!in% SL_normal$Var2 &
           Gene_CorrectedSynthetic_lethal$X2 %!in% SL_normal$Var1 &
           Gene_CorrectedSynthetic_lethal$X2 %!in% LethalRxns_Exp$RXNS & 
           Gene_CorrectedSynthetic_lethal$X1 %!in% LethalRxns_Exp$RXNS &
           Gene_CorrectedSynthetic_lethal$X2 %!in% LethalRxns_O$X1 & 
           Gene_CorrectedSynthetic_lethal$X1 %!in% LethalRxns_O$X1
           )
FN<-33

TN<-17157-28-236-136 #16894
ContMatrix_Normal<- as.data.frame(matrix(c(TP,FN,FP,TN),nrow=2))
ContMatrix_Normal
percNorm<- ((TN+TP)/possibleCombinations)*100 #97.83179 <- Prediction ratio
percNorm
##
## Statistics ## 

fisherNormal<-fisher.test(ContMatrix_Normal) #Fisher's Exact test
fisherNormal$p.value # 1.912942e-21 <- P-value
