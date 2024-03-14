#IDENTIFY ESSENTIAL REACTIONS FROM LETHAL GENES || June 7 2022 || EXPERIMENTAL DATA SCRIPT 2/3
#EssentialRxnsFromEssentialGenes_Exp
#Description: Get 'essential reactions' that come from the MATLAB list of rxns associated with the essential genes from the Keio Collection 
# NOTE: The final list may also contain genes that will be SL -> as the gene sometimes codes for more than one reaction in the pathway. 
#~~~~~~~~~~~~~~~~~~~~~~~~~~
# Libraries ####
library(readr);library(tidyverse)
'%!in%' <- Negate('%in%')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load files ####
SynthLethals_Exp <- read_delim("Desktop/SynthLethals_Exp.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE) #single lethal genes
SingleDelitions<-read_csv("Desktop/SingleLethalsGenes_Exp_7Jun.csv") #Single genes
SingleDelitionRxns_Exp <- read_csv("Desktop/PotentialLethalRxns_7Jun.csv") #potential lethal rxns

# 0-> Make list of all essential genes #####
AllLethalGenes1<-as.data.frame(unique(SynthLethals_Exp$`Query- strain-bnumberb`))
colnames(AllLethalGenes1)<- 'Gene'
AllLethalGenes2<-as.data.frame(unique(SynthLethals_Exp$`bnumber-of-recipient-mutant-strain`))
colnames(AllLethalGenes2)<- 'Gene'
AllLethalGenesF<- rbind(AllLethalGenes1,AllLethalGenes2)
AllLethalGenesF<- as.data.frame(unique(AllLethalGenesF$Gene))

SingleDelitions<-as.data.frame(SingleDelitions)
colnames(SingleDelitions) <-'Gene'
colnames(AllLethalGenesF) <-'Gene'
AllLethalGenesFINAL<- rbind(SingleDelitions,AllLethalGenesF)
AllLethalGenesFINAL<-unique(AllLethalGenesFINAL)
colnames(AllLethalGenesFINAL) <-'Gene' 
AllLethalGenesFINAL<-as.character(AllLethalGenesFINAL$Gene) #~ Final list of lethal genes ~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1-> Direct Essential Rxns #####
#These rxns are coded by only one gene in the essential genes list

EssentialRxns<-SingleDelitionRxns_Exp %>% filter(nchar(SingleDelitionRxns_Exp$cell5)==5) #~LethalRxns1 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2 -> Reactions in simple logical relations #####
#### 2.1 -> OR ###
#~~~~~ 
# G1 OR G2 -> RXN1 || If G1 = G2 = E -> RXN1 =E ---  G1 != G2 = E -> RXN1 !=E 
#~~~~ 
#Remove () - filter only relations with OR - remove whitespaces
RxnsInOR<- SingleDelitionRxns_Exp %>% filter(str_detect(SingleDelitionRxns_Exp$cell5, " or ")==TRUE) 
  RxnsInONLYOR<- RxnsInOR %>% filter(str_detect(RxnsInOR$cell5, " and ")==FALSE) #Filter out ANDs
  RxnsInOR2<-as.data.frame(str_replace_all(RxnsInONLYOR$cell5, "[(]", ""))
  RxnsInOR2<-rename(RxnsInOR2,Rules=`str_replace_all(RxnsInONLYOR$cell5, "[(]", "")`)
  RxnsInOR2<-as.data.frame(str_replace_all(RxnsInOR2$Rules, "[)]", ""))
  RxnsInOR2<-as.data.frame(str_split(RxnsInOR2$`str_replace_all(RxnsInOR2$Rules, "[)]", "")`," or ",simplify = TRUE))
  RxnsInOR2<-   RxnsInOR2 %>%  mutate(across(where(is.character), str_trim))
  RxnsInOR<-cbind(RxnsInONLYOR[,1],RxnsInOR2[,])
 ## Decide which rxns are Essential from GPRs (Done by #Genes in OR relation)
   ##Two genes 
    dos<-  RxnsInOR %>% filter(RxnsInOR$V3=="")
    Set<- dos
    NewLethals <-Set %>% filter(Set$V1 %in% AllLethalGenesFINAL & Set$V2 %in% AllLethalGenesFINAL)
    Removed <-Set %>% filter( Set$V2 %!in% AllLethalGenesFINAL)
   ## Three genes - no rxns matched all the criteria
    tres<-  RxnsInOR %>% filter(RxnsInOR$V4=="") 
    tres<- tres %>% filter(tres$cell1 %!in%  dos$cell1)
    set<-tres
    NewLethals3 <-Set %>% filter(Set$V1 %in% AllLethalGenesFINAL & Set$V2 %in% AllLethalGenesFINAL & Set$V3 %in% AllLethalGenesFINAL)
  ## Four genes - no rxns matched all the criteria
    four<-  RxnsInOR %>% filter(RxnsInOR$V4!="")
    Set<-four
    NewLethals4 <-Set %>% filter(Set$V1 %in% AllLethalGenesFINAL & Set$V2 %in% AllLethalGenesFINAL & Set$V3 %in% SingleDelitions & Set$V4 %in% SingleDelitions)
   
##Combine new lethal rxns 
NewLethalsOR<- NewLethals  #~LethalRxns2

##### 2.2 -> AND ###
#~~~~~ 
# G1 AND G2 -> RXN1 || If G1 OR G2 = E -> RXN1 =E 
#~~~~
#Remove () - filter only relations with OR - remove whitespaces
RxnsInAND<- SingleDelitionRxns_Exp %>% filter(str_detect(SingleDelitionRxns_Exp$cell5, " and ")==TRUE)  
    RxnsInONLYAND<-RxnsInAND %>% filter(str_detect(RxnsInAND$cell5, " or ")==FALSE) #Filter out OR
    RxnsInONLYAND2<-as.data.frame(str_split(RxnsInONLYAND$cell5," and ",simplify = TRUE))
    RxnsInONLYAND<- cbind(RxnsInONLYAND[,1],RxnsInONLYAND2)
  #Get essential rxns from GPRS (Done by #of genes in AND relation)
    #Two genes
    dos<-  RxnsInONLYAND %>% filter(RxnsInONLYAND$V3=="")
    Set<- dos
    NewLethals <-Set %>% filter(Set$V1 %in% AllLethalGenesFINAL | Set$V2 %in% AllLethalGenesFINAL)
    #No Rxn has 3 genes
    #Four genes 
    cuatro<-  RxnsInONLYAND %>% filter(RxnsInONLYAND$V5=="") 
    cuatro<- cuatro %>% filter(cuatro$cell1 %!in%  dos$cell1)
    Set<-cuatro
    NewLethals3 <-Set %>% filter(Set$V1 %in% AllLethalGenesFINAL | Set$V2 %in% AllLethalGenesFINAL | Set$V3 %in% AllLethalGenesFINAL)
    #Five genes
    cinco<-  RxnsInONLYAND %>% filter(RxnsInONLYAND$V5 != "")  
    Set<-cinco
    NewLethals5 <-Set %>% filter(Set$V1 %in% AllLethalGenesFINAL | Set$V2 %in% AllLethalGenesFINAL | Set$V3 %in% AllLethalGenesFINAL | Set$V5 %in% Set$V1 %in% AllLethalGenesFINAL)
    
 #Join results
  NewLethals5<-as.data.frame(NewLethals5[,1]);colnames(NewLethals5)<-'Rxn'
  NewLethals3<-as.data.frame(NewLethals5[,1]);colnames(NewLethals3)<-'Rxn'
  NewLethals<-as.data.frame(NewLethals5[,1]);colnames(NewLethals)<-'Rxn'
  
  NewLethalsAND<-rbind(NewLethals,NewLethals3,NewLethals5)  #~LethalRxns3
   
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3 -> Reactions in complicated logical relations (ORs and ANDs together) ##### 
#Filter ANDs&ORs - remove () - remove whitespaces
RxnsInOR<- SingleDelitionRxns_Exp %>% filter(str_detect(SingleDelitionRxns_Exp$cell5, " or ")==TRUE) #get all ORs
RxnsInSpecialOR<- RxnsInOR %>% filter(str_detect(RxnsInOR$cell5, " and ")!=FALSE) #Filter in only ORs & ANDs
RxnsInSpecialOR<-RxnsInSpecialOR[,-(2:4)]  

UniqueGPRS<-as.data.frame(unique(RxnsInSpecialOR$cell5))
UniqueGPRS<-as.data.frame(str_replace_all(UniqueGPRS$`unique(RxnsInSpecialOR$cell5)`, "[(]", ""))
UniqueGPRS<-as.data.frame(str_replace_all(UniqueGPRS[,1], "[)]", ""))
colnames(UniqueGPRS)<- c('Rules') 

### Get unique GPRs(7)
UniqueGPRS<-as.data.frame(str_split(UniqueGPRS$Rules," or ",simplify = TRUE))
#~~~~~
#Split the GPR is sets 
# So from (G1 AND G2) OR (G3 AND G4) we get GenSet1 -> (G1,G2) and GenSet2 -> (G3,G4)
#If at least 1 gene from every geneset is E, the Rxn is E
#~~~~~~~

##GPR 1
geneset1=as.data.frame(str_split(UniqueGPRS[1,1]," and ",simplify = TRUE))
geneset2=as.data.frame(str_split(UniqueGPRS[1,2]," and ",simplify = TRUE))
geneset3=as.data.frame(str_split(UniqueGPRS[1,3]," and ",simplify = TRUE))
geneset4=as.data.frame(str_split(UniqueGPRS[1,4]," and ",simplify = TRUE))


  a=geneset1 %in% AllLethalGenesFINAL
  a=sum(a,na.rm = TRUE)
  flagA=1
  
  b=geneset2 %in% AllLethalGenesFINAL 
  b=sum(b,na.rm = TRUE)
  flagB=0 #:(
### GPR 2
geneset1=as.data.frame(str_split(UniqueGPRS[2,1]," and ",simplify = TRUE))
geneset2=as.data.frame(str_split(UniqueGPRS[2,2]," and ",simplify = TRUE))
geneset3=as.data.frame(str_split(UniqueGPRS[2,3]," and ",simplify = TRUE))
geneset4=as.data.frame(str_split(UniqueGPRS[2,4]," and ",simplifyUniqueGPRS))

  a=geneset1 %in% AllLethalGenesFINAL
  a=sum(a,na.rm = TRUE)
  a
  flagA=1
  
  b=geneset2 %in% AllLethalGenesFINAL 
  b=sum(b,na.rm = TRUE)
  b
  flagB=1 #:)

### GPR 3 
geneset1=as.data.frame(str_split(UniqueGPRS[3,1]," and ",simplify = TRUE))
geneset2=as.data.frame(str_split(UniqueGPRS[3,2]," and ",simplify = TRUE))
geneset3=as.data.frame(str_split(UniqueGPRS[3,3]," and ",simplify = TRUE))
geneset4=as.data.frame(str_split(UniqueGPRS[3,4]," and ",simplify = TRUE))

  a=geneset1 %in% AllLethalGenesFINAL
  a=sum(a,na.rm = TRUE)
  a
  flagA=1
                                   
  b=geneset2 %in% AllLethalGenesFINAL 
  b=sum(b,na.rm = TRUE)
  b 
  flagB=0# :( 

### GPR 4
geneset1=as.data.frame(str_split(UniqueGPRS[4,1]," and ",simplify = TRUE))
geneset2=as.data.frame(str_split(UniqueGPRS[4,2]," and ",simplify = TRUE))
geneset3=as.data.frame(str_split(UniqueGPRS[4,3]," and ",simplify = TRUE))
geneset4=as.data.frame(str_split(UniqueGPRS[4,4]," and ",simplify = TRUE))
                                
  a=geneset1 %in% AllLethalGenesFINAL
  a=sum(a,na.rm = TRUE)
  a
  flagA=1
                                 
  b=geneset2 %in% AllLethalGenesFINAL 
  b=sum(b,na.rm = TRUE)
  b # :( 
### GPR 5
geneset1=as.data.frame(str_split(UniqueGPRS[5,1]," and ",simplify = TRUE))
geneset2=as.data.frame(str_split(UniqueGPRS[5,2]," and ",simplify = TRUE))
geneset3=as.data.frame(str_split(UniqueGPRS[5,3]," and ",simplify = TRUE))
geneset4=as.data.frame(str_split(UniqueGPRS[5,4]," and ",simplify = TRUE))
                                 
  a=geneset1 %in% AllLethalGenesFINAL
  a=sum(a,na.rm = TRUE)
  a
  flagA=1
                                 
  b=geneset2 %in% AllLethalGenesFINAL 
  b=sum(b,na.rm = TRUE)
  b 
  flagB=1

  c=geneset3 %in% AllLethalGenesFINAL 
  c=sum(b,na.rm = TRUE)
  c 
  flagC=1 #:)
### GPR 6
geneset1=as.data.frame(str_split(UniqueGPRS[6,1]," and ",simplify = TRUE))
geneset2=as.data.frame(str_split(UniqueGPRS[6,2]," and ",simplify = TRUE))
geneset3=as.data.frame(str_split(UniqueGPRS[6,3]," and ",simplify = TRUE))
geneset4=as.data.frame(str_split(UniqueGPRS[6,4]," and ",simplify = TRUE))
                                 
  a=geneset1 %in% AllLethalGenesFINAL
  a=sum(a,na.rm = TRUE)
  a
  flagA=1
                                 
  b=geneset2 %in% AllLethalGenesFINAL 
  b=sum(b,na.rm = TRUE)
  b 
  flagB=2 #:)
### GPR 7
geneset1=as.data.frame(str_split(UniqueGPRS[7,1]," and ",simplify = TRUE))
geneset2=as.data.frame(str_split(UniqueGPRS[7,2]," and ",simplify = TRUE))
geneset3=as.data.frame(str_split(UniqueGPRS[7,3]," and ",simplify = TRUE))
geneset4=as.data.frame(str_split(UniqueGPRS[7,4]," and ",simplify = TRUE))
                                 
  a=geneset1 %in% AllLethalGenesFINAL
  a=sum(a,na.rm = TRUE)
  a
  flagA=1
                                 
 b=geneset2 %in% AllLethalGenesFINAL 
 b=sum(b,na.rm = TRUE)
 b 
 flagB=2 #:)                                 
####
#The GPRS 7,6,5 & 2 belong to essential rxns 
 NewLethalsSp<-RxnsANDOR[-c(9:14),1] #~~ LethalRxnsList4
                                 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                                 
# 4 ->Combine all rxns  #####
EssentialRxns<-EssentialRxns[,1]
colnames(EssentialRxns)<-"Rxns"
NewLethalsOR<-NewLethalsOR[,1]
colnames(NewLethalsOR)<-"Rxns"
colnames(NewLethalsAND)<-"Rxns"
colnames(NewLethalsSp)<-"Rxns"
FinalEssentialRxnsList<-rbind(EssentialRxns,NewLethalsOR,NewLethalsAND,NewLethalsSp)
FinalEssentialRxnsList<-unique(FinalEssentialRxnsList) #~~~~~~ FinalList
#save
write_csv(FinalEssentialRxnsList,'LethalRxns_EXP_7JUN.csv')
