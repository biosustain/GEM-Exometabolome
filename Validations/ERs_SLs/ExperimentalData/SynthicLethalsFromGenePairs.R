#IDENTIFY SYNTHETIC LETHALS FROM SYNTHETIC LETHAL GENES || EXPERIMENTAL DATA SCRIPT 3/3
#SynthicLethalsFromGenePairs
#Description: Get synthetic lethals from the MATLAB list of rxns associated with the pair of lethal genes. 
#NOTE:Some reactions that should be in this list ended up in the list of single lethal reactions, as they were associated 
# and thus identified with/from an essential gene. 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Libraries ####
library(readr); library(tidyverse)
# 0 -> Load data ####
SynthLethals_Exp <- read_delim("Desktop/SingleLethalsGenes_Exp_7Jun.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
tableFinalGene1 <- read_csv("Desktop/tableFinalGene1.csv")
tableFinalGene2 <- read_csv("Desktop/tableFinalGene2.csv")
#Filter data to only genes 
tableFinalGene1<- tableFinalGene1[,-c(2:4)]
tableFinalGene2<- tableFinalGene2[,-c(2:4)]
SynthLethals_Exp<- SynthLethals_Exp[,c(2,4)]

# 1 -> get reaction pairs from list 
combo<-empty.dump()
combonew<-empty.dump()
specialCaseGen1<- empty.dump()
specialCaseGen2<- empty.dump()

for (i in c(1:138)) #for every pair in the list 
{
  genA<-as.character(SynthLethals_Exp[i,1]) #Get genes
  genB<-as.character(SynthLethals_Exp[i,2])
  
  
  GPR1Gen1=tableFinalGene1 %>% filter(str_detect(tableFinalGene1$cell5,genA)==TRUE) #Get all rxns for those genes
  GPR1Gen2=tableFinalGene2 %>% filter(str_detect(tableFinalGene2$cell25,genB)==TRUE)
  
  x=dim(GPR1Gen1 )#check != 0 --- this is necessary as some genes are not in the model and here will be 0 and might cause an error afterwards
  y=dim(GPR1Gen2)
  
  if((x[1]!=0) & (y[1]!=0)){
    RxnsGen1<-unique(GPR1Gen1[,1]) #Get unique rxns
    RxnsGen2<-unique(GPR1Gen2[,1])
    
    if (str_detect(GPR1Gen2$cell25," or ")==FALSE & str_detect(GPR1Gen1$cell5," or ")==FALSE){ # is no ORs are detected
      combonew<-expand.grid(RxnsGen2$cell21,RxnsGen1$cell1) #combine all rxns in the same pair of genes
      combo<-rbind(combo,combonew)
    } else{   
      specialCaseGen1New= RxnsGen1 #If ORs are detected, they go here so they can be looked in detail 
      specialCaseGen2New= RxnsGen2 
      
      specialCaseGen1=rbind(specialCaseGen1New,specialCaseGen1)
      specialCaseGen2=rbind(specialCaseGen2New,specialCaseGen2)
      }
  }
}

##NOTE: None of the reations with ORs were included in the final list, 
# they had to meet the G1 OR G2 -> RXN1 || If G1 = G2 = E -> RXN1 =E ---  G1 != G2 = E -> RXN1 !=E  rule
# and didnt.  
combo #~~ Final List of synthetic lethals. 
write_csv(combo,'Synthetic lethal_Rxns_Exp.csv')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
#Test code ####
genA<-as.character(SynthLethals_Exp[3,1])
genB<-as.character(SynthLethals_Exp[3,2])


GPR1Gen1=tableFinalGene1 %>% filter(str_detect(tableFinalGene1$cell5,genA)==TRUE)
GPR1Gen2=tableFinalGene2 %>% filter(str_detect(tableFinalGene2$cell25,genB)==TRUE)


PR1Gen1=tableFinalGene1 %>% filter(str_detect(tableFinalGene1$cell5,"b1680")==TRUE)
GPR1Gen2=tableFinalGene2 %>% filter(str_detect(tableFinalGene2$cell25,"b2530")==TRUE)

RxnsGen2<-unique(GPR1Gen2[,1])
RxnsGen2<-as.character(RxnsGen2$cell21)
RxnsGen1<-unique(GPR1Gen1[,1])

B1680_B2530<-expand.grid(RxnsGen2$cell21,RxnsGen1$cell1)

