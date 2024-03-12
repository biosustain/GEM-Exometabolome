### RESULTS ANALYSIS - SL METHOD VALIDATION -- PLOTS|| JUNE 2022
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~ ORIGINAL DATA 
names<-c("1. ParallelSL+Original","2. ParallelSL+Loopless","3. ParallelSL+Compacted+Loopless","4.ParallelSL+ROOM+Compacted+Loopless",
         "5. ROOM+Compacted+Loopless",
         "1. ParallelSL+Original","2. ParallelSL+Loopless","3. ParallelSL+Compacted+Loopless","4.ParallelSL+ROOM+Compacted+Loopless",
         "5. ROOM+Compacted+Loopless")

TP_plot<-c(94,95,88,47,49,25,28,27,12,0.01)
FP_plot<-c(175,179,170,113,106,236,247,256,48,0)
Rxn<-c('Essential reactions','Essential reactions','Essential reactions','Essential reactions','Essential reactions',
       'Synthetic lethals','Synthetic lethals','Synthetic lethals','Synthetic lethals','Synthetic lethals')

data<-as.data.frame(cbind(names,Rxn,TP_plot,FP_plot))
data$TP_plot<- as.numeric(data$TP_plot)
data$FP_plot<- as.numeric(data$FP_plot)
data<-data%>% pivot_longer(cols = c(3:4),names_to = "Type") 

colors<- c("#B8C9E4","#95AAD3")

ggplot(data)+geom_bar(aes(x=names,y=value,fill=Type),position="stack", stat="identity")+ylim(0,300)+
  ylab("Predictions")+xlab(NULL)+scale_fill_manual(values=colors,labels=c("False Positive","True Positive"))+
  facet_grid(cols=vars(Rxn))+theme_minimal()+theme( panel.spacing = unit(2, "lines"),
                                                    axis.text.x = element_text(size=10,angle=20,hjust =0.7),
                                                    strip.background = element_rect(color="black", fill="gray88", size=1, linetype="solid"),
                                                    strip.text.x = ggtext::element_markdown(face='bold',size =10),plot.margin = margin(1,1,1,1,"cm"))

#~~~~~~~~~~~~ GENE CORRECTED ~~~~~~~~~~~~~~~

TP_plot<-c(94,94,87,46,43,28,28,27,13,0)
FP_plot<-c(175,180,171,114,107,233,247,256,48,0)

data<-as.data.frame(cbind(names,Rxn,TP_plot,FP_plot))

data$TP_plot<- as.numeric(data$TP_plot)
data$FP_plot<- as.numeric(data$FP_plot)
data<-data%>% pivot_longer(cols = c(3:4),names_to = "Type") 

colors<- c("#B8C9E4","#95AAD3")

ggplot(data)+geom_bar(aes(x=names,y=value,fill=Type),position="stack", stat="identity")+ylim(0,300)+
  ylab("Predictions")+xlab(NULL)+scale_fill_manual(values=colors,labels=c("False Positive","True Positive"))+
  facet_grid(cols=vars(Rxn))+theme_minimal()+
  theme( strip.text.x = ggtext::element_markdown(, face='bold',size =10),
         strip.background = element_rect(color="black", fill="gray88", size=1, linetype="solid"),
         legend.position = "bottom",legend.title = ggtext::element_markdown(family = 'Times', face='bold'),
         axis.text.x = ggtext::element_markdown(face='bold',,size =8,,angle=20),
         axis.title.y = ggtext::element_markdown(face='bold'),plot.margin = margin(1,1,1,1,"cm"),
         panel.spacing = unit(2, "lines"))
 


