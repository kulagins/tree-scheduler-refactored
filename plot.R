setwd("~/ChangjiangGou/MemCom/journal/code/result")
library(ggplot2)
library(reshape2)
library(tikzDevice)
rm(list=ls())

#figures in the journal paper
#-----------------makespan of stage 1, figure 6----------------------
makespan.long<-read.table("./stage1/result_phase1.txt",header=TRUE)

long.copy<-makespan.long
makespan.wide<-dcast(makespan.long,TreeName+NPR+CCR+AmountProcessors~Stage1,value.var = "Makespan")
makespan.wide$ASAP<-makespan.wide$ASAP/makespan.wide$Sequence
makespan.wide$AvoidChain<-makespan.wide$AvoidChain/makespan.wide$Sequence
makespan.wide$ImprovedSplit<-makespan.wide$ImprovedSplit/makespan.wide$Sequence
makespan.wide$SplitSubtrees<-makespan.wide$SplitSubtrees/makespan.wide$Sequence

temp.wide<-makespan.wide[makespan.wide$NPR==10000,]
temp.wide<-temp.wide[temp.wide$CCR%in%c(0.1,1,10),]
sum(temp.wide$ImprovedSplit<=1)/nrow(temp.wide)
#on 56% cases, ImprovedSplit is better than or equal to SplitSubtrees

long<-melt(makespan.wide,id.vars = c("TreeName","NPR","CCR","AmountProcessors"),variable.name = "Heuristics",value.name = "makespan" )
long<-long[long$NPR%in%c(100,1000,10000),]
long<-long[long$CCR%in%c(0.1,1,10),]
long$NPR<-1/long$NPR
long$NPR<-as.factor(long$NPR)
long$CCR<-as.factor(long$CCR)

levels(long$NPR)[levels(long$NPR)=="1e-04"]<-"PNR = 1e-04"
levels(long$NPR)[levels(long$NPR)=="0.001"]<-"PNR = 0.001"
levels(long$NPR)[levels(long$NPR)=="0.01"]<-"PNR = 0.01"

long<-long[which(long$Heuristics!="Sequence"),]

aggregate(long$makespan,list(long$NPR,long$Heuristics),quantile,c(0.5))

hues = seq(15, 375, length = length(levels(long$Heuristics)) + 1)
hcl(h = hues, l = 65, c = 100)[1:length(levels(long$Heuristics))]
levels(long$Heuristics)
cb_palette <- c(ASAP="#F8766D", AvoidChain="#B79F00",
                ImprovedSplit="#00BA38", SplitSubtrees="#619CFF")

#tikz("~/ChangjiangGou/MemCom/journal/figure_makespan1.tex",width = 7.2,height = 3)
ggplot(long,aes(x=CCR,y=makespan,fill=Heuristics))+geom_boxplot(outlier.size = 0.05)+
  scale_x_discrete(name="CCR")+scale_y_continuous(name="Makespan normalized to \\textbf{Sequence}")+
  guides(fill=FALSE)+
  scale_fill_manual(labels=c("ASAP","ASAPnochain","ImprovedSplit","SplitSubtrees"),values = cb_palette)+
  labs(fill="")+theme(legend.position = c(0.799,0.85),legend.background = element_blank(),legend.key = element_blank())+
  facet_grid(.~NPR)#+
  #geom_hline(yintercept=0.55)
dev.off()

long<-long.copy
long<-long[long$CCR%in%c(0.1,1,10),]
long<-long[long$NPR%in%c(100,1000,10000),]
long$NPR<-1/long$NPR
long$NPR<-as.factor(long$NPR)
long$CCR<-as.factor(long$CCR)
long$AmountProcessors[long$AmountProcessors<=2]=3
long$perc<-long$AmountSubtrees/long$AmountProcessors
long<-long[long$Stage1!="Sequence",]

levels(long$NPR)[levels(long$NPR)=="1e-04"]<-"PNR = 1e-04"
levels(long$NPR)[levels(long$NPR)=="0.001"]<-"PNR = 0.001"
levels(long$NPR)[levels(long$NPR)=="0.01"]<-"PNR = 0.01"

temp<-aggregate(long$perc,list(long$NPR,long$Stage1),quantile,c(0.5))
View(temp)

temp<-aggregate(long$perc,list(long$NPR,long$Stage1),mean)
View(temp)

wide<-dcast(long,TreeName+NPR+CCR~Stage1,value.var = "perc")
wide$AvoidChain<-wide$AvoidChain-wide$ASAP
mean(wide$AvoidChain)

levels(long$Stage1)[levels(long$Stage1)=="AvoidChain"]<-"ASAPnochain"

cb_palette <- c(ASAP="#F8766D", ASAPnochain="#B79F00",
                ImprovedSplit="#00BA38", SplitSubtrees="#619CFF")#,
                #FirstFit="#FFFFFF")

#-------------------------amount of subtrees after stage1, figure 7
#tikz("~/ChangjiangGou/MemCom/journal/figure_amountSubtrees_stage1.tex",width = 7.2,height = 3)
ggplot(long,aes(x=CCR,y=perc,fill=Stage1))+geom_boxplot(outlier.size = 0.05)+
  scale_x_discrete(name="CCR")+
  scale_y_continuous(name="Number of subtrees to processors")+
  labs(fill="")+
  guides(fill=guide_legend(nrow = 1))+
  scale_fill_manual(values = cb_palette)+
  theme(legend.position = "top",legend.background = element_blank(),legend.key = element_blank(),
        axis.text.y = element_text(angle = 40, vjust=0.5))+
  facet_grid(.~NPR)#+
  #geom_hline(yintercept=0.16)
dev.off()

#---------------------how many trees does not fit in memory, the big table---------------------
rm(list=ls())

rawdata<-read.table("./stage23/result_twothreeStage.txt",header = TRUE)
long<-rawdata[is.na(rawdata$Stage3),]#select stage 1 and stage 2
long<-long[long$MemoryConstraint==1|is.na(long$MemoryConstraint),]#select memory constraint as the maxoutd or NA

#------------ branch here 
CCRoption=10
long<-long[!is.na(long$Stage2),]#select data after stage2
long<-long[,c("TreeName","NPR","CCR","AmountSubtrees","AmountProcessors","Makespan","Stage1","Stage2")]
long<-unique(long)
long<-long[long$CCR==CCRoption,]
long$AmountProcessors[which(long$AmountProcessors<3)]=3
long$ratio<-long$AmountSubtrees/long$AmountProcessors

ratio.mean<-aggregate(long$ratio,list(long$NPR,long$Stage1,long$Stage2),mean)
colnames(ratio.mean)<-c("NPR","Stage1","Stage2","ratio")
ratio.mean$NPR<-1/ratio.mean$NPR
colnames(ratio.mean)[1]<-"PNR"
ratio.mean$ratio<-round(ratio.mean$ratio,digits=2)
View(ratio.mean)#the left column of all

#------------ branch here
long<-rawdata[is.na(rawdata$Stage3),]#select stage 1 and stage 2
long<-long[long$MemoryConstraint==1|is.na(long$MemoryConstraint),]#select memory constraint as the maxoutd or NA
long$Heuristic<-paste(long$Stage1,"+",long$Stage2,sep="")
long<-unique(long)
colnames(long)
wide<-dcast(long,TreeName+NPR+CCR~Heuristic,value.var = "Makespan")

#decrease in makespan
wide$`AvoidChain+FirstFit`<-wide$`AvoidChain+FirstFit`/wide$`AvoidChain+NA`
wide$`AvoidChain+Immediately`<-wide$`AvoidChain+Immediately`/wide$`AvoidChain+NA`
wide$`AvoidChain+LargestFirst`<-wide$`AvoidChain+LargestFirst`/wide$`AvoidChain+NA`

wide$`ImprovedSplit+FirstFit`<-wide$`ImprovedSplit+FirstFit`/wide$`ImprovedSplit+NA`
wide$`ImprovedSplit+Immediately`<-wide$`ImprovedSplit+Immediately`/wide$`ImprovedSplit+NA`
wide$`ImprovedSplit+LargestFirst`<-wide$`ImprovedSplit+LargestFirst`/wide$`ImprovedSplit+NA`

wide$`SplitSubtrees+FirstFit`<-wide$`SplitSubtrees+FirstFit`/wide$`SplitSubtrees+NA`
wide$`SplitSubtrees+Immediately`<-wide$`SplitSubtrees+Immediately`/wide$`SplitSubtrees+NA`
wide$`SplitSubtrees+LargestFirst`<-wide$`SplitSubtrees+LargestFirst`/wide$`SplitSubtrees+NA`

long<-melt(wide,id.vars = c("TreeName","NPR","CCR"),variable.name = "Heuristic",value.name = "ratio")
long<-long[!long$Heuristic%in%c("AvoidChain+NA","ImprovedSplit+NA","SplitSubtrees+NA"),]
long.sub<-long[long$CCR==CCRoption,]

heuristic<-strsplit(as.character(long.sub$Heuristic),"\\+")
stage1<-unlist(lapply(seq(1:nrow(long.sub)),function(x){heuristic[[x]][[1]]}))
stage2<-unlist(lapply(seq(1:nrow(long.sub)),function(x){heuristic[[x]][[length(heuristic[[x]])]]}))
long.sub$stage1<-stage1
long.sub$stage2<-stage2

decrease.mean<-aggregate(long.sub$ratio,list(long.sub$stage1,long.sub$stage2,long.sub$NPR),mean)
decrease.mean$x<-1-decrease.mean$x
decrease.mean$Group.3<-1/decrease.mean$Group.3
colnames(decrease.mean)<-c("Stage1","Stage2","PNR","decrease")
decrease.mean$decrease<-round(decrease.mean$decrease,digits = 2)
View(decrease.mean)#the right column of all
#done

rm(list=ls())

df<-read.table("./sequence/result_sequence.txt",header = TRUE)
colnames(df)[8]<-"Heuristic"
memoryOption<-1 #1 is the maxoutd, the strict case; 3 is the minmem, the loose case
CCRoption=10
long<-df[df$MemoryConstraint==memoryOption|is.na(df$MemoryConstraint),]

colnames(long)
wide.subtrees<-dcast(long,TreeName+NPR+CCR+AmountProcessors~Heuristic,value.var = "AmountSubtrees")
wide.makespan<-dcast(long,TreeName+NPR+CCR~Heuristic,value.var = "Makespan")

wide.subtrees$FIRST_FIT<-wide.subtrees$FIRST_FIT/wide.subtrees$AmountProcessors
wide.subtrees$LARGEST_FIT<-wide.subtrees$LARGEST_FIT/wide.subtrees$AmountProcessors
wide.subtrees$IMMEDIATELY<-wide.subtrees$IMMEDIATELY/wide.subtrees$AmountProcessors

long<-melt(wide.subtrees,id.vars = c("TreeName","NPR","CCR","AmountProcessors"),variable.name = "Heuristic",value.name = "Ratio")
long<-long[long$Heuristic%in%c("FIRST_FIT","LARGEST_FIT","IMMEDIATELY"),]
long<-long[long$CCR==CCRoption,]

decrease.mean<-aggregate(long$Ratio,list(long$Heuristic,long$NPR),mean)
colnames(decrease.mean)<-c("heuristic","PNR","mean")
decrease.mean$PNR<-1/decrease.mean$PNR
decrease.mean$mean<-round(decrease.mean$mean,digits = 2)
View(decrease.mean)#left column of Sequence
#done

wide.makespan$FIRST_FIT<-wide.makespan$FIRST_FIT/wide.makespan$Sequence
wide.makespan$LARGEST_FIT<-wide.makespan$LARGEST_FIT/wide.makespan$Sequence
wide.makespan$IMMEDIATELY<-wide.makespan$IMMEDIATELY/wide.makespan$Sequence

long<-melt(wide.makespan,id.vars = c("TreeName","NPR","CCR"),variable.name = "Heuristic",value.name = "Ratio")
long<-long[long$Heuristic%in%c("FIRST_FIT","LARGEST_FIT","IMMEDIATELY"),]
long<-long[long$CCR==CCRoption,]

decrease.mean<-aggregate(long$Ratio,list(long$Heuristic,long$NPR),mean)
colnames(decrease.mean)<-c("heuristic","PNR","mean")
decrease.mean$PNR<-1/decrease.mean$PNR
decrease.mean$mean<-1-decrease.mean$mean
decrease.mean$mean<-round(decrease.mean$mean,digits = 2)
View(decrease.mean)#right column of Sequence
#done

#---------------------makespan of stage three--------------------------

#-------the performance of Merge or SplitAgain, figure 7/12/14/17-----------
rm(list=ls())
rawdata<-read.table("./stage23/result_twothreeStage.txt",header = TRUE,stringsAsFactors = FALSE)
sequencedata<-read.table("./sequence/result_sequence.txt",header = TRUE, stringsAsFactors = FALSE)
long<-rawdata

memoryOption<-1 #1 is the maxoutd, the strict case; 3 is the minmem, the loose case
stage2Option<-"Immediately" # stage 2 method: FirstFit, Immediately, LargestFirst
CCRoption<-0.1

index_stage2<-which(long$MemoryConstraint==memoryOption&long$Stage2==stage2Option&is.na(long$Stage3))
index_stage3<-which(long$MemoryConstraint==memoryOption&long$Stage2==stage2Option&!is.na(long$Stage3))

index<-c(index_stage2,index_stage3)
long<-long[index,]
long<-long[long$CCR==CCRoption,]
long$Heuristic<-NA

colnames(sequencedata)
sequence.subset<-subset(sequencedata,CCR==CCRoption&MemoryConstraint==memoryOption)

heuristics<-strsplit(as.character(sequence.subset$Heuristic),"\\+")
stage2<-unlist(lapply(seq(1:nrow(sequence.subset)),function(x){heuristics[[x]][[1]]}))
stage3<-unlist(lapply(seq(1:nrow(sequence.subset)),function(x){
  if(length(heuristics[[x]])==2){heuristics[[x]][[2]]}
  else return(NA)}))

sequence.subset$Stage1<-"Sequence"
sequence.subset$Stage2<-"default"
sequence.subset$Stage2<-stage2
sequence.subset$Stage3<-stage3

sequence.subset$Stage3[which(sequence.subset$Stage3=="NA")]="Nothing"
sequence.subset$Stage2[which(sequence.subset$Stage2=="FIRST_FIT")]="FirstFit"
sequence.subset$Stage2[which(sequence.subset$Stage2=="LARGEST_FIT")]="LargestFirst"
sequence.subset$Stage2[which(sequence.subset$Stage2=="IMMEDIATELY")]="Immediately"

sequence.subset<-subset(sequence.subset,Stage2==stage2Option)

long<-rbind(long,sequence.subset)

long$Heuristic<-paste(long$Stage1,long$Stage2,long$Stage3,sep="+")

table(long$Stage3)
table(long$Heuristic)

long$s3method<-!is.na(long$Stage3)
long$s3method<-paste(long$Stage1,long$Stage2,long$s3method,sep="+")

mergeCases<-long$Makespan[which(long$Stage3=="Merge")] 
sum(mergeCases<0)/length(mergeCases)#failure ratio of using merge

wide.ms<-dcast(long,TreeName+NPR+CCR~s3method,value.var = "Makespan")

if(stage2Option=="FirstFit"){
  wide.ms$`AvoidChain+FirstFit+TRUE`<-wide.ms$`AvoidChain+FirstFit+TRUE`/wide.ms$`AvoidChain+FirstFit+FALSE`
  wide.ms$`ImprovedSplit+FirstFit+TRUE`<-wide.ms$`ImprovedSplit+FirstFit+TRUE`/wide.ms$`ImprovedSplit+FirstFit+FALSE`
  wide.ms$`SplitSubtrees+FirstFit+TRUE`<-wide.ms$`SplitSubtrees+FirstFit+TRUE`/wide.ms$`SplitSubtrees+FirstFit+FALSE`
  wide.ms$`Sequence+FirstFit+TRUE`<-wide.ms$`Sequence+FirstFit+TRUE`/wide.ms$`Sequence+FirstFit+FALSE`
}

if(stage2Option=="Immediately"){
  wide.ms$`AvoidChain+Immediately+TRUE`<-wide.ms$`AvoidChain+Immediately+TRUE`/wide.ms$`AvoidChain+Immediately+FALSE`
  wide.ms$`ImprovedSplit+Immediately+TRUE`<-wide.ms$`ImprovedSplit+Immediately+TRUE`/wide.ms$`ImprovedSplit+Immediately+FALSE`
  wide.ms$`SplitSubtrees+Immediately+TRUE`<-wide.ms$`SplitSubtrees+Immediately+TRUE`/wide.ms$`SplitSubtrees+Immediately+FALSE`
  wide.ms$`Sequence+Immediately+TRUE`<-wide.ms$`Sequence+Immediately+TRUE`/wide.ms$`Sequence+Immediately+FALSE`
}

if(stage2Option=="LargestFirst"){
  wide.ms$`AvoidChain+LargestFirst+TRUE`<-wide.ms$`AvoidChain+LargestFirst+TRUE`/wide.ms$`AvoidChain+LargestFirst+FALSE`
  wide.ms$`ImprovedSplit+LargestFirst+TRUE`<-wide.ms$`ImprovedSplit+LargestFirst+TRUE`/wide.ms$`ImprovedSplit+LargestFirst+FALSE`
  wide.ms$`SplitSubtrees+LargestFirst+TRUE`<-wide.ms$`SplitSubtrees+LargestFirst+TRUE`/wide.ms$`SplitSubtrees+LargestFirst+FALSE`
  wide.ms$`Sequence+LargestFirst+TRUE`<-wide.ms$`Sequence+LargestFirst+TRUE`/wide.ms$`Sequence+LargestFirst+FALSE`
}

long.ms<-melt(wide.ms,id.vars = c("TreeName","NPR","CCR"),variable.name = "Heuristic",value.name = "makespan" )

if(stage2Option=="FirstFit"){
  long.ms<-long.ms[long.ms$Heuristic%in%c("AvoidChain+FirstFit+TRUE","ImprovedSplit+FirstFit+TRUE","SplitSubtrees+FirstFit+TRUE","Sequence+FirstFit+TRUE"),]
}

if(stage2Option=="Immediately"){
  long.ms<-long.ms[long.ms$Heuristic%in%c("AvoidChain+Immediately+TRUE","ImprovedSplit+Immediately+TRUE","SplitSubtrees+Immediately+TRUE","Sequence+Immediately+TRUE"),]
}

if(stage2Option=="LargestFirst"){
  long.ms<-long.ms[long.ms$Heuristic%in%c("AvoidChain+LargestFirst+TRUE","ImprovedSplit+LargestFirst+TRUE","SplitSubtrees+LargestFirst+TRUE","Sequence+LargestFirst+TRUE"),]
}

long.ms<-long.ms[order(long.ms$NPR),]
long.ms$makespan[which(long.ms$makespan<0)]<-NA

long.ms$x<-1:nrow(long.ms)
long.ms$y<-1:nrow(long.ms)
for(i in 0:2){#NPR
  for(k in 0:3){#heuristic
    for(j in 0:2){#row
      long.ms$x[124*i+8*j+(1:8)+31*k]<-i*8+c(1:8)
      long.ms$y[124*i+8*j+(1:8)+31*k]<-k*4+j+1
    }
    long.ms$x[124*i+24+(1:7)+31*k]<-i*8+c(1:7)
    long.ms$y[124*i+24+(1:7)+31*k]<-k*4+4
  }
}

ggplot(long.ms,aes(x=x,y=y,fill=makespan))+geom_raster()+
  scale_x_reverse(name="Processor to Node Ratio",breaks=c(4.5,12.5,20.5),labels=c(0.01,0.001,1e-04))+
  scale_y_continuous(breaks=c(2.5,6.5,10.5,14.5),labels = c("ASAPnochain","ImprovedSplit","Sequence","SplitSubtrees"),name="Heuristic")+
  scale_fill_gradient2(low = "green", midpoint=1, mid="white",high = "red",limits=c(0,6),breaks=c(0.2,0.4,0.8,1,2,4,6),guide = "legend")+
  theme_bw()+
  theme(axis.text.y = element_text(angle=70, hjust=0.5, vjust=0),panel.border = element_blank())+
  labs(fill="Performance")+
  annotate("text",x=long.ms$x[is.na(long.ms$makespan)],y=long.ms$y[is.na(long.ms$makespan)],label="F")+
  geom_hline(yintercept = c(4.5,8.5,12.5),linetype="dashed",size=0.2)+
  geom_vline(xintercept = c(8.5,16.5),linetype="dashed",size=0.2)

if(memoryOption==1){m<-"strictM"}
if(memoryOption==3){m<-"looseM"}

output.name<-paste("../../figure_makespanReduced_",m,"_",stage2Option,".png",sep = "")
output.name
ggsave(output.name,width = 200, height=110,units = "mm",dpi=400)

#-----------------makespan of all heuristics, figure 8/ 10/ 15/ 18---------
#merge all data
rm(list=ls())
sequencedata<-read.table("./sequence/result_sequence.txt",header = TRUE, stringsAsFactors = FALSE)

colnames(sequencedata)[8]<-"Heuristic"

memoryOption<-1 #1 is the maxoutd, the strict case; 3 is the minmem, the loose case
stage2Option<-"Immediately" # stage 2 method: FirstFit, Immediately, LargestFirst
sequence.subset<-subset(sequencedata,MemoryConstraint==memoryOption&CCR%in%c(0.1,1,10))

if(stage2Option=="LargestFirst"){
  if(memoryOption==1)
    sequence.subset<-subset(sequence.subset,Heuristic%in%c("FIRST_FIT","LARGEST_FIT+Merge","LARGEST_FIT+NA","LARGEST_FIT+SplitAgain"))
  else
    sequence.subset<-subset(sequence.subset,Heuristic%in%c("LARGEST_FIT+Merge","LARGEST_FIT+NA","LARGEST_FIT+SplitAgain"))
  
  sequence.subset$Heuristic[which(sequence.subset$Heuristic=="LARGEST_FIT+Merge")]<-"LARGEST_FIT+#"
  sequence.subset$Heuristic[which(sequence.subset$Heuristic=="LARGEST_FIT+SplitAgain")]<-"LARGEST_FIT+#"
  sequence.subset$Heuristic[which(sequence.subset$Heuristic=="LARGEST_FIT+NA")]<-"LARGEST_FIT+#"
}

if(stage2Option=="FirstFit"){
  if(memoryOption==1)
    sequence.subset<-subset(sequence.subset,Heuristic%in%c("FIRST_FIT","FIRST_FIT+Merge","FIRST_FIT+NA","FIRST_FIT+SplitAgain"))
  else
    sequence.subset<-subset(sequence.subset,Heuristic%in%c("FIRST_FIT+Merge","FIRST_FIT+NA","FIRST_FIT+SplitAgain"))
  
  sequence.subset$Heuristic[which(sequence.subset$Heuristic=="FIRST_FIT+Merge")]<-"FIRST_FIT+#"
  sequence.subset$Heuristic[which(sequence.subset$Heuristic=="FIRST_FIT+SplitAgain")]<-"FIRST_FIT+#"
  sequence.subset$Heuristic[which(sequence.subset$Heuristic=="FIRST_FIT+NA")]<-"FIRST_FIT+#"
}

if(stage2Option=="Immediately"){
  if(memoryOption==1)
    sequence.subset<-subset(sequence.subset,Heuristic%in%c("FIRST_FIT","IMMEDIATELY+Merge","IMMEDIATELY+NA","IMMEDIATELY+SplitAgain"))
  else
    sequence.subset<-subset(sequence.subset,Heuristic%in%c("IMMEDIATELY+Merge","IMMEDIATELY+NA","IMMEDIATELY+SplitAgain"))
  
  sequence.subset$Heuristic[which(sequence.subset$Heuristic=="IMMEDIATELY+Merge")]<-"IMMEDIATELY+#"
  sequence.subset$Heuristic[which(sequence.subset$Heuristic=="IMMEDIATELY+SplitAgain")]<-"IMMEDIATELY+#"
  sequence.subset$Heuristic[which(sequence.subset$Heuristic=="IMMEDIATELY+NA")]<-"IMMEDIATELY+#"
}

#temp<-df[df$Heuristic=="FirstFit+#",]
#sum(temp$Makespan<0)/nrow(temp) #failure rate of FirstFit+#

makespan.others<-read.table("./stage23/result_twothreeStage.txt",header = TRUE)

if(memoryOption==3){
  splitsubtrees<-subset(makespan.others,Stage1=="SplitSubtrees"&is.na(Stage2))
  splitsubtrees<-unique(splitsubtrees)
  splitsubtrees$Stage1<-"s1"
}

makespan.others<-makespan.others[makespan.others$CCR%in%c(0.1,1,10),]
index<-which(makespan.others$MemoryConstraint==memoryOption&!is.na(makespan.others$Stage3))
makespan.others<-makespan.others[index,]
makespan.others<-makespan.others[makespan.others$Stage2==stage2Option,]

makespan.others$Heuristic<-makespan.others$Stage1;

makespan.others<-makespan.others[,c("TreeName","NPR","CCR","MemoryConstraint","AmountProcessors","Makespan","Heuristic","TimeConsuming")]
sequence.subset<-sequence.subset[,c("TreeName","NPR","CCR","MemoryConstraint","AmountProcessors","Makespan","Heuristic","TimeConsuming")]

makespan.all<-rbind(sequence.subset,makespan.others)

if(memoryOption==3){
  colnames(splitsubtrees)[8]<-"Heuristic"
  splitsubtrees<-splitsubtrees[,c("TreeName","NPR","CCR","MemoryConstraint","AmountProcessors","Makespan","Heuristic","TimeConsuming")]
  makespan.all<-rbind(makespan.all,splitsubtrees)
}

wide.makespan<-dcast(makespan.all,TreeName+NPR+CCR~Heuristic,value.var = "Makespan")

if(memoryOption==1){
  if(stage2Option=="LargestFirst")
    wide.makespan$Select<-apply(wide.makespan[,c('AvoidChain','FIRST_FIT','ImprovedSplit','LARGEST_FIT+#','SplitSubtrees')],1,min)
  if(stage2Option=="FirstFit")
    wide.makespan$Select<-apply(wide.makespan[,c('AvoidChain','FIRST_FIT','ImprovedSplit','FIRST_FIT+#','SplitSubtrees')],1,min)
  if(stage2Option=="Immediately")
    wide.makespan$Select<-apply(wide.makespan[,c('AvoidChain','FIRST_FIT','ImprovedSplit','IMMEDIATELY+#','SplitSubtrees')],1,min)
}

if(memoryOption==3){
  wide.makespan$Select<-apply(wide.makespan[,c('AvoidChain','ImprovedSplit','LARGEST_FIT+#','SplitSubtrees')],1,min)
}

if(memoryOption==1){
  wide.makespan$AvoidChain<-wide.makespan$AvoidChain/wide.makespan$FIRST_FIT
  wide.makespan$ImprovedSplit<-wide.makespan$ImprovedSplit/wide.makespan$FIRST_FIT
  wide.makespan$SplitSubtrees<-wide.makespan$SplitSubtrees/wide.makespan$FIRST_FIT
  wide.makespan$Select<-wide.makespan$Select/wide.makespan$FIRST_FIT
  
  if(stage2Option=="LargestFirst")
    wide.makespan$`LARGEST_FIT+#`<-wide.makespan$`LARGEST_FIT+#`/wide.makespan$FIRST_FIT
  if(stage2Option=="FirstFit")
    wide.makespan$`FIRST_FIT+#`<-wide.makespan$`FIRST_FIT+#`/wide.makespan$FIRST_FIT
  if(stage2Option=='Immediately')
    wide.makespan$`IMMEDIATELY+#`<-wide.makespan$`IMMEDIATELY+#`/wide.makespan$FIRST_FIT
}

if(memoryOption==3){
  wide.makespan$AvoidChain<-wide.makespan$AvoidChain/wide.makespan$s1
  wide.makespan$ImprovedSplit<-wide.makespan$ImprovedSplit/wide.makespan$s1
  wide.makespan$SplitSubtrees<-wide.makespan$SplitSubtrees/wide.makespan$s1
  wide.makespan$`LARGEST_FIT+#`<-wide.makespan$`LARGEST_FIT+#`/wide.makespan$s1
  wide.makespan$Select<-wide.makespan$Select/wide.makespan$s1
}

long.makespan<-melt(wide.makespan,id.vars = c("TreeName","NPR","CCR"),variable.name = "Heuristic",value.name = "makespan" )
#long.makespan<-long.makespan[!long.makespan$Heuristic%in%c("FirstFit","Select"),]
long.makespan<-long.makespan[!long.makespan$Heuristic%in%c("FIRST_FIT","s1"),]
long.makespan$NPR<-1/long.makespan$NPR

long.makespan$CCR<-as.factor(long.makespan$CCR)
long.makespan$NPR<-factor(long.makespan$NPR)
levels(long.makespan$NPR)[levels(long.makespan$NPR)=="1e-04"]<-"PNR = 1e-04"
levels(long.makespan$NPR)[levels(long.makespan$NPR)=="0.001"]<-"PNR = 0.001"
levels(long.makespan$NPR)[levels(long.makespan$NPR)=="0.01"]<-"PNR = 0.01"

library(naniar)
long.makespan$makespan[which(long.makespan$makespan<0)]=-1
long.makespan <- replace_with_na(long.makespan,replace = list(makespan=-1))
long.makespan<-long.makespan[!is.na(long.makespan$makespan),]

#levels(long.makespan$Heuristic)[levels(long.makespan$Heuristic=="FirstFit+#")]<-"FirstFit+\\&"
if(memoryOption==1){
  levels(long.makespan$Heuristic)[levels(long.makespan$Heuristic=="AvoidChain")]<-"AvoidChain+\\#"
  levels(long.makespan$Heuristic)[levels(long.makespan$Heuristic=="ImprovedSplit")]<-"ImprovedSplit+\\#"
  levels(long.makespan$Heuristic)[levels(long.makespan$Heuristic=="SplitSubtrees")]<-"SplitSubtrees+\\#"
}
if(memoryOption==3){
  levels(long.makespan$Heuristic)[levels(long.makespan$Heuristic=="AvoidChain")]<-"AvoidChain+\\&"
  levels(long.makespan$Heuristic)[levels(long.makespan$Heuristic=="ImprovedSplit")]<-"ImprovedSplit+\\&"
  levels(long.makespan$Heuristic)[levels(long.makespan$Heuristic=="SplitSubtrees")]<-"SplitSubtrees+\\&"
}

if(memoryOption==1){
  m<-"strictM"
}
if(memoryOption==3){
  m<-"looseM"
}
l<-c("ASAPnochain","ImprovedSplit","Sequence","SplitSubtrees","Select")

if(stage2Option=='LargestFirst'){
  cb_palette <- c(AvoidChain="#B79F00",
                  ImprovedSplit="#00BA38", SplitSubtrees="#619CFF",
                  'LARGEST_FIT+#'="#FFFFFF","Select"="#FF7417")
  long.makespan$Heuristic<-factor(long.makespan$Heuristic,levels = c("AvoidChain","ImprovedSplit","LARGEST_FIT+#","SplitSubtrees","Select"))
}

if(stage2Option=='FirstFit'){
  cb_palette <- c(AvoidChain="#B79F00",
                  ImprovedSplit="#00BA38", SplitSubtrees="#619CFF",
                  'FIRST_FIT+#'="#FFFFFF","Select"="#FF7417")
  long.makespan$Heuristic<-factor(long.makespan$Heuristic,levels = c("AvoidChain","ImprovedSplit","FIRST_FIT+#","SplitSubtrees","Select"))
}
  
if(stage2Option=='Immediately'){
  cb_palette <- c(AvoidChain="#B79F00",
                  ImprovedSplit="#00BA38", SplitSubtrees="#619CFF",
                  'IMMEDIATELY+#'="#FFFFFF","Select"="#FF7417")
  long.makespan$Heuristic<-factor(long.makespan$Heuristic,levels = c("AvoidChain","ImprovedSplit","IMMEDIATELY+#","SplitSubtrees","Select"))
}

if(memoryOption==1){
  yaxis_title="Makespan normalized to \\textbf{FirstFit}"
}
if(memoryOption==3){
  yaxis_title="Makespan normalized to \\textbf{SplitSubtrees}"
}

temp<-max(long.makespan$makespan[long.makespan$Heuristic=="LARGEST_FIT+#"])
long.makespan[which(long.makespan$makespan==temp),]
temp<-max(long.makespan$makespan[long.makespan$Heuristic=="SplitSubtrees"])
long.makespan[which(long.makespan$makespan==temp),]

aggregate(long.makespan$makespan,list(long.makespan$NPR,long.makespan$Heuristic),mean)

output.name<-paste("~/ChangjiangGou/MemCom/journal/figure_makespan_final_",m,"_",stage2Option,".tex",sep = "")
output.name
#tikz(output.name,width = 7.2,height = 3)
ggplot(long.makespan,aes(x=factor(CCR),y=makespan,fill=factor(Heuristic)))+
  geom_boxplot(outlier.size = 0.05)+
  scale_x_discrete(name="CCR")+
  #coord_cartesian(ylim=c(0,2))+
  scale_y_continuous(name=yaxis_title)+
  scale_fill_manual(labels=l, values = cb_palette)+
  labs(fill="")+
  theme(legend.position = c(0.8,0.82),legend.background = element_blank(),legend.key = element_blank())+
  facet_grid(.~NPR)#+
  #geom_hline(yintercept=0.37)
dev.off()

#---------scheduling time of heuristics, figure 9/ 13/ 16/ 19----------------------
rm(list=ls())
result.firstStage<-read.table("./stage1/result_phase1.txt",header = TRUE)
result.firstStage$Heuristic<-paste(result.firstStage$Stage1,result.firstStage$Stage2,result.firstStage$Stage3,sep = "+")
time.firstStage<-dcast(result.firstStage,TreeName+NPR+CCR+AmountProcessors~Heuristic,value.var = "TimeConsuming")
time.firstStage$`ASAP+AvoidChain`<-time.firstStage$`AvoidChain+NA+NA`+time.firstStage$`ASAP+NA+NA`
time.firstStage<-time.firstStage[time.firstStage$CCR%in%c(0.1,1,10),]

result.twothreeStage<-read.table("./stage23/result_twothreeStage.txt",header = TRUE)
result.twothreeStage$afterS3<-!is.na(result.twothreeStage$Stage3)
result.twothreeStage<-result.twothreeStage[!is.na(result.twothreeStage$MemoryConstraint),]
result.twothreeStage$Heuristic<-paste(result.twothreeStage$Stage1,result.twothreeStage$Stage2,result.twothreeStage$afterS3,sep="+")
result.twothreeStage<-result.twothreeStage[result.twothreeStage$CCR%in%c(0.1,1,10),]
table(result.twothreeStage$Heuristic)

memoryOption<-3 #1 is the maxoutd, the strict case; 3 is the minmem, the loose case
stage2Option<-"LargestFirst" # stage 2 method: FirstFit, Immediately, LargestFirst

result.twothreeStage<-result.twothreeStage[result.twothreeStage$MemoryConstraint==memoryOption,]
table(result.twothreeStage$Heuristic)
result.twothreeStage<-result.twothreeStage[,c("TreeName","NPR","CCR","TimeConsuming","Heuristic")]
time.twothreeStage<-dcast(result.twothreeStage,TreeName+NPR+CCR~Heuristic,value.var = "TimeConsuming")

time<-merge(time.firstStage,time.twothreeStage,by=c("TreeName","NPR","CCR"))
colnames(time)

if(stage2Option=="LargestFirst"){
    time$`ASAP+AvoidChain+LargestFirst`<-time$`AvoidChain+LargestFirst+FALSE`+time$`ASAP+AvoidChain`
    time$`ASAP+AvoidChain+LargestFirst+Merge/LarSav`<-time$`ASAP+AvoidChain+LargestFirst`+time$`AvoidChain+LargestFirst+TRUE`
    time$`IS+LargestFirst`<-time$`ImprovedSplit+LargestFirst+FALSE`+time$`ImprovedSplit+NA+NA`
    time$`IS+LargestFirst+Merge/LarSav`<-time$`IS+LargestFirst`+time$`ImprovedSplit+LargestFirst+TRUE`
    time$`Split+LargestFirst`<-time$`SplitSubtrees+LargestFirst+FALSE`+time.firstStage$`SplitSubtrees+NA+NA`
    time$`Split+LargestFirst+Merge/LarSav`<-time$`Split+LargestFirst`+time$`SplitSubtrees+LargestFirst+TRUE`
}

if(stage2Option=="FirstFit"){
  time$`ASAP+AvoidChain+FirstFit`<-time$`AvoidChain+FirstFit+FALSE`+time$`ASAP+AvoidChain`
  time$`ASAP+AvoidChain+FirstFit+Merge/LarSav`<-time$`ASAP+AvoidChain+FirstFit`+time$`AvoidChain+FirstFit+TRUE`
  time$`IS+FirstFit`<-time$`ImprovedSplit+FirstFit+FALSE`+time$`ImprovedSplit+NA+NA`
  time$`IS+FirstFit+Merge/LarSav`<-time$`IS+FirstFit`+time$`ImprovedSplit+FirstFit+TRUE`
  time$`Split+FirstFit`<-time$`SplitSubtrees+FirstFit+FALSE`+time.firstStage$`SplitSubtrees+NA+NA`
  time$`Split+FirstFit+Merge/LarSav`<-time$`Split+FirstFit`+time$`SplitSubtrees+FirstFit+TRUE`
}

if(stage2Option=="Immediately"){
  time$`ASAP+AvoidChain+Immediately`<-time$`AvoidChain+Immediately+FALSE`+time$`ASAP+AvoidChain`
  time$`ASAP+AvoidChain+Immediately+Merge/LarSav`<-time$`ASAP+AvoidChain+Immediately`+time$`AvoidChain+Immediately+TRUE`
  time$`IS+Immediately`<-time$`ImprovedSplit+Immediately+FALSE`+time$`ImprovedSplit+NA+NA`
  time$`IS+Immediately+Merge/LarSav`<-time$`IS+Immediately`+time$`ImprovedSplit+Immediately+TRUE`
  time$`Split+Immediately`<-time$`SplitSubtrees+Immediately+FALSE`+time.firstStage$`SplitSubtrees+NA+NA`
  time$`Split+Immediately+Merge/LarSav`<-time$`Split+Immediately`+time$`SplitSubtrees+Immediately+TRUE`
}

filelist<-list.files("./version3/",pattern = "result*", full.names = TRUE)
temp<-lapply(filelist,read.table,header=TRUE,stringsAsFactors=FALSE)
result.firstfit<-Reduce(rbind,temp)

result.firstfit<-result.firstfit[result.firstfit$MemoryConstraint==memoryOption,]
result.firstfit<-result.firstfit[result.firstfit$CCR%in%c(0.1,1,10),]
colnames(result.firstfit)
colnames(result.firstfit)[8]<-"Heuristic"
result.firstfit$Heuristic[which(result.firstfit$Heuristic!="FirstFit")]<-"FirstFit+Merge/LarSav"
time.firstfit<-dcast(result.firstfit,TreeName+NPR+CCR+AmountProcessors~Heuristic,value.var = "TimeConsuming")
time.firstfit$`FirstFit+Merge/LarSav`<-time.firstfit$FirstFit+time.firstfit$`FirstFit+Merge/LarSav`

sequence.data<-read.table("./sequence/result_sequence.txt",header = TRUE)
data.subset<-subset(sequence.data,MemoryConstraint==memoryOption&CCR%in%c(0.1,1,10))
if(stage2Option=="LargestFirst"){
  data.subset<-subset(data.subset,Heuristic%in%c("LARGEST_FIT",
               "LARGEST_FIT+Merge","LARGEST_FIT+NA","LARGEST_FIT+SplitAgain"))
}
if(stage2Option=="FirstFit"){
  data.subset<-subset(data.subset,Heuristic%in%c("FIRST_FIT","FIRST_FIT+Merge",
    "FIRST_FIT+NA","FIRST_FIT+SplitAgain"))
}
if(stage2Option=="IMMEDIATELY"){
  data.subset<-subset(data.subset,Heuristic%in%c("IMMEDIATELY","IMMEDIATELY+Merge",
                                                 "IMMEDIATELY+NA","IMMEDIATELY+SplitAgain"))
}

data.wide<-dcast(data.subset,TreeName+NPR+CCR~Heuristic,value.var = "TimeConsuming")

if(stage2Option=="LargestFirst"){
  if(memoryOption==1){
    data.wide$`LARGEST_FIT+Merge`[is.na(data.wide$`LARGEST_FIT+Merge`)]<-0
    data.wide$`LARGEST_FIT+NA`[is.na(data.wide$`LARGEST_FIT+NA`)]<-0
    data.wide$`LARGEST_FIT+SplitAgain`[is.na(data.wide$`LARGEST_FIT+SplitAgain`)]<-0
    data.wide$'Sequence+LargestFirst+SplitAgain/Merge'<-data.wide$LARGEST_FIT+data.wide$`LARGEST_FIT+Merge`+
      data.wide$`LARGEST_FIT+NA`+data.wide$`LARGEST_FIT+SplitAgain`
  }
  if(memoryOption==3){
    data.wide$'Sequence+LargestFirst+SplitAgain/Merge'<-data.wide$LARGEST_FIT+data.wide$`LARGEST_FIT+SplitAgain`
  }
  
  time.sequence<-data.wide[,c("TreeName","NPR","CCR","Sequence+LargestFirst+SplitAgain/Merge")]
}

if(stage2Option=="FirstFit"){
  if(memoryOption==1){
    data.wide$`FIRST_FIT+Merge`[is.na(data.wide$`FIRST_FIT+Merge`)]<-0
    data.wide$`FIRST_FIT+NA`[is.na(data.wide$`FIRST_FIT+NA`)]<-0
    data.wide$`FIRST_FIT+SplitAgain`[is.na(data.wide$`FIRST_FIT+SplitAgain`)]<-0
    data.wide$'Sequence+FirstFit+SplitAgain/Merge'<-data.wide$FIRST_FIT+data.wide$`FIRST_FIT+Merge`+
      data.wide$`FIRST_FIT+NA`+data.wide$`FIRST_FIT+SplitAgain`
  }
  if(memoryOption==3){
    data.wide$'Sequence+FirstFit+SplitAgain/Merge'<-data.wide$FIRST_FIT+data.wide$`FIRST_FIT+SplitAgain`
  }
  time.sequence<-data.wide[,c("TreeName","NPR","CCR","Sequence+FirstFit+SplitAgain/Merge")]
}

if(stage2Option=="Immediately"){
  if(memoryOption==1){
    data.wide$`IMMEDIATELY+Merge`[is.na(data.wide$`IMMEDIATELY+Merge`)]<-0
    data.wide$`IMMEDIATELY+NA`[is.na(data.wide$`IMMEDIATELY+NA`)]<-0
    data.wide$`IMMEDIATELY+SplitAgain`[is.na(data.wide$`IMMEDIATELY+SplitAgain`)]<-0
    data.wide$'Sequence+Immediately+SplitAgain/Merge'<-data.wide$IMMEDIATELY+data.wide$`IMMEDIATELY+Merge`+
      data.wide$`IMMEDIATELY+NA`+data.wide$`IMMEDIATELY+SplitAgain`
  }
  if(memoryOption==3){
    data.wide$'Sequence+Immediately+SplitAgain/Merge'<-data.wide$IMMEDIATELY+data.wide$`IMMEDIATELY+SplitAgain`
  }
  time.sequence<-data.wide[,c("TreeName","NPR","CCR","Sequence+Immediately+SplitAgain/Merge")]
}

time<-merge(time,time.sequence)

if(stage2Option=="LargestFirst"){
  time$Select<-apply(time[,c('ASAP+AvoidChain+LargestFirst+Merge/LarSav','IS+LargestFirst+Merge/LarSav','Split+LargestFirst+Merge/LarSav','Sequence+LargestFirst+SplitAgain/Merge')],1,sum)
}

if(stage2Option=="FirstFit"){
  time$Select<-apply(time[,c('ASAP+AvoidChain+FirstFit+Merge/LarSav','IS+FirstFit+Merge/LarSav','Split+FirstFit+Merge/LarSav','Sequence+FirstFit+SplitAgain/Merge')],1,sum)
}

if(stage2Option=="Immediately"){
  time$Select<-apply(time[,c('ASAP+AvoidChain+Immediately+Merge/LarSav','IS+Immediately+Merge/LarSav','Split+Immediately+Merge/LarSav','Sequence+Immediately+SplitAgain/Merge')],1,sum)
}

CLOCKS_PER_SEC = 1e6
CLOCKS_PER_MIN = CLOCKS_PER_SEC*60
CLOCKS_PER_HOUR = CLOCKS_PER_SEC*3600

if(stage2Option=="LargestFirst"){
  time$`ASAP+AvoidChain+LargestFirst+Merge/LarSav`<-time$`ASAP+AvoidChain+LargestFirst+Merge/LarSav`/CLOCKS_PER_MIN
  time$`IS+LargestFirst+Merge/LarSav`<-time$`IS+LargestFirst+Merge/LarSav`/CLOCKS_PER_MIN
  time$`Split+LargestFirst+Merge/LarSav`<-time$`Split+LargestFirst+Merge/LarSav`/CLOCKS_PER_MIN
  time$`Sequence+LargestFirst+SplitAgain/Merge`<-time$`Sequence+LargestFirst+SplitAgain/Merge`/CLOCKS_PER_MIN
}

if(stage2Option=="FirstFit"){
  time$`ASAP+AvoidChain+FirstFit+Merge/LarSav`<-time$`ASAP+AvoidChain+FirstFit+Merge/LarSav`/CLOCKS_PER_MIN
  time$`IS+FirstFit+Merge/LarSav`<-time$`IS+FirstFit+Merge/LarSav`/CLOCKS_PER_MIN
  time$`Split+FirstFit+Merge/LarSav`<-time$`Split+FirstFit+Merge/LarSav`/CLOCKS_PER_MIN
  time$`Sequence+FirstFit+SplitAgain/Merge`<-time$`Sequence+FirstFit+SplitAgain/Merge`/CLOCKS_PER_MIN
}

if(stage2Option=="Immediately"){
  time$`ASAP+AvoidChain+Immediately+Merge/LarSav`<-time$`ASAP+AvoidChain+Immediately+Merge/LarSav`/CLOCKS_PER_MIN
  time$`IS+Immediately+Merge/LarSav`<-time$`IS+Immediately+Merge/LarSav`/CLOCKS_PER_MIN
  time$`Split+Immediately+Merge/LarSav`<-time$`Split+Immediately+Merge/LarSav`/CLOCKS_PER_MIN
  time$`Sequence+Immediately+SplitAgain/Merge`<-time$`Sequence+Immediately+SplitAgain/Merge`/CLOCKS_PER_MIN
}

time$Select<-time$Select/CLOCKS_PER_MIN

time$NPR<-1/time$NPR

long<-melt(time,id.vars = c("TreeName","NPR","CCR","AmountProcessors"),variable.name = "Heuristic",value.name = "time" )

if(stage2Option=="LargestFirst"){
  long.subset<-subset(long,Heuristic%in%c('Sequence+LargestFirst+SplitAgain/Merge','ASAP+AvoidChain+LargestFirst+Merge/LarSav','IS+LargestFirst+Merge/LarSav','Split+LargestFirst+Merge/LarSav','Select'))
}

if(stage2Option=="FirstFit"){
  long.subset<-subset(long,Heuristic%in%c('Sequence+FirstFit+SplitAgain/Merge','ASAP+AvoidChain+FirstFit+Merge/LarSav','IS+FirstFit+Merge/LarSav','Split+FirstFit+Merge/LarSav','Select'))
}

if(stage2Option=="Immediately"){
  long.subset<-subset(long,Heuristic%in%c('Sequence+Immediately+SplitAgain/Merge','ASAP+AvoidChain+Immediately+Merge/LarSav','IS+Immediately+Merge/LarSav','Split+Immediately+Merge/LarSav','Select'))
}

long.subset$NPR<-as.factor(long.subset$NPR)
levels(long.subset$NPR)[levels(long.subset$NPR)=="1e-04"]<-"PNR = 1e-04"
levels(long.subset$NPR)[levels(long.subset$NPR)=="0.001"]<-"PNR = 0.001"
levels(long.subset$NPR)[levels(long.subset$NPR)=="0.01"]<-"PNR = 0.01"

#long.subset$time[which(long.subset$time==0)]=1 

if(stage2Option=="LargestFirst"){
  cb_palette <- c("ASAP+AvoidChain+LargestFirst+Merge/LarSav"="#B79F00",
                  "IS+LargestFirst+Merge/LarSav"="#00BA38", "Split+LargestFirst+Merge/LarSav"="#619CFF",
                  "Sequence+LargestFirst+SplitAgain/Merge"="#FFFFFF","Select"="#FF7417")
  long.subset$Heuristic<-factor(long.subset$Heuristic,levels = c("ASAP+AvoidChain+LargestFirst+Merge/LarSav","IS+LargestFirst+Merge/LarSav","Sequence+LargestFirst+SplitAgain/Merge","Split+LargestFirst+Merge/LarSav","Select"))
}

if(stage2Option=="FirstFit"){
  cb_palette <- c("ASAP+AvoidChain+FirstFit+Merge/LarSav"="#B79F00",
                  "IS+FirstFit+Merge/LarSav"="#00BA38", "Split+FirstFit+Merge/LarSav"="#619CFF",
                  "Sequence+FirstFit+SplitAgain/Merge"="#FFFFFF","Select"="#FF7417")
  long.subset$Heuristic<-factor(long.subset$Heuristic,levels = c("ASAP+AvoidChain+FirstFit+Merge/LarSav","IS+FirstFit+Merge/LarSav","Sequence+FirstFit+SplitAgain/Merge","Split+FirstFit+Merge/LarSav","Select"))
}

if(stage2Option=="Immediately"){
  cb_palette <- c("ASAP+AvoidChain+Immediately+Merge/LarSav"="#B79F00",
                  "IS+Immediately+Merge/LarSav"="#00BA38", "Split+Immediately+Merge/LarSav"="#619CFF",
                  "Sequence+Immediately+SplitAgain/Merge"="#FFFFFF","Select"="#FF7417")
  long.subset$Heuristic<-factor(long.subset$Heuristic,levels = c("ASAP+AvoidChain+Immediately+Merge/LarSav","IS+Immediately+Merge/LarSav","Sequence+Immediately+SplitAgain/Merge","Split+Immediately+Merge/LarSav","Select"))
}

if(memoryOption==1){
  m<-"strictM"
}
if(memoryOption==3){
  m<-"looseM"
}

l<-c("ASAPnochain","ImprovedSplit","Sequence","SplitSubtrees","Select")

output.name<-paste("~/ChangjiangGou/MemCom/journal/figure_schedulingTime_",m,"_",stage2Option,".tex",sep = "")
output.name
#tikz(output.name,width = 7.2,height = 4)
ggplot(long.subset,aes(x=factor(CCR),y=time,fill=Heuristic))+geom_boxplot(outlier.size = 0.3)+
  scale_x_discrete(name="CCR")+
  scale_y_continuous(name="Logarithmic scheduling time of each case (minute)",trans = "log10")+
  scale_fill_manual(labels=l,values = cb_palette)+
  guides(fill=guide_legend(nrow = 1))+
  theme(legend.position = "top",axis.text.y = element_text(angle=70, hjust=0.5, vjust=0),legend.background = element_blank(),legend.key = element_blank())+
  labs(fill="")+
  facet_grid(.~NPR)#+
  #guides(fill=FALSE)
dev.off()

#-----------------scheduling time, as a function of tree size, figure 10/ 12 /19

tsize<-read.table("tree_size.txt",header = TRUE)
tsize<-tsize[,c("name","size")]
long.subset<-merge(long.subset,tsize,by.x = "TreeName",by.y = "name")

if(memoryOption==1){
  m<-"strictM"
}
if(memoryOption==3){
  m<-"looseM"
}

l<-c("ASAPnochain","ImprovedSplit","Sequence","SplitSubtrees","Select")

if(stage2Option=="LargestFirst"){
  cb_palette <- c("ASAP+AvoidChain+LargestFirst+Merge/LarSav"="#B79F00",
                  "IS+LargestFirst+Merge/LarSav"="#00BA38", "Split+LargestFirst+Merge/LarSav"="#619CFF",
                  "Sequence+LargestFirst+SplitAgain/Merge"="#FFD300","Select"="#FF7417")
}

if(stage2Option=="FirstFit"){
  cb_palette <- c("ASAP+AvoidChain+FirstFit+Merge/LarSav"="#B79F00",
                  "IS+FirstFit+Merge/LarSav"="#00BA38", "Split+FirstFit+Merge/LarSav"="#619CFF",
                  "Sequence+FirstFit+SplitAgain/Merge"="#FFD300","Select"="#FF7417")
}

if(stage2Option=="Immediately"){
  cb_palette <- c("ASAP+AvoidChain+Immediately+Merge/LarSav"="#B79F00",
                  "IS+Immediately+Merge/LarSav"="#00BA38", "Split+Immediately+Merge/LarSav"="#619CFF",
                  "Sequence+Immediately+SplitAgain/Merge"="#FFD300","Select"="#FF7417")
}

output.name<-paste("~/ChangjiangGou/MemCom/journal/figure_schedulingTime_size_",m,"_",stage2Option,".tex",sep = "")
output.name

#tikz(output.name,width = 7.2,height = 4)
ggplot(long.subset,aes(x=size,y=time,colour=Heuristic))+geom_point(size=0.7)+
  scale_x_continuous(name="Tree size",breaks = c(0,200000,500000,900000),labels = c(0,200*10^3,500*10^3,900*10^3))+
  #scale_y_continuous(name="Execution time normalized to \\textbf{Select}")+
  scale_y_continuous(name="Logarithmic scheduling time of each case (minute)", trans = "log10")+
  scale_color_manual(labels=l,values = cb_palette)+
  labs(colour="")+
  theme(legend.position = "top",axis.text.y = element_text(angle=70, hjust=0.5, vjust=0),axis.text.x = element_text(angle = 30, hjust = 1,vjust = 1),legend.background = element_blank(),legend.key = element_blank())+
  facet_grid(.~NPR)+
  guides(colour=guide_legend(nrow = 1))
dev.off()

#####################figure 9####################
output.name<-paste("~/ChangjiangGou/MemCom/journal/figure_schedulingTime_",m,"_",stage2Option,".tex",sep = "")
output.name
#tikz(output.name,width = 7.2,height = 4)
ggplot(long.subset,aes(x=size,y=time,colour=Heuristic))+geom_point(size = 0.7, position = position_dodge(width = 0.008))+
  scale_x_continuous(name="Logarithmic tree size",trans = "log10")+
  scale_y_continuous(name="Logarithmic scheduling time of each case (minute)",trans = "log10")+
  scale_color_manual(labels=l,values = cb_palette)+
  labs(colour="")+
  theme(legend.position = "top",axis.text.y = element_text(angle=70, hjust=0.5, vjust=0),axis.text.x = element_text(angle = 30, hjust = 1,vjust = 1),
        legend.background = element_blank(),legend.key = element_blank())+
  guides(fill=guide_legend(nrow = 1))
dev.off()

----------------------------------------------------------------------------------------------------
---------------------------------archive------------------------------------------------------------

#-------------makespan reduced and subtrees generated after stage2  
rm(list = ls())
rawdata<-read.table("./stage23/result_twothreeStage.txt",header = TRUE)
long<-rawdata[is.na(rawdata$Stage3),]#select stage 1 and stage 2
long<-long[long$MemoryConstraint==1|is.na(long$MemoryConstraint),]#select memory constraint as the maxoutd or NA

long$Heuristic<-paste(long$Stage1,"+",long$Stage2,sep="")
long<-long[,-c(4,6,8,9,10,11)]#remove coloum 
long.unique<-unique(long)

wide.ms<-dcast(long.unique,TreeName+NPR+CCR~Heuristic,value.var = "Makespan")
wide.amount<-dcast(long.unique,TreeName+NPR+CCR~Heuristic,value.var = "AmountSubtrees")

wide.amount$`AvoidChain+LargestFirst`<-wide.amount$`AvoidChain+LargestFirst`-wide.amount$`AvoidChain+NA`
wide.amount$`AvoidChain+Immediately`<-wide.amount$`AvoidChain+Immediately`-wide.amount$`AvoidChain+NA`
wide.amount$`AvoidChain+FirstFit`<-wide.amount$`AvoidChain+FirstFit`-wide.amount$`AvoidChain+NA`

wide.amount$`SplitSubtrees+LargestFirst`<-wide.amount$`SplitSubtrees+LargestFirst`-wide.amount$`SplitSubtrees+NA`
wide.amount$`SplitSubtrees+Immediately`<-wide.amount$`SplitSubtrees+Immediately`-wide.amount$`SplitSubtrees+NA`
wide.amount$`SplitSubtrees+FirstFit`<-wide.amount$`SplitSubtrees+FirstFit`-wide.amount$`SplitSubtrees+NA`

wide.amount$`ImprovedSplit+FirstFit`<-wide.amount$`ImprovedSplit+FirstFit`-wide.amount$`ImprovedSplit+NA`
wide.amount$`ImprovedSplit+Immediately`<-wide.amount$`ImprovedSplit+Immediately`-wide.amount$`ImprovedSplit+NA`
wide.amount$`ImprovedSplit+LargestFirst`<-wide.amount$`ImprovedSplit+LargestFirst`-wide.amount$`ImprovedSplit+NA`
#---
wide.ms$`AvoidChain+LargestFirst`<-wide.ms$`AvoidChain+LargestFirst`/wide.ms$`AvoidChain+NA`
wide.ms$`AvoidChain+Immediately`<-wide.ms$`AvoidChain+Immediately`/wide.ms$`AvoidChain+NA`
wide.ms$`AvoidChain+FirstFit`<-wide.ms$`AvoidChain+FirstFit`/wide.ms$`AvoidChain+NA`

wide.ms$`SplitSubtrees+LargestFirst`<-wide.ms$`SplitSubtrees+LargestFirst`/wide.ms$`SplitSubtrees+NA`
wide.ms$`SplitSubtrees+Immediately`<-wide.ms$`SplitSubtrees+Immediately`/wide.ms$`SplitSubtrees+NA`
wide.ms$`SplitSubtrees+FirstFit`<-wide.ms$`SplitSubtrees+FirstFit`/wide.ms$`SplitSubtrees+NA`

wide.ms$`ImprovedSplit+FirstFit`<-wide.ms$`ImprovedSplit+FirstFit`/wide.ms$`ImprovedSplit+NA`
wide.ms$`ImprovedSplit+Immediately`<-wide.ms$`ImprovedSplit+Immediately`/wide.ms$`ImprovedSplit+NA`
wide.ms$`ImprovedSplit+LargestFirst`<-wide.ms$`ImprovedSplit+LargestFirst`/wide.ms$`ImprovedSplit+NA`
#---
long.amount<-melt(wide.amount,id.vars = c("TreeName","NPR","CCR"),variable.name = "Heuristics",value.name = "amountIncreased" )
long.amount<-long.amount[!long.amount$Heuristics%in%c("AvoidChain+NA","ImprovedSplit+NA","SplitSubtrees+NA"),]

heuristics<-strsplit(as.character(long.amount$Heuristics),"\\+")
stage1<-unlist(lapply(seq(1:nrow(long.amount)),function(x){heuristics[[x]][[1]]}))
stage2<-unlist(lapply(seq(1:nrow(long.amount)),function(x){heuristics[[x]][[length(heuristics[[x]])]]}))
long.amount$stage1<-stage1
long.amount$stage2<-stage2
#---
long.ms<-melt(wide.ms,id.vars = c("TreeName","NPR","CCR"),variable.name = "Heuristics",value.name = "msRatio" )
long.ms<-long.ms[!long.ms$Heuristics%in%c("AvoidChain","ImprovedSplit","SplitSubtrees"),]

heuristics<-strsplit(as.character(long.ms$Heuristics),"\\+")
stage1<-unlist(lapply(seq(1:nrow(long.ms)),function(x){heuristics[[x]][[1]]}))
stage2<-unlist(lapply(seq(1:nrow(long.ms)),function(x){heuristics[[x]][[length(heuristics[[x]])]]}))
long.ms$stage1<-stage1
long.ms$stage2<-stage2
#---
long<-merge(long.ms,long.amount)
long$stage2<-factor(long$stage2,levels=c("FirstFit", "LargestFirst", "Immediately"))

tikz("~/ChangjiangGou/MemCom/journal/figure_stage1plus2_npr100.tex",width = 3.6,height = 3.6)
#tikz("~/ChangjiangGou/MemCom/journal/figure_stage1plus2_npr10000.tex",width = 3.6,height = 3.6)
ggplot(long[long$NPR==100,],aes(x=amountIncreased,y=msRatio,colour=Heuristics))+geom_point(size=0.7,shape=21)+
  facet_grid(stage1~stage2)+
  scale_x_continuous(name="Subtrees generated by fitting memory")+#,breaks = c(0,25,50,100))+
  scale_y_continuous(name="Makespan normalized to before fitting memory ",breaks = c(0.25,0.50,0.75,1.00))+
  coord_cartesian(ylim=c(0.15,1.1))+
  theme(axis.text.y = element_text(angle = 40, vjust=0.5))+
  guides(colour=FALSE)
#labs(colour="")+
#theme(legend.position = c(0.18,0.18),legend.background = element_blank(),legend.key = element_blank())
dev.off()


#--------------------------figure 8
tikz("~/ChangjiangGou/MemCom/journal/figure_amountSubtrees_density.tex",width = 3.6,height = 3)
ggplot(long[long$CCR==0.001,],aes(x=AmountSubtrees,colour=stage2))+geom_line(stat = "density")+
  scale_color_discrete(limits=c("FirstFit","LargestFirst","Immediately"))+
  labs(colour="",x="Number of subtrees")+
  theme(legend.position = c(0.15,0.8),legend.background = element_blank(),legend.key = element_blank())+
  facet_grid(.~NPR)
dev.off()

amountsubtrees<-aggregate(long$AmountSubtrees,list(long$stage2,long$NPR),sum)
colnames(amountsubtrees)<-c("stage2","npr","amount")
wide<-dcast(amountsubtrees,stage2~npr,value.var = "amount")
View(wide)#table 2 in the paper

#decrease in makespan
long<-rawdata[is.na(rawdata$Stage3),]#select stage 1 and stage 2
long<-long[long$MemoryConstraint==1|is.na(long$MemoryConstraint),]#select memory constraint as the maxoutd or NA

long$Heuristic<-paste(long$Stage1,"+",long$Stage2,sep="")
colnames(long)
long<-long[,c("TreeName","NPR","CCR","AmountSubtrees","AmountProcessors","Makespan","Heuristic")]
long<-unique(long)

wide.ms<-dcast(long,TreeName+NPR+CCR+AmountProcessors~Heuristic,value.var = "Makespan")
wide.ms<-wide.ms[wide.ms$CCR==0.001,]#select ccr 0.001

wide.ms$`AvoidChain+LargestFirst`<-wide.ms$`AvoidChain+LargestFirst`/wide.ms$`AvoidChain+NA`
wide.ms$`AvoidChain+Immediately`<-wide.ms$`AvoidChain+Immediately`/wide.ms$`AvoidChain+NA`
wide.ms$`AvoidChain+FirstFit`<-wide.ms$`AvoidChain+FirstFit`/wide.ms$`AvoidChain+NA`

wide.ms$`SplitSubtrees+LargestFirst`<-wide.ms$`SplitSubtrees+LargestFirst`/wide.ms$`SplitSubtrees+NA`
wide.ms$`SplitSubtrees+Immediately`<-wide.ms$`SplitSubtrees+Immediately`/wide.ms$`SplitSubtrees+NA`
wide.ms$`SplitSubtrees+FirstFit`<-wide.ms$`SplitSubtrees+FirstFit`/wide.ms$`SplitSubtrees+NA`

wide.ms$`ImprovedSplit+FirstFit`<-wide.ms$`ImprovedSplit+FirstFit`/wide.ms$`ImprovedSplit+NA`
wide.ms$`ImprovedSplit+Immediately`<-wide.ms$`ImprovedSplit+Immediately`/wide.ms$`ImprovedSplit+NA`
wide.ms$`ImprovedSplit+LargestFirst`<-wide.ms$`ImprovedSplit+LargestFirst`/wide.ms$`ImprovedSplit+NA`

colnames(wide.ms)
wide.ms<-wide.ms[,c("TreeName","NPR","AvoidChain+FirstFit", "AvoidChain+Immediately",
                    "AvoidChain+LargestFirst","ImprovedSplit+FirstFit","ImprovedSplit+Immediately",
                    "ImprovedSplit+LargestFirst","SplitSubtrees+FirstFit",
                    "SplitSubtrees+Immediately","SplitSubtrees+LargestFirst")]

long.ms<-melt(wide.ms,id.vars = c("TreeName","NPR"),variable.name = "Heuristics",value.name = "msRatio")

heuristics<-strsplit(as.character(long.ms$Heuristics),"\\+")
stage1<-unlist(lapply(seq(1:nrow(long.ms)),function(x){heuristics[[x]][[1]]}))
stage2<-unlist(lapply(seq(1:nrow(long.ms)),function(x){heuristics[[x]][[length(heuristics[[x]])]]}))
long.ms$stage1<-stage1
long.ms$stage2<-stage2

MS_ratio_mean<-aggregate(long.ms$msRatio,list(long.ms$NPR,long.ms$stage2),mean)
MS_ratio_mean$mean<-1-MS_ratio_mean$mean
colnames(MS_ratio_mean)<-c("NPR","Heuristic","mean")
View(MS_ratio_mean)# table 3 in the paper

#------
rm(list=ls())
rawdata<-read.table("./stage23/result_twothreeStage.txt",header = TRUE)
long<-rawdata

memoryOption<-1 #1 is the maxoutd, the strict case; 3 is the minmem, the loose case
stage2Option<-"Immediately" # stage 2 method: FirstFit, Immediately, LargestFirst

long$perc<-long$AmountSubtrees/long$AmountProcessors
long$NPR<-1/long$NPR
long<-long[!is.na(long$Stage2),]
long<-long[long$Stage2==stage2Option,]
long<-long[long$MemoryConstraint==memoryOption,]#choose the maxoutd option
long<-long[is.na(long$Stage3),]
long$Heuristic<-paste(long$Stage1,long$Stage2,sep="")

#long<-long[long$CCR==0.001,]

long$NPR<-as.factor(long$NPR)
levels(long$NPR)[levels(long$NPR)=="1e-04"]<-"PNR = 1e-04"
levels(long$NPR)[levels(long$NPR)=="0.001"]<-"PNR = 0.001"
levels(long$NPR)[levels(long$NPR)=="0.01"]<-"PNR = 0.01"

#-------------------amount of subtrees after stage2, figure 8/20/25
if(memoryOption==1){m<-"strictM"}
if(memoryOption==3){m<-"looseM"}

if(stage2Option=="LargestFirst"){
  cb_palette <- c(AvoidChainLargestFirst="#7CAE00",
                  ImprovedSplitLargestFirst="#00BFC4", SplitSubtreesLargestFirst="#C77CFF")
  l<-c("AvoidChain+LargestFirst","ImprovedSplit+LargestFirst","SplitSubtrees+LargestFirst")
}
if(stage2Option=="FirstFit"){
  cb_palette <- c(AvoidChainFirstFit="#7CAE00",
                  ImprovedSplitFirstFit="#00BFC4", SplitSubtreesFirstFit="#C77CFF")
  l<-c("AvoidChain+FirstFit","ImprovedSplit+FirstFit","SplitSubtrees+FirstFit")
}
if(stage2Option=="Immediately"){
  cb_palette <- c(AvoidChainImmediately="#7CAE00",
                  ImprovedSplitImmediately="#00BFC4", SplitSubtreesImmediately="#C77CFF")
  l<-c("AvoidChain+Immediately","ImprovedSplit+Immediately","SplitSubtrees+Immediately")
}

output.name<-paste("../../figure_amountSubtrees_",m,"_",stage2Option,".tex",sep = "")
output.name
#tikz(output.name,width = 3.6,height = 3)
ggplot(long,aes(x=factor(CCR),y=perc,fill=Heuristic))+geom_boxplot(outlier.size = 0.5)+
  scale_x_discrete(name="CCR")+
  scale_y_continuous(name="Ratio of number of subtrees to processors")+
  labs(fill="")+
  scale_fill_manual(labels=l,values = cb_palette)+
  theme(legend.position = c(0.72,0.85),legend.background = element_blank(),legend.key = element_blank())+
  coord_cartesian(ylim=c(0,4))+
  facet_grid(.~NPR)
dev.off()

