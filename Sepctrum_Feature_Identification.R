{
  library(base)
  library(vegan)
  library(fossil)
  library(ape)
  library(GUniFrac)
  library(data.table)
  library(bnlearn)
  library(Rgraphviz)
  library(Metrics)
  library(tidyverse)
  library(ggplot2)
  library(e1071)
  library(rpart)
  library(rpart.plot)
  library(maptree)
  library(xgboost)
}
#---------------------

### feature identification
{
  setwd("C:/Users/CZ100/Desktop/Project_HU/Spectrum")
  list<-list.files()
  ###The features include number of peaks, roughness, and the highest peak of 8 ranges
  ftable<-matrix(nrow=1,ncol=8*3+7,0)
  ftable_notrans<-matrix(nrow=1,ncol=437+7,0)
    
  coren<-list
  coren<-gsub('[ ]','_',coren)
  coren<-sub('.csv','',coren)
  cn<-rep("",8*3+7)
  cn_notrans<-rep("",437+7)
  for (i in 1:8)
  {
    cn[i*3-2]=paste("Npeaks_R",i,sep="")
    cn[i*3-1]=paste("Roughness_R",i,sep="")
    cn[i*3]=paste("Peakhign_R",i,sep="")
  }
  cn[25:31]<-coren
  cn_notrans[1:437]=1:437
  cn_notrans[438:444]=coren
  colnames(ftable)<-cn
  
  ftable<-ftable[-1,,drop=F]
  ftable_notrans<-ftable_notrans[-1,,drop=F]
  
  ## The features include number of peaks, roughness, and the highest peak of 8 ranges
  ## R1:940-980   R2:980-1020   R3:1120-1200   R4:1260-1340
  ## R5:400-1480  R6:1500-1540  R7:1580-1620   R8:1620-1680
  range_lower=c(940,980,1120,1260,1400,1500,1580,1620)
  range_upper=c(980,1020,1200,1340,1480,1540,1620,1680)
  js_shift=0
  for(k in 1:8)
    js_shift=js_shift+length(which(table_spec$V1>=range_lower[k] & table_spec$V1<=range_upper[k]))
  
  
  for (i in 1:7)
  {
    #i=1
    table_spec<-read.csv(file=list[i],header=F,check.names=F)
    
    for (j in 1:60)
    {
      #j=1
      rec<-matrix(nrow=1,ncol=8*3+7,0)
      rec_notrans<-NULL
      for (k in 1:8)
      {
        #k=8
        rec_range=table_spec[which(table_spec$V1>=range_lower[k] & table_spec$V1<=range_upper[k]),j+1,drop=F]
        rec_notrans<-cbind(rec_notrans,t(rec_range))
        ### peak identification
        l<-nrow(rec_range)
        js_peak=0
        for (q in 6:(l-5))
          if (rec_range[q,1]>rec_range[q-1,1] & rec_range[q,1]>rec_range[q+1,1])
          {
            diff_left<-rec_range[q,1]-rec_range[q-5,1]
            diff_right<-rec_range[q,1]-rec_range[q+5,1]
            if (diff_left>200 & diff_right>200)
              js_peak=js_peak+1
          }
            
        rec[1,k*3-2]=js_peak
        ### roughness (STD)
        rec[1,k*3-1]=sd(rec_range[,1])
        ### highest peak
        rec[1,k*3]=max(rec_range[,1])
      }
      ### give label
      rec[1,24+i]=1
      if (length(rec_notrans)<437) rec_notrans<-cbind(rec_notrans,rep(rec_notrans[1,length(rec_notrans)],437-length(rec_notrans))) 
      rec_notrans=cbind(rec_notrans,rec[,25:31,drop=F])
      ftable<-rbind(ftable,rec)
      ftable_notrans<-rbind(ftable_notrans,rec_notrans)
    }
  }
  colnames(ftable)<-cn
  colnames(ftable_notrans)<-cn_notrans
  row.names(ftable_notrans)<-1:420
} 

setwd("C:/Users/CZ100/Desktop/Project_HU/")
write.csv(ftable,"feature_table.csv")
write.csv(ftable_notrans,"feature_notransformation_table.csv")


### Set dividing

{
  ftable<-as.data.frame(ftable)
  ftable_notrans<-as.data.frame(ftable_notrans)
  
  ftable$Species=0
  ftable_notrans$Species=0
  
  for (i in 1:420)
  {
    ftable$Species[i]=coren[which(ftable[i,25:31]==1)]
    ftable_notrans$Species[i]=coren[which(ftable[i,25:31]==1)]
  }
  ftable_ml<-ftable[,-(25:31)]
  ftable_notrans_ml<-ftable_notrans[,-(438:444)]
  
  ## Scaling
  ftable_ml_input<-ftable_ml[,-25]
  par1=seq(from=1,to=22,by=3)
  par2=seq(from=2,to=23,by=3)
  par3=seq(from=3,to=24,by=3)
  ## caution: detect negative values of the peak of Erysiphe cichoracearum (4th core)
  {
    maxs <- max(ftable_ml_input[,par1]) 
    mins <- min(ftable_ml_input[,par1])
    for (i in 1:8)
      ftable_ml_input[,par1[i]] <- as.data.frame(scale(ftable_ml_input[,par1[i]], center = mins, scale = maxs - mins)) # scaling to 0 - 1
    
    maxs <- max(ftable_ml_input[,par2]) 
    mins <- min(ftable_ml_input[,par2])
    for (i in 1:8)
      ftable_ml_input[,par2[i]] <- as.data.frame(scale(ftable_ml_input[,par2[i]], center = mins, scale = maxs - mins)) # scaling to 0 - 1
    
    maxs <- max(ftable_ml_input[,par3]) 
    mins <- min(ftable_ml_input[,par3])
    for (i in 1:8)
      ftable_ml_input[,par3[i]] <- as.data.frame(scale(ftable_ml_input[,par3[i]], center = mins, scale = maxs - mins)) # scaling to 0 - 1
  }
  
  ftable_ml[,1:24]<-ftable_ml_input
  
  ## train:test=4:1
  loc_base=c(49:60)
  loc_testing=loc_base
  for (i in 2:7)
    loc_testing=c(loc_testing,loc_base+60*(i-1))
  
  ftable_test=ftable_ml[loc_testing,]
  ftable_train=ftable_ml[-loc_testing,]
  
  ftable_notrans_test=ftable_notrans_ml[loc_testing,]
  ftable_notrans_train=ftable_notrans_ml[-loc_testing,]
  
  nofs_train<-336
  nofs_testing<-84
}
