##Basic script describing the data

library(dplyr)
library(tidyr)
library(ggfortify)
source("R/pdcFileProcessing.R")



#' plot means from values
compareMeansToVals<-function(withmeans,thresh=3){
  require(ggplot2)
     #pivot_longer(cols=c(meanDiffs,logratio),names_to='Statistic',values_to='value')
  res=subset(withmeans,abs(meanDiffs)>thresh)
#  ggplot(withmeans,aes(x=disease,y=value,fill=Statistic))+geom_boxplot()+scale_y_log10()
  ggplot(res,aes(x=logratio,y=meanDiffs,col=disease))+geom_point()
  res
  }

#' cluster samples by disease
plotPCAbyDisease<-function(distab){
  library(viridis)
  
  mat=distab2Matrix(distab)
  dat<-distab%>%
    dplyr::select(Patient,disease)%>%distinct()
  
  gres<-autoplot(prcomp(mat,scale=T,center=T),data=dat,colour='disease')
  gres
}

distab2Matrix<-function(distab){
  mat<-distab%>%
    subset(!Gene%in%c('Mean','Median','StdDev'))%>%
    dplyr::select(-NCBIGeneID)%>%
    dplyr::select(Gene,Patient,disease,logratio)%>%
    tidyr::pivot_wider(names_from=Gene,values_from=logratio,values_fill=list(logratio=0))%>%
    tibble::column_to_rownames('Patient')
  xm=apply(mat[,-1],2,as.numeric)
  rownames(xm)<-rownames(mat)
  return(xm)
}

#'try to find predictive proteins using lasso regression
#'@param distab is a tidied datset of gene, logratio, and disease values
#'@return a tidied data frame of selected genes and the disease or which they are preidctive
#'@export
getSelectedProteins<-function(distab){
  library(glmnet)
  library(purrr)
  mat<-distab%>%
    subset(!Gene%in%c('Mean','Median','StdDev'))%>%
    dplyr::select(-NCBIGeneID)%>%
    dplyr::select(Gene,Patient,disease,logratio)%>%
    tidyr::pivot_wider(names_from=Gene,values_from=logratio,values_fill=list(logratio=0))%>%
    tibble::column_to_rownames('Patient')
  xm=apply(mat[,-1],2,as.numeric)
  rownames(xm)<-rownames(mat)
  cvfit=cv.glmnet(x=xm,y=mat$disease,family='multinomial')
  #extract proteins with non-zero co-efficiets
  coefs<-coef(cvfit,s='lambda.min')
  res<-map_df(coefs,as.matrix)%>%
    mutate(Gene=c('Intercept',colnames(xm)))%>%
    pivot_longer(cols=names(coefs),values_to='coef',names_to='disease')%>%
    subset(coef!=0)
   return(res)
  
}


plotHeatmapOfProts<-function(distab,prots,prefix=''){
  xm=distab2Matrix(distab)
  library(pheatmap)
  adf<-dplyr::select(distab,c(Patient,disease))%>%
    distinct()%>%column_to_rownames('Patient')
  
  pheatmap(t(xm[,intersect(colnames(xm),prots)]),
           annotation_col=adf,
           clustering_distance_cols='correlation',
           clustering_distance_rows='correlation',
           annotation_names_cols=FALSE,
           clustering_method = 'ward.D2',
           show_colnames=F,
           cellheight=9,
           filename=(paste0(prefix,paste(unique(distab$disease),collapse='_'),'.pdf')))
}

#####Build strawman set of proteins that distinguish a class 
#####Testing on ability to predict cancer subtype



testRegressionFromProteinAbundance<-function(){
  distab<-readInAllFiles()
  #lasso regression can select great protein markers
  res1=getSelectedProteins(distab)
  plotHeatmapOfProts(distab,res1$Gene,prefix='lassoSelectedProtFor')
  p1<-plotPCAbyDisease(subset(distab,Gene%in%res1$Gene))
  ggsave(p1,'lassoSelectedProts.png')
  
  #now if we filter by quantile
  quantiled=distab%>%group_by(Patient)%>%
    mutate(topThresh=quantile(logratio,0.99))%>%
    mutate(highExp=logratio>topThresh)%>%
    ungroup()
  
  qprots<-subset(quantiled,highExp==TRUE)%>%
    select(Gene)%>%distinct()
  
  quantiled%>%
    group_by(disease)%>%
    summarize(numProts=n_distinct(Gene))
  
  res2=getSelectedProteins(subset(quantiled,Gene%in%qprots$Gene))
  plotHeatmapOfProts(distab,res2$Gene,prefix='lassoSelectedTopQuantProtFor')
  p2<-plotPCAbyDisease(subset(distab,Gene%in%res2$Gene))
  ggsave(p2,'lassoSelectedTopQuantProts.png')
  }
#Fourth: Look they clsuter by subtype now!