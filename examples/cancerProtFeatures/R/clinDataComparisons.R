##compare communities to clinical variables
require(tidyr)
require(dplyr)

clin.dat<- read.csv('data/PDC_clinical_manifest_07182020_150316.csv',header=T,stringsAsFactors=F)%>%
  rename(Sample="ï..Case.ID")


comDist<-read.csv('panCancerDistances.csv',header=T)%>%
  subset(net2_type=='community')%>%distinct()%>%
  rename(Community='hyp2',Sample='net1')



createDf<-function(clin.dat,comDist,disease){
  dists=comDist%>%subset(Community==disease)%>%
  pivot_wider(-c(X,net1_type,net2_type),values_from=distance,names_from=net2)%>%
    left_join(clin.dat)
  
}