##read and tidy files

library(dplyr)
library(tidyr)
#First: read in files

#' do the tidying of the data to create a basic patient/protein log ratio model
processSingleFile<-function(fname){
  tab<-read.csv2(fname,sep='\t',header=T,stringsAsFactors = F,check.names=F)
  unshared=grep('Unshared',colnames(tab))
  if(length(unshared)>0)
    tab<-tab[,-unshared]
  
  df<-tab%>%pivot_longer(-c(Gene,NCBIGeneID),names_to='Patient',values_to="logratio")
  print(paste('Processed file with',length(unique(df$Patient)),'patient samples'))
  df
  
}

readInAllFiles<-function(){
  dislist<-list(BRCA=processSingleFile('data/CPTAC2_Breast_Prospective_Collection_BI_Proteome.tmt10.tsv'),
                LUAD=processSingleFile('data/CPTAC3_Lung_Adeno_Carcinoma_Proteome.tmt10.tsv'),
                GBM = processSingleFile('data/CPTAC3_Glioblastoma_Multiforme_Proteome.tmt11.tsv'),
                COAD = processSingleFile('data/CPTAC2_Colon_Prospective_Collection_PNNL_Proteome.tmt10.tsv'))
  
  distab <-do.call(rbind,lapply(names(dislist),function(x) mutate(dislist[[x]],disease=x)))
  distab$logratio<-as.numeric(distab$logratio)
  distab<-subset(distab,!is.na(logratio))%>%
    subset(!Patient%in%c('POOL Log Ratio','Authority','Organism','Locus','Description','Chromosome'))%>%
    subset(!Gene%in%c('Mean','Median','StdDev'))
}
