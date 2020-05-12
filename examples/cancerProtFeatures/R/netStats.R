
library(ggplot2)
library(dplyr)
visualizeCommunitySize<-function(fname){
  tab<-read.csv(fname)%>%
    group_by(Partition)%>%
    summarize(comSize=n_distinct(Node))
  print(tab)
  tab

}

forestGOStats<-function(topQuant='GOenrichmentPerPatient.tsv',
                        forestVals=list(luad='luadenrichedForestGoTerms.csv',
                                        brca='brcaenrichedForestGoTerms.csv',
                                        gbm='gbmenrichedForestGoTerms.csv',
                                        coad='coadenrichedForestGoTerms.csv')){
  forestGo=do.call(rbind,lapply(names(forestVals),function(x){
    print(x)
    res<-read.csv2(forestVals[[x]],sep=',',header=T,stringsAsFactors=F)
    res$disease=rep(x,nrow(res))
    return(res)
  }))%>%
    subset(namespace='biological_process')
  
  quantGo=read.csv2(topQuant,sep='\t',header=T,stringsAsFactors = F)
  
  patStats<-quantGo%>%select(disease,Patient,quantTerm='ID')%>%
    rowwise()%>%
    mutate(Patient=stringr::str_replace(Patient," Log Ratio",""))

  
  netStats<-forestGo%>%select(disease,Patient,netTerm='term')%>%
    mutate(disease=toupper(disease))
  ##now compare how many go terms are unique,distinct and shared for each patient in each disease
  
  fullStats<-patStats%>%full_join(netStats,by=c('Patient','disease'))
  
  go.ver<-fullStats%>%
    group_by(disease,Patient)%>%
    summarize(quantOnly=length(setdiff(quantTerm,netTerm)),netOnly=length(setdiff(netTerm,quantTerm)),agreement=length(intersect(quantTerm,netTerm)))%>%select(disease,Patient,quantOnly,netOnly,agreement)%>%
    distinct()
  
  ggplot(go.ver)+
    geom_point(aes(x=quantOnly,y=netOnly,size=agreement,col=disease))
  ggsave('goEnrichmentOverlap.pdf',useDingbats=FALSE)
}


forStats<-function(forestFiles=list(luad='luad_forestStats.csv',
                                      brca='brca_forestStats.csv', 
                                        coad='coad_forestStats.csv',
                                             gbm='gbm_forestStats.csv')){
      
    full.stats<-do.call(rbind,lapply(names(forestFiles),function(x){
      read.csv2(forestFiles[[x]],sep=',',header=T,stringsAsFactors=F)%>%
        mutate(disease=x)%>%
        mutate(ForestPrize=as.numeric(ForestPrize))%>%
        mutate(DisWeight=as.numeric(DisWeight))
    }))

    ##first compute the correlation of weights
    corStats<-full.stats%>%
      group_by(disease,Patient)%>%
      mutate(correlation=cor(as.numeric(DisWeight),as.numeric(ForestPrize)))%>%
      ungroup()%>%
      select(disease,Patient,correlation)%>%distinct()
    
    ggplot(corStats,aes(fill=disease))+
      geom_histogram(aes(x=correlation),position='dodge')+
        ggtitle('Correlation of disease weights to forest weights') 
    ggsave('proteinNetworkCorrelation.pdf')
    
    
    nodeType<-full.stats%>%#subset(ForestPrize!=0)%>%
      mutate(isSteiner=case_when(DisWeight==0 ~ TRUE,DisWeight>0 ~ FALSE))
    
    steins<-subset(nodeType,isSteiner==TRUE)%>%
      group_by(disease,Patient)%>%
      summarize(numSteiner=n_distinct(Gene))%>%
      select(disease,Patient,numSteiner)
    terms=subset(nodeType,isSteiner==FALSE)%>%
      subset(ForestPrize!=0)%>%
      group_by(disease,Patient)%>%
      summarize(numTerms=n_distinct(Gene))%>%
      select(disease,Patient,numTerms)
  
    nodeTypeVals<-steins%>%inner_join(terms,by=c('disease','Patient'))
    
    ggplot(nodeTypeVals,aes(col=disease))+geom_point(aes(x=numTerms,y=numSteiner))+ggtitle("Number of steiner vs terminals")
    ggsave('steinerVsTerminal.pdf',useDingbats=FALSE)
    
      
    full.stats<-corStats%>%
      left_join(nodeTypeVals)
    
    return(full.stats)
}

compareCommunitiesToForests<-function(communityFiles=list(luad='luadcommunities.csv',
                                                   gbm='gbmcommunities.csv',brca='brcacommunities.csv',
                                                   coad='coadcommunities.csv')){
    
  
  
}


plotNetworkDistances<-function(communityDistanceFile){
  library(cowplot)
  library(ggplot2)
  
  tab<-read.csv(communityDistanceFile,header=T,stringsAsFactors=F)
  
  mod.tab<-tab%>%rowwise()%>%
    mutate(net1_name=stringr::str_c(hyp1,net1,sep='_'))%>%
    mutate(net2_name=stringr::str_c(hyp2,net2,sep='_'))
  
  annotes<-mod.tab%>%
    dplyr::select(net1_name,net2_name,net1_type,net2_type,hyp1,hyp2)%>%
    distinct()
  
  ###we can probably remove this or alter it once we fix code
  red.annote<-annotes%>%
    dplyr::select(sample='net2_name',graph='net2_type',disease='hyp2')%>%
    distinct()
  
  ##remove this oncec we fix code
  missed<-subset(mod.tab,net2_name%in%setdiff(annotes$net2_name,annotes$net1_name))%>%
    dplyr::select(net1_name,net2_name,distance)%>%
    rename(net2_name='net1_name',net1_name='net2_name')


  ##first do each disease
  res=lapply(unique(mod.tab$hyp1),function(disease){
    red.tab<-subset(mod.tab,hyp1==disease)%>%
      subset(hyp2==disease)
    
    red.missed<-subset(red.tab,net2_name%in%setdiff(annotes$net2_name,annotes$net1_name))%>%
      dplyr::select(net1_name,net2_name,distance)%>%
      rename(net2_name='net1_name',net1_name='net2_name')
    
    as.mat <-red.tab%>%
      dplyr::select(net1_name,net2_name,distance)%>%
      rbind(red.missed)%>%
      tidyr::pivot_wider(values_from=distance,names_from=net2_name,values_fn=list(distance=mean),values_fill=list(distance=1.0))%>%
      tibble::column_to_rownames('net1_name')%>%
      as.matrix()
    
    
    ##calculate how many nodes in neighborhood
    neighborhood=red.tab%>%
        group_by(net2_name)%>%
        mutate(noVals=(distance==1.0))%>%
        subset(noVals==FALSE)%>%
        summarize(neighborhood=n_distinct(net1_name))%>%
    rename(sample='net2_name')
    
    res<-cmdscale(as.dist(as.mat))
    
    colnames(res)<-c('Dim1','Dim2')
    
    full.tab<-res%>%
      as.data.frame()%>%
      tibble::rownames_to_column("sample")%>%
      left_join(red.annote)%>%
      left_join(neighborhood)
    
    p <-ggplot(full.tab,aes(x=Dim1,y=Dim2,col=graph,shape=graph,size=neighborhood))+
      geom_point()+ggtitle(paste(disease,'Hyphae'))+theme_minimal()
    return(p)
  })
   
  ##now do all diseases combined
  cowplot::plot_grid(plotlist=res)
  ggsave('individualPlots.pdf',useDingbats=FALSE)
    
 
  
  ##then do combined
  as.mat<-mod.tab%>%
    dplyr::select(net1_name,net2_name,distance)%>%
    rbind(missed)%>%
    tidyr::pivot_wider(values_from=distance,names_from=net2_name,values_fn=list(distance=mean),values_fill=list(distance=1.0))%>%
    tibble::column_to_rownames('net1_name')%>%
    as.matrix()
  
  neighborhood=mod.tab%>%
    group_by(net2_name)%>%
    mutate(noVals=(distance==1.0))%>%
    subset(noVals==FALSE)%>%
    summarize(neighborhood=n_distinct(net1_name))%>%
    rename(sample='net2_name')
  
  
 
  res<-cmdscale(as.dist(as.mat))
  colnames(res)<-c('Dim1','Dim2')
  
  
  red.annote<-annotes%>%
    dplyr::select(sample='net2_name',graph='net2_type',disease='hyp2')%>%
    distinct()
  full.tab<-res%>%
    as.data.frame()%>%
    tibble::rownames_to_column("sample")%>%
    left_join(red.annote)%>%
    left_join(neighborhood)
  
  library(ggplot2)
  ggplot(full.tab,aes(x=Dim1,y=Dim2,col=disease,shape=graph,size=neighborhood))+geom_point()
  ggsave('allPlotsTogether.pdf',useDingbats=FALSE)
}


forStats() ##how well do the forests recapitulate biology?
plotNetworkDistances('panCancerDistances.csv') #how well do the communities summarize the diversity of the forests?
forestGOStats()##we get more function for patients, and it agrees with other stuff