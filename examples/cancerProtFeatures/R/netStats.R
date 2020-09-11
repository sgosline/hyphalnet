
library(ggplot2)
library(dplyr)
visualizeCommunitySize<-function(fname){
  tab<-read.csv(fname)%>%
    group_by(Partition)%>%
    summarize(comSize=n_distinct(Node))
  print(tab)
  tab

}

#'forestGOStats
#'Compares GO enrihcment of forests compared to patient data
forestGOStats<-function(topQuant='GOenrichmentPerPatientVsNormal.tsv'){
  
    forestGo<-subset(goStats,networkType=='Patient')

    quantGo=read.csv2(topQuant,sep='\t',header=T,stringsAsFactors = F)%>%dplyr::select(disease,Patient,quantTerm='ID')%>%
    rowwise()%>%
    mutate(Patient=stringr::str_replace(Patient," Log Ratio",""))

  
  netStats<-forestGo%>%dplyr::select(disease,Patient,netTerm='term')%>%
    mutate(disease=toupper(disease))
  ##now compare how many go terms are unique,distinct and shared for each patient in each disease
  
  fullStats<-quantGo%>%full_join(netStats,by=c('Patient','disease'))
  
  go.ver<-fullStats%>%
    group_by(disease,Patient)%>%
    summarize(quantOnly=length(setdiff(quantTerm,netTerm)),netOnly=length(setdiff(netTerm,quantTerm)),agreement=length(intersect(quantTerm,netTerm)))%>%
    dplyr::select(disease,Patient,quantOnly,netOnly,agreement)%>%
    distinct()
  
  ggplot(go.ver)+
    geom_point(aes(x=quantOnly,y=netOnly,size=agreement,col=disease))+scale_color_viridis_d()
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
      dplyr::select(disease,Patient,correlation)%>%distinct()
    
    ggplot(corStats,aes(fill=disease))+
      geom_histogram(aes(x=correlation),position='dodge')+
        ggtitle('Correlation of disease weights to forest weights') 
    ggsave('proteinNetworkCorrelation.pdf')
    
    
    nodeType<-full.stats%>%#subset(ForestPrize!=0)%>%
      mutate(isSteiner=case_when(DisWeight==0 ~ TRUE,DisWeight>0 ~ FALSE))
    
    steins<-subset(nodeType,isSteiner==TRUE)%>%
      group_by(disease,Patient)%>%
      summarize(numSteiner=n_distinct(Gene))%>%
      dplyr::select(disease,Patient,numSteiner)
    terms=subset(nodeType,isSteiner==FALSE)%>%
      subset(ForestPrize!=0)%>%
      group_by(disease,Patient)%>%
      summarize(numTerms=n_distinct(Gene))%>%
      dplyr::select(disease,Patient,numTerms)
  
    nodeTypeVals<-steins%>%inner_join(terms,by=c('disease','Patient'))
    
    ggplot(nodeTypeVals,aes(col=disease))+geom_point(aes(x=numTerms,y=numSteiner))+ggtitle("Number of steiner vs terminals")+scale_color_viridis_d()
    ggsave('steinerVsTerminal.pdf',useDingbats=FALSE)
    
      
    full.stats<-corStats%>%
      left_join(nodeTypeVals)
    
    return(full.stats)
}


assignClosestNodes<-function(communityDistanceFile){
  #process the distances
  tab<-read.csv(communityDistanceFile,header=T,stringsAsFactors=F)
  minDists=tab%>%subset(net2_type=='community')%>%
      group_by(net1)%>%
      filter(distance==min(distance))%>%
    dplyr::select(graph='net1',disease='hyp1',community='net2',commDis='hyp2',distance)
    
  minDists%>%
      group_by(disease,commDis)%>%
      summarize(numNets=n_distinct(graph))%>%
      ggplot(aes(x=disease,y=numNets,fill=commDis))+
    geom_bar(stat='identity',position='dodge')+scale_fill_viridis_d()
  write.csv(minDists,file='networksAndClosestCommunity.csv')
  ggsave('networksAssignedToClosestCommunity.pdf')  
  return(minDists)
}


getGoValues<-function(communityVals=list(luad='luad0.01enrichedCommunityGOterms.csv',
                                        brca='brca0.01enrichedCommunityGOterms.csv',
                                        gbm='gbm0.01enrichedCommunityGOterms.csv',
                                        coad='coad0.01enrichedCommunityGOterms.csv'),
                     forestVals=list(luad='luad0.01enrichedForestGoTerms.csv',
                                     brca='brca0.01enrichedForestGoTerms.csv',
                                     gbm='gbm0.01enrichedForestGoTerms.csv',
                                     coad='coad0.01enrichedForestGoTerms.csv')){
  
  commGo=do.call(rbind,lapply(names(communityVals),function(x){
    print(x)
    res<-read.csv2(communityVals[[x]],sep=',',header=T,stringsAsFactors=F)
    res$disease=rep(x,nrow(res))
    return(res)
  }))%>%
    subset(namespace='biological_process')%>%
    select(term,q,name,Community,disease)%>%
    tidyr::pivot_longer(Community,names_to='networkType',values_to='sample')
  
  ##read in all enrichment files - those from individual forests, communities and actual diffex
  forestGo=do.call(rbind,lapply(names(forestVals),function(x){
    print(x)
    res<-read.csv2(forestVals[[x]],sep=',',header=T,stringsAsFactors=F)
    res$disease=rep(x,nrow(res))
    return(res)
  }))%>%
    subset(namespace='biological_process')%>%
    select(term,q,name,Patient,disease)%>%
    tidyr::pivot_longer(Patient,names_to='networkType',values_to='sample')
  
  return(rbind(commGo,forestGo))
}

compareGOtoDistance<-function(communityDistanceFile){
  
   #process the distances
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
  
  assigned.comm<-mod.tab%>%subset(net1_type=='forest')%>%subset(net2_type=='community')%>%group_by(net1,hyp1)%>%filter(distance==min(distance))
  
}


distanceRidgelines<-function(communityDistanceFile){
  library(ggplot2)
  library(ggridges)
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
    dplyr::rename(net2_name='net1_name',net1_name='net2_name')
  
  ##first do each disease
  res=lapply(unique(mod.tab$hyp1),function(disease){
    
    red.tab<-subset(mod.tab,hyp1==disease)%>%
      subset(hyp2==disease)%>%
      rename(Community='net2')%>%
      mutate(Similarity=1-distance)
    
  #  red.missed<-subset(red.tab,net2_name%in%setdiff(annotes$net2_name,annotes$net1_name))%>%
  #    dplyr::select(net1_name,net2_name,distance)%>%
  #    dplyr::rename(net2_name='net1_name',net1_name='net2_name')
    comm.tab<-subset(red.tab,net2_type=='community')
    mean.dist<-comm.tab%>%group_by(Community)%>%summarize(meanSim=mean(Similarity))%>%
      ungroup()%>%
      arrange(meanSim)
    comStat<-read.csv(paste0(disease,'_communityStats.csv'))%>%
        select(Community,Nodes)%>%
      mutate(Community=as.character(Community))%>%
      mutate(numProteins=as.numeric(Nodes))
    comm.tab<-comm.tab%>%left_join(comStat)
    p<-ggplot(comm.tab,
              aes(x=Similarity,y=Community,fill=numProteins))+
      scale_y_discrete(limits = mean.dist$Community)+
              geom_density_ridges()+ggtitle('Community Forest Overlap')+scale_fill_viridis_b()
    print(p)
    ggsave(paste0(disease,'_ridgelines.pdf'),height=11)
    
    })
  
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
    dplyr::rename(net2_name='net1_name',net1_name='net2_name')

  ##first do each disease
  res=lapply(unique(mod.tab$hyp1),function(disease){
    red.tab<-subset(mod.tab,hyp1==disease)%>%
      subset(hyp2==disease)
    
    red.missed<-subset(red.tab,net2_name%in%setdiff(annotes$net2_name,annotes$net1_name))%>%
      dplyr::select(net1_name,net2_name,distance)%>%
      dplyr::rename(net2_name='net1_name',net1_name='net2_name')
    
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
    dplyr::rename(sample='net2_name')
    
    res<-cmdscale(as.dist(as.mat))
    
    colnames(res)<-c('Dim1','Dim2')
    
    full.tab<-res%>%
      as.data.frame()%>%
      tibble::rownames_to_column("sample")%>%
      left_join(red.annote)%>%
      left_join(neighborhood)
    
    p <-ggplot(full.tab,aes(x=Dim1,y=Dim2,col=graph,shape=graph,size=neighborhood))+
      geom_point()+ggtitle(paste(disease,'Hyphae'))+theme_minimal()+scale_color_viridis_d()
    return(p)
  })
   
  ##now do all diseases combined
  cowplot::plot_grid(plotlist=res)
  ggsave('individualPlots.pdf', useDingbats=FALSE)
  
  ##then do combined
  as.mat<-mod.tab%>%
    dplyr::select(net1_name,net2_name,distance)%>%
    rbind(missed)%>%
    tidyr::pivot_wider(values_from=distance,names_from=net2_name,
                       values_fn=list(distance=mean),values_fill=list(distance=1.0))%>%
    tibble::column_to_rownames('net1_name')%>%
    as.matrix()
  
  neighborhood=mod.tab%>%
    group_by(net2_name)%>%
    mutate(noVals=(distance==1.0))%>%
    subset(noVals==FALSE)%>%
    summarize(neighborhood=n_distinct(net1_name))%>%
    dplyr::rename(sample='net2_name')

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
  ggplot(full.tab,aes(x=Dim1,y=Dim2,col=disease,shape=graph,size=neighborhood))+geom_point()+scale_color_viridis_d()
  ggsave('allPlotsTogether.pdf',useDingbats=FALSE)
}

#'for a given community, get the go enrichment terms, and how many of the
#'closest forests have those terms
#'@param tumorType
#'@param community
getCommunityAndClosestsForests<-function(tumorType='brca',comm='0'){
  fors<-subset(minDists,disease==tumorType)%>%
    subset(community==comm)
  
  comGo<-subset(goStats,disease==tumorType)%>%
    filter(networkType=='Community')
  
  forGo<-subset(goStats,disease==tumorType)%>%
    subset(networkType=='Patient')%>%
    mutate(isClosest=sample%in%fors$graph)
  
  hasCom=comGo%>%group_by(name)%>%subset(networkType=='Community')%>%summarize(numComm=n())
  hasForest=comGo%>%group_by(name)%>%subset(networkType=='Patient')%>%summarize(numForests=n())
}


getAlluvialAcrossHyphae<-function(nodefiles=list(brca='brca0.01_nodelist.csv',
                                       coad='coad0.01_nodelist.csv',
                                       luad='luad0.01_nodelist.csv',
                                       hcc='hcc0.01_nodelist.csv',
                                       gbm='gbm0.01_nodelist.csv',
                                       hnscc='hnscc0.01_nodelist.csv')){

    library(ggalluvial)
  all.comms<-do.call(rbind,lapply(names(nodefiles),function(x){
    read.csv(nodefiles[[x]])%>%dplyr::select(Node,Community)%>%mutate(Hypha=x)}))
  all.comms$Community<-as.factor(all.comms$Community)
    p<-all.comms%>%#%>%subset(Gene%in%changing$Gene)%>%
    ggplot(aes(x=Hypha,stratum=Community,alluvium=Node,fill=Community,label=Hypha))+
    geom_flow(stat='alluvium',lode.guidance='frontback')+
    geom_stratum()+
   # ggtitle(paste("Changing",prefix,"proteins in",cl))+
    theme_minimal()+
    viridis::scale_fill_viridis(3,discrete=T)
  
}

#' getPanCanDistanceStats
#' Evaluates the distances between individual trees and the samples in other hyphae
#' @param fname fie name of distace csv
getPanCanDistanceStats<-function(fname='panCancerDistances.csv'){
  library(pheatmap)
  library(ggplot2)
  tab<-read.csv(fname)%>%subset(net2_type=='community')%>%select(-X)
  
  dtab<-tab%>%subset(hyp2=='panCan')%>%subset(hyp1!='panCan')
  norms<-read.csv("data/PDC_biospecimen_manifest_07182020_151323.csv")$Aliquot.Submitter.ID
 dmat<-dtab%>%select(net2,distance,net1)%>%distinct()%>%
    tidyr::pivot_wider(values_from=distance,names_from=net2)%>%
   subset(!net1%in%norms)%>%
    tibble::column_to_rownames('net1')
 
  annotes<-dtab%>%select(net1,disease="hyp1")%>%distinct()%>%
    tibble::column_to_rownames('net1')
  annotes$disease<-as.factor(as.character(annotes$disease))
    pheatmap(dmat,annotation_row = annotes,clustering_distance_rows = 'correlation', 
           clustering_distance_cols = 'correlation',clustering_method = 'ward.D2',
           labels_row=rep("",nrow(dmat)),filename='panCancerDistanceClusters.pdf')
  
    
    ##lets check to see that the communities can predict cancer
    library(glmnet)
  res<-cv.glmnet(x=as.matrix(dmat),y=as.factor(annotes[rownames(dmat),1]),family='multinomial')
  res2<-glmnet(x=as.matrix(dmat),y=as.factor(annotes[rownames(dmat),1]),family='multinomial')
    
  comms<-lapply(res2$beta,function(x) {
    vals<-x[,which(res$cvm==min(res$cvm))]
    names(vals)[vals!=0]})
  
  ##plot confusion matrix
  con<-confusion.glmnet(res,as.matrix(dmat),annotes[rownames(dmat),1],family='multinomial')
  pheatmap(con,cluster_rows = F,cluster_cols = F,filename='confusionMatrix.pdf')
  
  #lastly plot distances based on distance to clsuteres
  dists<-cmdscale(dist(dmat))
  colnames(dists)<-c('Dim1',"Dim2")
  ddf<-data.frame(dists,disease=annotes[rownames(dists),'disease'])
  p<-ggplot(ddf,aes(x=Dim1,y=Dim2,color=disease))+geom_point()
  ggsave('allPatientDistanceCluster.pdf',p)
}
#goStats<<-getGoValues()
getPanCanDistanceStats()
minDists<<-assignClosestNodes('panCancerDistances.csv')
res=distanceRidgelines('panCancerDistances.csv')
#forStats() ##how well do the forests recapitulate biology?
#plotNetworkDistances('panCancerDistances.csv') #how well do the communities summarize the diversity of the forests?
#forestGOStats()##we get more function for patients, and it agrees with other stuff

