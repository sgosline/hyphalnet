#!/usr/bin/env Rscript

#' Command line R script designed to create figures for
#' HyphalNet analysis
#' @author Sara JC Gosline
#'
#'
#'

if(!require("optparse"))
    install.packages('optparse')
if(!require(ggplot2))
    install.packages("ggplot2")
if(!require(dplyr))
    install.packages("dplyr")
if(!require('ggridges'))
    install.packages('ggridges')
if(!require(cowplot))
   install.packages('cowplot')
if(!require('stringr'))
    install.packages("stringr")
if(!require(tidyr))
    install.packages('tidyr')
if(!require(reticulate))
  install.packages("reticulate")

option_list = list(
    make_option(c("-d", "--distFile"), type="character", default=NULL,
                help="Distance file name", metavar="character"),
    make_option(c("-c", "--commFile"), type='character', default=NULL,
                help='Community stats file name', metavar='character'),
    make_option(c("s", "--synapseProj"), type='character', default=NULL,
                help="Synapse id of project to store table in", metavar='character'))
   # make_option(c("-o", "--out"), type="character", default="out.txt",
   #             help="output file name [default= %default]", metavar="character")
   # ;

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


#' Build basic predictor of multiple hyphae
#' @param communityDistanceFile
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

#' distanceRidglines - calculates distances to each community
#' @param communityDistanceFile
#'
distanceRidgelines<-function(communityDistanceFile, commStats){

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
    comStat<-read.csv(commStats)%>%
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

#' plotNetworkDistances
#' plots network distances to each other and to the communities
#' @param communityDistanceFile
plotNetworkDistances<-function(communityDistanceFile){

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
    if(require(cowplot)){
        cowplot::plot_grid(plotlist=res)
        ggsave('individualPlots.pdf', useDingbats=FALSE)
}
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

#' storeDistances
#' @param distFile filepath
#' @param synapseProj id of project
storeDistances<-function(distFile,synapseProj){
  require(dplyr)
  require(reticulate)
  syn=reticulate::import('synapseclient')
  sync=syn$login()
  tab<-read.csv(distFile,header=T,stringsAsFactors=F)
  
  mod.tab<-tab%>%rowwise()%>%
    mutate(net1_name=stringr::str_c(hyp1,net1,sep='_'))%>%
    mutate(net2_name=stringr::str_c(hyp2,net2,sep='_'))
  
  new.tab<-mod.tab%>%select(net1,net1_type,net2,net2_type,distance)
  #table <- synapser::synBuildTable("HyphalNetwork Distances", synapseProj, mod.tab)
  write.table(new.tab,file='tmp.csv',sep=',',row.names=FALSE,quote=FALSE)
  tab<-syn$build_table('HyphalNetwork Distances',synapseProj,'tmp.csv')
  sync$store(tab)  
  
}
                                        #goStats<<-getGoValues()

minDists<<-assignClosestNodes(opt$distFile)
res=distanceRidgelines(opt$distFile, opt$commFile)
#forStats() ##how well do the forests recapitulate biology?
plotNetworkDistances(opt$distFile) #how well do the communities summarize the diversity of the forests?

if(!is.null(opt$synapseProj))
  storeDistances(opt$distFile,opt$synapseProj)
