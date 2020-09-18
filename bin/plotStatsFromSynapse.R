if(!require(ggplot2))
    install.packages("ggplot2")
if(!require(dplyr))
    install.packages("dplyr")
if(!require(reticulate))
    install.packages('reticulate')
if(!require(argparse))
    install.packages('argparse')
if(!require(glmnet))
    install.packages('glmnet')
if(!require(pheatmap))
    install.packages('pheatmap')
if(!require(ggridges))
    install.packages('ggridges')
if(!require(tidyr))
    install.packages('tidyr')
if(!require(stringr))
    install.packages('stringr')

visualizeCommunitySize<-function(fname){
  tab<-read.csv(fname)%>%
    group_by(Partition)%>%
    summarize(comSize=n_distinct(Node))
  print(tab)
  tab

}

assignClosestNodes<-function(tab){
  #process the distances
 # tab<-read.csv(communityDistanceFile,header=T,stringsAsFactors=F)
  minDists=tab%>%subset(net2_type=='community')%>%
      group_by(net1)%>%
      filter(distance==min(distance))%>%
    dplyr::select(graph='net1',disease='hyp1',community='net2',commDis='hyp2',distance)

  minDists%>%
      group_by(disease,commDis)%>%
      summarize(numNets=n_distinct(graph))%>%
      ggplot(aes(x=disease,y=numNets,fill=commDis))+
      geom_bar(stat='identity',position='dodge')+scale_fill_viridis_d()+
      theme(axis.text.x = element_text(angle = 90, vjust=0.5,hjust=1))
  write.csv(minDists,file='networksAndClosestCommunity.csv')
  ggsave('networksAssignedToClosestCommunity.pdf')
  return(minDists)
}


distanceRidgelines<-function(tab){

#  tab<-read.csv(communityDistanceFile,header=T,stringsAsFactors=F)

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


plotNetworkDistances<-function(tab){
  library(cowplot)
  library(ggplot2)

#  tab<-read.csv(communityDistanceFile,header=T,stringsAsFactors=F)

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


#' getPanCanDistanceStats
#' Evaluates the distances between individual trees and the samples in other hyphae
#' @param fname fie name of distace csv
getPanCanDistanceStats<-function(tab,baseline){
                                        # tab<-read.csv(fname)%>%subset(net2_type=='community')%>%select(-X)
    tab<-tab%>%subset(net2_type=='community')
    dtab<-tab%>%subset(hyp2==baseline)%>%subset(hyp1!=baseline)
    print(head(dtab))
                                        # norms<-read.csv("data/PDC_biospecimen_manifest_07182020_151323.csv")$Aliquot.Submitter.ID
    dmat<-dtab%>%select(net2,distance,net1)%>%distinct()%>%
        tidyr::pivot_wider(values_from=distance,names_from=net2)%>%
                                        # subset(!net1%in%norms)%>%
        tibble::column_to_rownames('net1')
                                        #  print(dmat)
    rownames(dtab)<-c()
    annotes<-dtab%>%select(net1,disease="hyp1")%>%
        distinct()%>%
        tibble::column_to_rownames('net1')

    ##lets check to see that the communities can predict cancer
    res<-cv.glmnet(x=as.matrix(dmat),y=as.factor(annotes[rownames(dmat),1]),family='multinomial')
    res2<-glmnet(x=as.matrix(dmat),y=as.factor(annotes[rownames(dmat),1]),family='multinomial')

    comms<-lapply(res2$beta,function(x) {
        vals<-x[,which(res$cvm==min(res$cvm))]
        names(vals)[vals!=0]})


    annotes$disease<-as.factor(as.character(annotes$disease))
    pheatmap(dmat,annotation_row = annotes,clustering_method = 'ward.D2',
             labels_row=rep("",nrow(dmat)),filename=paste0(baseline,'distanceClusters.pdf'))

    pheatmap(dmat[,unique(unlist(comms))],annotation_row = annotes,clustering_method = 'ward.D2',
             labels_row=rep("",nrow(dmat)),filename=paste0(baseline,'distanceClustersSelectedComms.pdf'))



    ##plot confusion matrix
    con<-confusion.glmnet(res,as.matrix(dmat),annotes[rownames(dmat),1],family='multinomial')
    pheatmap(con,cluster_rows = F,cluster_cols = F,filename='confusionMatrix.pdf')

                                        #lastly plot distances based on distance to clsuteres
    dists<-cmdscale(dist(dmat))
    colnames(dists)<-c('Dim1',"Dim2")
    ddf<-data.frame(dists,disease=annotes[rownames(dists),'disease'])
    p<-ggplot(ddf,aes(x=Dim1,y=Dim2,color=disease))+geom_point()
    ggsave('allSamplsDistanceCluster.pdf',p)

    ##now get only selected proteins
    rdists<-cmdscale(dist(dmat[,unique(unlist(comms))]))
    colnames(rdists)<-c('Dim1',"Dim2")
    ddf<-data.frame(rdists,disease=annotes[rownames(dists),'disease'])
    p<-ggplot(ddf,aes(x=Dim1,y=Dim2,color=disease))+geom_point()
    ggsave('allSamplsDistanceClusterSelectedComms.pdf',p)


}
parser <- argparse::ArgumentParser()
parser$add_argument('-s','--synId',dest='synid',default="",help='synapse id of table to pull')
parser$add_argument('-b','--baseline',dest='baseline',default='panCan0.01',help='Community distance to compare to')
args <- parser$parse_args()

if(args$synid!=""){
  #condaenv="C:\\Users\\gosl241\\OneDrive - PNNL\\Documents\\GitHub\\amlresistancenetworks\\renv\\python\\r-reticulate\\"
  library(reticulate)
  #reticulate::use_condaenv(condaenv)
  syn=reticulate::import('synapseclient')
  sync=syn$login()

  fname=sync$tableQuery(paste("select * from",args$synid))$asDataFrame()


  getPanCanDistanceStats(fname,args$baseline)
  minDists<<-assignClosestNodes(fname)
 # res=distanceRidgelines(fname)
}
