
library(dplyr)
library(pheatmap)
source("R/pdcFileProcessing.R")
#' Old plot using clusterProfiler
#' 
plotOldGSEA<-function(genes.with.values,prefix='',doPlot=FALSE){
 # print(head(genes.with.values))
  require(org.Hs.eg.db)
  mapping<-as.data.frame(org.Hs.egALIAS2EG)%>%
    dplyr::rename(Gene='alias_symbol')
  
  genes.with.values<-genes.with.values%>%
    dplyr::left_join(mapping,by='Gene')%>%
    arrange(desc(value))
  
  genelist=genes.with.values$value
  names(genelist)=genes.with.values$Gene

  gr<-clusterProfiler::gseGO(genelist[!is.na(genelist)],ont="BP",keyType="SYMBOL",
                             OrgDb=org.Hs.eg.db,pAdjustMethod = 'BH')#,eps=1e-10)
  #gr<-clusterProfiler::gseKEGG(genelist[!is.na(genelist)],organism='hsa',keyType="kegg",
  #OrgDb=org.Hs.eg.db,
  #                           pAdjustMethod = 'BH')#,eps=1e-10)
  
  # if(nrow(as.data.frame(gr))==0){
  #    gr<-clusterProfiler::gseGO(genelist[!is.na(genelist)],ont="BP",keyType="SYMBOL",
  #                             OrgDb=org.Hs.eg.db,pAdjustMethod = 'BH',pvalueCutoff = 0.1)#,eps=1e-10)
  #  }
  
  if(doPlot==TRUE){
    enrichplot::ridgeplot(gr,showCategory = 50,fill='pvalue')+ggplot2::ggtitle(paste0("GO BP Terms for ",prefix))
    ggplot2::ggsave(paste0(prefix,'_KEGG.pdf'),width=10,height=10)
  }
  
  return(summary(gr))
}


checkGSEAtermEnrichment<-function(){
  disTab<-readInAllFiles()%>%dplyr::rename(value='logratio')
  
  ##let's get disease means
  disNormTerms = disTab%>%subset(isNorm)%>%
    group_by(disease,Gene)%>%
    summarize(normMean=mean(value))
  
  disPats<-disTab%>%subset(!isNorm)%>%
    left_join(disNormTerms,by=c('disease','Gene'))
  
  ##only enrichment for gbm!
  disPatTerms = disPats%>%
    mutate(logFC=value-normMean)%>%
    group_by(disease,Patient)%>%
    dplyr::select(Gene,logFC)%>%
    dplyr::rename(value='logFC')%>%
    group_modify(plotOldGSEA)%>%
    ungroup()
  
  write.table(disPatTerms,file = 'GOenrichmentPerPatientVsNormal.tsv',sep='\t',row.names=F,col.names=T)
  disCounts<-disPats%>%group_by(disease)%>%
    summarize(numPats=n_distinct(Patient))
  
  disPatCounts = disPatTerms%>%
    group_by(Description,disease)%>%summarize(numTerms=n())%>%
    left_join(disCounts)%>%
    mutate(fracSamps=numTerms/numPats)
  
  topTerms = disPatCounts%>%
    group_by(disease)%>%
    top_n(20,wt=fracSamps)%>%
    dplyr::select(Description)
  dmat = disPatCounts%>%
    subset(Description%in%topTerms$Description)%>%
    pivot_wider(-c(numTerms,numPats),values_from = fracSamps,names_from=disease,values_fill=list(fracSamps=0.0))%>%
    tibble::column_to_rownames('Description')
  
  pheatmap(dmat,cellheight=10,cellwidth=10,filename = 'top20GOtermsFoundInpatientVsNormal.pdf')

}

checkGSEAtermEnrichment()