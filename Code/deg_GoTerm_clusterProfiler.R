deg_GoTerm_clusterProfiler<-function(clusterProfile_result){
  #function: creat a file for go term result
  
  path_result<-clusterProfile_result@result
  GO_ID<-as.matrix(path_result$ID)
  term_name<-as.matrix(path_result$Description)
  pvalue<-as.matrix(path_result$pvalue)
  p.adjust<-as.matrix(path_result$p.adjust)
  qvalue<-as.matrix(path_result$qvalue)
  #Count<-as.matrix(path_result$Count)
  geneID<-as.matrix(path_result$geneID)
  path_out_table<-cbind(GO_ID,term_name,pvalue,p.adjust,qvalue,geneID)
  colnames(path_out_table)<-c("GO_ID","term_name","pvalue","p.adjust","qvalue","geneID")
  return(path_out_table)
}
