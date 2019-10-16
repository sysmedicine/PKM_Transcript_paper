DEGs_DEseq_table_overlap_quantile_th<-function(deg_DEseq_1,deg_DEseq_2,quantile_th){
  #function: extract DEGs based on quantile value and DEGs overlapping analysis
  #DEG_DEseq_1 is the output txt file of DEseq2
  #eg: deg_DEGseq_1<-as.matrix(read.csv("DESeq_output_fdr_005_NoRemBatch.txt",header=T,sep="\t"))
  #fdr<-0.05
  gene_s1<-as.matrix(deg_DEseq_1[,1])
  gene_s2<-as.matrix(deg_DEseq_2[,1])
  int_gene<-as.matrix(intersect(gene_s1,gene_s2))
  bg<-dim(int_gene)[1]
  ########################################################################
  #first cancer
  cutoff_1 = quantile(as.numeric(deg_DEseq_1[,"padj"]), quantile_th, na.rm = T)
  index_1<-which(as.numeric(deg_DEseq_1[,"padj"])<cutoff_1)
  deg_DEseq_1<-deg_DEseq_1[index_1,]
  
  index_up_1=which(as.numeric(deg_DEseq_1[,col.names="log2FoldChange"])>0)
  index_down_1=which(as.numeric(deg_DEseq_1[,col.names="log2FoldChange"])<0)
  
  deg_up_1<-as.matrix(deg_DEseq_1[index_up_1,1])
  deg_down_1<-as.matrix(deg_DEseq_1[index_down_1,1])
  deg_all_1<-rbind(deg_up_1,deg_down_1)
  num_deg_1<-dim(deg_all_1)[1]
  
  deg_up_11<-as.matrix(intersect(deg_up_1,int_gene))
  deg_down_11<-as.matrix(intersect(deg_down_1,int_gene))
  deg_all_11<-rbind(deg_up_11,deg_down_11)
  num_deg_11<-dim(deg_all_11)[1]
  #####################################################################
  #second cancer
  cutoff_2 = quantile(as.numeric(deg_DEseq_2[,"padj"]), quantile_th, na.rm = T)
  index_2<-which(as.numeric(deg_DEseq_2[,"padj"])<cutoff_2)
  deg_DEseq_2<-deg_DEseq_2[index_2,]
  
  index_up_2=which(as.numeric(deg_DEseq_2[,col.names="log2FoldChange"])>0)
  index_down_2=which(as.numeric(deg_DEseq_2[,col.names="log2FoldChange"])<0)
  deg_up_2<-as.matrix(deg_DEseq_2[index_up_2,1])
  deg_down_2<-as.matrix(deg_DEseq_2[index_down_2,1])
  deg_all_2<-rbind(deg_up_2,deg_down_2)
  num_deg_2<-dim(deg_all_2)[1]
  
  deg_up_22<-as.matrix(intersect(deg_up_2,int_gene))
  deg_down_22<-as.matrix(intersect(deg_down_2,int_gene))
  deg_all_22<-rbind(deg_up_22,deg_down_22)
  num_deg_22<-dim(deg_all_22)[1]
  #####################################################################
  deg_inform_yuan<-rbind(c(num_deg_1,dim(deg_up_1)[1],dim(deg_down_1)[1]),c(num_deg_2,dim(deg_up_2)[1],dim(deg_down_2)[1]))
  #colnames(deg_inform_yuan)<-c("all","up","down")
  #rownames(deg_inform_yuan)<-c("deg_1_yuan","deg_2_yuan")
  
  deg_inform_common<-rbind(c(num_deg_11,dim(deg_up_11)[1],dim(deg_down_11)[1]),c(num_deg_22,dim(deg_up_22)[1],dim(deg_down_22)[1]))
  #colnames(deg_inform_common)<-c("all","up","down")
  #rownames(deg_inform_common)<-c("deg_1_co","deg_2_co")
  deg_inform<-rbind(deg_inform_yuan,deg_inform_common)
  colnames(deg_inform)<-c("all","up","down")
  rownames(deg_inform)<-c("deg_1_yuan","deg_2_yuan","deg_1_common","deg_2_common")
  #####################################################################
  over_all<-as.matrix(intersect(deg_all_11,deg_all_22))
  over_up<-as.matrix(intersect(deg_up_11,deg_up_22))
  over_down<-as.matrix(intersect(deg_down_11,deg_down_22))
  
  union_unique_deg_over<-as.matrix(unique(rbind(deg_all_11,deg_all_22)))
  union_deg_num<-dim(union_unique_deg_over)[1]
  
  #p_binom=1-pbinom(dim(over_up)[1]+dim(over_down)[1]-1,dim(over_all)[1],0.5)
  
  p_hyper=1-phyper(dim(over_up)[1]+dim(over_down)[1]-1,num_deg_11,bg-num_deg_11,num_deg_22)
  
  over_result_yuan<-t(as.matrix(c(num_deg_1,num_deg_2,union_deg_num,dim(over_all)[1],dim(over_up)[1]+dim(over_down)[1],
                                  (dim(over_up)[1]+dim(over_down)[1])/dim(over_all)[1],
                                  (dim(over_up)[1]+dim(over_down)[1])/union_deg_num,
                                  dim(over_up)[1],dim(over_down)[1],p_hyper)))
  #colnames(over_result_yuan)<-c("deg_1_yuan","deg_2_yuan","overlap","consistent","percentage","co_up","co_down","p_binom","p_hyper")
  
  over_result_co<-t(as.matrix(c(num_deg_11,num_deg_22,union_deg_num,dim(over_all)[1],dim(over_up)[1]+dim(over_down)[1],
                                (dim(over_up)[1]+dim(over_down)[1])/dim(over_all)[1],
                                (dim(over_up)[1]+dim(over_down)[1])/union_deg_num,
                                dim(over_up)[1],dim(over_down)[1],p_hyper)))
  
  #colnames(over_result_co)<-c("deg_1_co","deg_2_co","overlap","consistent","percentage","co_up","co_down","p_binom","p_hyper")
  over_result<-rbind(over_result_yuan,over_result_co)
  colnames(over_result)<-c("deg_1","deg_2","union_deg_num","overlap","consistent","ratio_1","ratio_2","co_up","co_down","p_hyper")
  rownames(over_result)<-c("over_result_yuan","over_result_co")
  
  output<-list(deg_inform=deg_inform,
               over_result=over_result,
               over_up=over_up,
               over_down=over_down)
  return(output)
}