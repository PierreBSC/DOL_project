library(pagoda2)
library(fifer)
library(uwot)
library(pheatmap)
library(GEOquery)

color_convertion=function(x,max_scale=NULL) {
  f <- colorRamp(c("white","yellow","orange","red"))
  x=as.numeric(x)
  if (is.null(max_scale)) {
    max_scale=quantile(x,0.999,na.rm = T)
  }
  x_prime=ifelse(x>max_scale,max_scale,x)
  x_prime=x_prime/max_scale
  x_color=f(x_prime)/255
  x_color[!complete.cases(x_color),]=c(0,0,0)
  x_color=rgb(x_color)
  return(x_color)
}


#I) Data loading and processing 

#A)Loading 

list_files = list.files("Desktop/DOL_meta_analysis/Wheeler/GSE129609_RAW///",full.names = T)
names(list_files) = list.files("Desktop/DOL_meta_analysis/Wheeler/GSE129609_RAW/",full.names = F)

data_raw = read.delim(list_files[1],sep = "\t",row.names = 1)
data_raw = as(as.matrix(data_raw),"dgCMatrix")
condition = rep(names(list_files)[1],ncol(data_raw))


for (k in 2:length(list_files)){
  print(k)
  u = read.delim(list_files[k],sep = "\t",row.names = 1)
  u = as(as.matrix(u),"dgCMatrix")
  l = intersect(rownames(u),rownames(data_raw))
  data_raw = cbind(data_raw[l,],u[l,])
  condition = c(condition,rep(names(list_files)[k],ncol(u)))
}

condition_bis = strsplit((condition),split = "_",fixed = T)
condition_bis = unlist(lapply(condition_bis,FUN = function(x) {x[1]}))
#B)Filtering

Lib_size = colSums(data_raw)
hist(log10(Lib_size),100)
Gene_size = rowSums(data_raw)
hist(log10(Gene_size+1),100)

data_count = data_raw[Gene_size>100,Lib_size>500]
colnames(data_count) = 1:ncol(data_count)

#C)Metadata

Annotation = getGEO("GSE129609")
Annotation_samples = pData(Annotation$GSE129609_series_matrix.txt.gz)
List_sample =  strsplit(names(list_files),split = "_",fixed = T)
List_sample = unlist(lapply(List_sample,FUN = function(x) {x[1]}))

Type_sample = Annotation_samples[List_sample,1]
Type_sample[grepl(Type_sample,pattern = "Priming")] = "Priming"
Type_sample[grepl(Type_sample,pattern = "CFA")] = "CFA"
Type_sample[grepl(Type_sample,pattern = "Acute")] = "Acute"
Type_sample[grepl(Type_sample,pattern = "Naive")] = "Naive"
Type_sample[grepl(Type_sample,pattern = "Remitting")] = "Remitting"
names(Type_sample) = List_sample

Cell_condition = Type_sample[condition_bis]

Cell_condition_count = Cell_condition[Lib_size>500]

#D)Gene selection

Zero_proporion = Matrix::rowSums((data_count)==0) / ncol(data_count)
Total_expression = Matrix::rowSums(data_count)
plot(log10(Total_expression),Zero_proporion[names(Total_expression)])
Regression_zero = loess(Zero_proporion~log10(Total_expression),subset = Total_expression>0)
Regression_zero = Regression_zero$residuals
Regression_zero = Regression_zero[order(Regression_zero,decreasing = T)]
Selected_genes = names(Regression_zero[1:1000])
plot(Regression_zero[1:1000],type="l")


colnames(data_count) = 1:ncol(data_count)
r <- Pagoda2$new(data_count,log.scale=F)
r$adjustVariance(plot=T,gam.k=10)
r$calculatePcaReduction(nPcs=50,odgenes = Selected_genes )
r$makeKnnGraph(k=30,type='PCA',distance = "cosine")
r$getKnnClusters(method=multilevel.community,type='PCA')

r$getDifferentialGenes(type = "PCA",clusterType = "community",verbose = T,z.threshold = 3)

umap_plot = umap(r$reductions$PCA,n_neighbors = 30,spread = 4,
                 n_components = 2,metric = "cosine",verbose = T)

plot(umap_plot,xlab="UMAP 1",ylab="UMAP 2",pch=21,bg=string.to.colors(r$clusters$PCA$community))

#B)Meta clustering

Mean_expression_cluster = aggregate(as.matrix(r$counts[,r$getOdGenes(1000)]),FUN = mean,by=list(r$clusters$PCA$community))
rownames(Mean_expression_cluster) =Mean_expression_cluster$Group.1
Mean_expression_cluster = t(Mean_expression_cluster[,-1])
Meta_clustering = pheatmap(cor(Mean_expression_cluster,method = "spearman"),clustering_method = "ward")
Meta_clustering = Meta_clustering$tree_col
Order_cluster = Meta_clustering$order
boxplot(r$counts[,"Mog"]~factor(r$clusters$PCA$community,levels =Order_cluster ),outline=F)

Oligo_cluster = 18
plot(umap_plot,xlab="UMAP 1",ylab="UMAP 2",pch=21,bg=string.to.colors(r$clusters$PCA$community==Oligo_cluster))

#III)DE analaysis 

#A)Log2FC computation

DOL_genes = read.delim("Documents/DOL_genes.txt")
DOL_genes = DOL_genes$Genes


data_oligo = data_count[,r$clusters$PCA$community%in%Oligo_cluster]
Cell_condition_oligo = Cell_condition_count[r$clusters$PCA$community%in%Oligo_cluster]
Cell_condition_oligo = factor(Cell_condition_oligo,levels = c("Naive","CFA","Priming","Acute","Remitting"))

TPM_oligo = log2(1+t(t(data_oligo)/colSums(data_oligo)*10^6))
Mean_log2FC_oligo = aggregate(as.matrix(t(TPM_oligo)),FUN = mean ,by=list(Cell_condition_oligo))
rownames(Mean_log2FC_oligo) = Mean_log2FC_oligo[,1]
Mean_log2FC_oligo = t(Mean_log2FC_oligo[,-1])

plot(Mean_log2FC_oligo[,"Acute"]~Mean_log2FC_oligo[,"CFA"])
Delta_log2FC = loess(Mean_log2FC_oligo[,"Acute"]~Mean_log2FC_oligo[,"CFA"],degree = 1)
Delta_log2FC = Delta_log2FC$residuals


#B)Barplot of 'intersting' genes

barplot_gene_EAE = function(gene="Serpina3n") {
  x = Mean_log2FC_oligo[gene,]
  par(las=1,bty="l")
  barplot(x,col=c("grey70","grey80","orange","red3","darkred"),ylim=c(0,1.3*max(x)),
          ylab="Mean expression (Log2(1+TPM))",xlab="Condition",main=gene,
          cex.lab=1.3,cex.axis=1.15,cex.names = 1.15,cex.main = 2)
}

pdf("Desktop/DOL_meta_analysis/Wheeler/Barplot_genes.pdf",width = 8,height = 8,useDingbats = F)
for (k in intersect(DOL_genes,names(Delta_log2FC))) {
  barplot_gene_EAE(k)
}
dev.off()

#C)P-value computation for all genes for the Peak condition

Lib_size_count = Lib_size[Lib_size>500]
Lib_size_oligo = Lib_size_count[r$clusters$PCA$community%in%Oligo_cluster]
Lib_size_oligo_bis = Lib_size_oligo[Cell_condition_oligo%in%c("CFA","Acute")]
Data_oligo_bis = data_oligo[,Cell_condition_oligo%in%c("CFA","Acute")]
Cell_condition_bis = Cell_condition_oligo[Cell_condition_oligo%in%c("CFA","Acute")]

cloglog_regression = function(x,l,condition) {
  x_bin = x>0
  model_temp = glm(x_bin~1+offset(log(l))+condition,family = binomial(link="cloglog"))
  model_temp = summary(model_temp)
  return(model_temp$coefficients[,4])
}

Cloglog_result = apply(as.matrix(Data_oligo_bis),MARGIN = 1,FUN = function(x) {cloglog_regression(x,Lib_size_oligo_bis,Cell_condition_bis)})
Cloglog_result = t(Cloglog_result)
Cloglog_result[,2] = p.adjust(Cloglog_result[,2],method = "fdr")

P_value_cloglog = Cloglog_result[,2]
plot(Delta_log2FC,-log10(P_value_cloglog))


pdf("Desktop/DOL_meta_analysis/Wheeler//Volcanoplot.pdf",width = 7,height = 7,useDingbats = F)
par(las=1,bty="l")
plot(Delta_log2FC,-log10(P_value_cloglog),pch=21,xlim=c(-8,8),ylim=c(0,15),
     yaxs="i",xaxs="i",xlab="Log2FC (Acute vs CFA)",ylab="-Log10(FDR)",cex.lab=1.3,cex.axis=1.3,
     bg=string.to.colors(names(Delta_log2FC)%in%DOL_genes,colors = c("grey70","orange")),cex=1.5)
plot(Delta_log2FC,-log10(P_value_cloglog),pch=21,xlim=c(-8,8),ylim=c(0,15),
     yaxs="i",xaxs="i",xlab="Log2FC (Acute vs CFA)",ylab="-Log10(FDR)",cex.lab=1.3,cex.axis=1.3,
     bg=string.to.colors(names(Delta_log2FC)%in%DOL_genes,colors = c("grey70","orange")),cex=1.5)
text(Delta_log2FC,-log10(P_value_cloglog),labels = names(Delta_log2FC))
dev.off()



##D)P-value for all conditions but only for 6 genes 

List_gene = c("B2m","H2-D1","H2-K1","Serpina3n","C4b","Lbh")
List_condition = c("Priming","Acute","Remitting")
Table_p_value = matrix(0,nrow = length(List_gene),ncol=length(List_condition))
rownames(Table_p_value) = List_gene
colnames(Table_p_value) = List_condition


for (i in List_gene) {
  for (j in List_condition) {
    Lib_size_oligo_bis = Lib_size_oligo[Cell_condition_oligo%in%c("CFA",j)]
    Data_oligo_bis = data_oligo[i,Cell_condition_oligo%in%c("CFA",j)]
    Cell_condition_bis = Cell_condition_oligo[Cell_condition_oligo%in%c("CFA",j)]
    x = cloglog_regression(Data_oligo_bis,Lib_size_oligo_bis,Cell_condition_bis)
    Table_p_value[i,j]=x[2]
  }
}

Table_p_value = matrix(p.adjust(Table_p_value,method = "fdr"),nrow = length(List_gene),ncol=length(List_condition))
rownames(Table_p_value) = List_gene
colnames(Table_p_value) = List_condition

#E)GSEA analysis


Enrichment_DOL_signature = bulk.gsea(Delta_log2FC,set.list = list("DOL"=DOL_genes))
par(las=1)
gsea(Delta_log2FC,DOL_genes)





