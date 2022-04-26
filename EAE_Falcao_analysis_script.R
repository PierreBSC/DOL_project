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

list_files = list.files("Desktop/DOL_meta_analysis/Falcao/",pattern = "tab",full.names = T)
names(list_files) = list.files("Desktop/DOL_meta_analysis/Falcao/",pattern = "tab",full.names = F)

data_raw = read.delim(list_files[1],sep = "\t")
x = table(data_raw[,1])
genes_to_keep = names(which(x==1))
data_raw = data_raw[data_raw[,1]%in%genes_to_keep,]
rownames(data_raw) = data_raw[,1]
data_raw = as(as.matrix(data_raw[,-1]),"dgCMatrix")
condition = rep(names(list_files)[1],ncol(data_raw))

for (k in 2:length(list_files)){
  u = read.delim(list_files[k],sep = "\t")
  u = u[u[,1]%in%genes_to_keep,]
  rownames(u) = u[,1]
  u = as(as.matrix(u[,-1]),"dgCMatrix")
  data_raw = cbind(data_raw,u)
  condition = c(condition,rep(names(list_files)[k],ncol(u)))
}

#B)Filtering

Lib_size = colSums(data_raw)
hist(log10(Lib_size),100)
Gene_size = rowSums(data_raw)
hist(log10(Gene_size+1),100)

data_count = data_raw[Gene_size>10,Lib_size>10^5]
colnames(data_count) = 1:ncol(data_count)

#C)Metadata

Annotation = getGEO("GSE113973")
Annotation_samples = pData(Annotation$GSE113973_series_matrix.txt.gz)
list_plates = Annotation_samples$title
list_plates = strsplit(list_plates,split = "_")
list_plates = unlist(lapply(list_plates,FUN = function(x) {paste(x[1],x[2],x[3],sep = "_")}))
Condition_samples = data.frame(Sample = list_plates,Condition = Annotation_samples$description )
Condition_samples = unique(Condition_samples)

Condition_cells = apply(Condition_samples,MARGIN = 2,FUN = function(x) {rep(x,each=384)})
Condition_cells = as.data.frame(Condition_cells)

Condition_cell_counts = Condition_cells[Lib_size>10^5,]

#II) Data analysis

#A)Analysis per se

r <- Pagoda2$new(data_count,log.scale=T)
r$adjustVariance(plot=T,gam.k=10)
r$calculatePcaReduction(nPcs=50,n.odgenes = 3000)
r$makeKnnGraph(k=30,type='PCA',distance = "cosine")
r$getKnnClusters(method=multilevel.community,type='PCA')

r$getDifferentialGenes(type = "PCA",clusterType = "community",verbose = T,z.threshold = 3)

umap_plot = umap(r$reductions$PCA,n_neighbors = 30,spread = 4,
                 n_components = 2,metric = "cosine",verbose = T)

plot(umap_plot,xlab="UMAP 1",ylab="UMAP 2",pch=21,bg=string.to.colors(r$clusters$PCA$community))

#B)Meta clustering

Mean_expression_cluster = aggregate(as.matrix(r$counts[,r$getOdGenes(3000)]),FUN = mean,by=list(r$clusters$PCA$community))
rownames(Mean_expression_cluster) =Mean_expression_cluster$Group.1
Mean_expression_cluster = t(Mean_expression_cluster[,-1])
Meta_clustering = pheatmap(cor(Mean_expression_cluster,method = "spearman"),clustering_method = "ward")
Meta_clustering = Meta_clustering$tree_col
Order_cluster = Meta_clustering$order

boxplot(r$counts[,"Mog"]~factor(r$clusters$PCA$community,Order_cluster))
plot(umap_plot,xlab="UMAP 1",ylab="UMAP 2",pch=21,bg=string.to.colors(r$clusters$PCA$community==5))
plot(umap_plot,xlab="UMAP 1",ylab="UMAP 2",pch=21,bg=color_convertion(r$counts[,"eGFP"]))

#C)Isolation of oligodendrocytes and DE analysis

Oligo_cluster = c(7,4,9,3,1,10,5)

data_oligo = data_count[,r$clusters$PCA$community%in%Oligo_cluster]
Condition_oligo = Condition_cell_counts$Condition[r$clusters$PCA$community%in%Oligo_cluster]

TPM_oligo = log2(1+t(t(data_oligo)/colSums(data_oligo)*10^6))
Mean_log2FC_oligo = aggregate(as.matrix(t(TPM_oligo)),FUN = mean ,by=list(Condition_oligo))
rownames(Mean_log2FC_oligo) = Mean_log2FC_oligo[,1]
Mean_log2FC_oligo = t(Mean_log2FC_oligo[,-1])
Log2FC_oligo = loess(Mean_log2FC_oligo[,2]~Mean_log2FC_oligo[,1],degree = 2)
Log2FC_oligo = Log2FC_oligo$residuals


##Here we use DESeq2 and not the cloglog regression model used in for the other datasets as it is a Smart-seq2 dataset (not UMI based)
library(DESeq2)

des = DESeqDataSetFromMatrix(countData = data_oligo,
                             colData = Condition_cell_counts[r$clusters$PCA$community%in%Oligo_cluster,],
                             design = ~Condition)
des$Condition = relevel(des$Condition,ref = "Complete Freund’s adjuvant")

des = DESeq(des,fitType ='local')
plotDispEsts(des)
coef(des)
VST_oligo = varianceStabilizingTransformation(des)
shrinked_log2FC = lfcShrink(dds = des,coef = "Condition_MOG35.55.in.complete.Freund.s.adjuvant_vs_Complete.Freund.s.adjuvant",type = "normal")
shrinked_log2FC = as.data.frame(shrinked_log2FC)
shrinked_log2FC$padj[shrinked_log2FC$padj==0]=10^-300

DOL_genes = read.delim("Documents/DOL_genes.txt")
DOL_genes = DOL_genes$Genes
pdf("Desktop/DOL_meta_analysis/Falcao/Volcanoplot.pdf",width = 7,height = 7,useDingbats = F)
par(las=1,bty="l")
plot(Log2FC_oligo,-log10(shrinked_log2FC[names(Log2FC_oligo),6]),pch=21,xlim=c(-10,10),
     ylim=c(0,350),yaxs="i",xaxs="i",xlab="Log2FC (EAE vs CFA)",ylab="-Log10(FDR)",cex.lab=1.3,cex.axis=1.3,
     bg=string.to.colors(names(Log2FC_oligo)%in%DOL_genes,colors = c("grey70","orange")),cex=1.5)
plot(Log2FC_oligo,-log10(shrinked_log2FC[names(Log2FC_oligo),6]),pch=21,xlim=c(-10,10),
     ylim=c(0,350),yaxs="i",xaxs="i",xlab="Log2FC (EAE vs CFA)",ylab="-Log10(FDR)",cex.lab=1.3,cex.axis=1.3,
     bg=string.to.colors(names(Log2FC_oligo)%in%DOL_genes,colors = c("grey70","orange")),cex=1.5)
text(Log2FC_oligo,-log10(shrinked_log2FC[names(Log2FC_oligo),6]),labels = names(Log2FC_oligo))
dev.off()


#D)GSEA analysis


Enrichment_DOL_signature = bulk.gsea(Log2FC_oligo,set.list = list("DOL"=DOL_genes))
par(las=1)
gsea(Delta_log2FC,DOL_genes)

