library(pagoda2)
library(uwot)
library(Matrix)
library(MASS)
library(igraph)
library(pheatmap)
library(liger)



#I)Auxiliary functions

string.to.colors = function (string, colors = NULL) 
{
  if (is.factor(string)) {
    string = as.character(string)
  }
  if (!is.null(colors)) {
    if (length(colors) != length(unique(string))) {
      (break)("The number of colors must be equal to the number of unique elements.")
    }
    else {
      conv = cbind(unique(string), colors)
    }
  }
  else {
    conv = cbind(unique(string), rainbow(length(unique(string))))
  }
  unlist(lapply(string, FUN = function(x) {
    conv[which(conv[, 1] == x), 2]
  }))
}



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


data_1 =read.10x.matrices("Desktop/Aging_mouse_SVZ/Cellranger_aging_data/O1/filtered_gene_bc_matrices//mm10/",version = "V2")
data_2 =read.10x.matrices("Desktop/Aging_mouse_SVZ/Cellranger_aging_data/O2/filtered_gene_bc_matrices/mm10/",version = "V2")
data_3 =read.10x.matrices("Desktop/Aging_mouse_SVZ/Cellranger_aging_data/O3/filtered_gene_bc_matrices/mm10/",version = "V2")
data_4 =read.10x.matrices("Desktop/Aging_mouse_SVZ/Cellranger_aging_data/Y1//filtered_gene_bc_matrices/mm10/",version = "V2")
data_5 =read.10x.matrices("Desktop/Aging_mouse_SVZ/Cellranger_aging_data/Y2/filtered_gene_bc_matrices/mm10/",version = "V2")
data_6 =read.10x.matrices("Desktop/Aging_mouse_SVZ/Cellranger_aging_data/Y3/filtered_gene_bc_matrices/mm10/",version = "V2")


Plate = c(rep("Old_1",ncol(data_1)),rep("Old_2",ncol(data_2)),rep("Old_3",ncol(data_3)),
          rep("Young_1",ncol(data_4)),rep("Young_2",ncol(data_5)),rep("Young_3",ncol(data_6)))

Merged_data_set = cbind(data_1,data_2,data_3,data_4,data_5,data_6)
colnames(Merged_data_set) = paste("Cell_",1:ncol(Merged_data_set),sep="")

Condition = substr(Plate,1,3)
names(Condition) = colnames(Merged_data_set)
names(Plate) = colnames(Merged_data_set)
rm(data_1,data_2,data_3,data_4,data_5,data_6)

lib_size = colSums(Merged_data_set)
gene_size = rowSums(Merged_data_set)

par(las=1)
hist(log10(lib_size),100,xlab="Cell Lib size (Log10)")
abline(v=log10(1000),lwd=2,lty=2,)

hist(log10(gene_size),100,xlab="Gene size (Log10)")
abline(v=log10(100),lwd=2,lty=2,)

data_count = Merged_data_set[gene_size>10,lib_size>1000]
data_count = data_count[rownames(data_count)!="Smim20",]

Plate_count = Plate[lib_size>1000]
Condition_count = Condition[lib_size>1000]

Zero_proportion = rowSums(data_count==0)/ncol(data_count)
Mean_expression = rowMeans(data_count)

par(las=1)
plot(log10(Mean_expression),Zero_proportion,pch=21,bg="orange",xaxs='i',yaxs='i')
Zeros_excess_proportion = loess(Zero_proportion~log10(Mean_expression),degree = 2)
Zeros_excess_proportion = Zeros_excess_proportion$residuals
Zeros_excess_proportion = Zeros_excess_proportion[order(Zeros_excess_proportion,decreasing = T)]
Selected_genes = names(Zeros_excess_proportion[1:2000])


#II)Analysis per se

x = rownames(data_count)
x = table(x)
data_count = data_count[names(which(x==1)),]
Selected_genes = intersect(Selected_genes,rownames(data_count))


r <- Pagoda2$new(data_count,log.scale=F)
Selected_genes = intersect(Selected_genes,colnames(r$counts))
r$adjustVariance(plot=T,gam.k=10)
r$calculatePcaReduction(nPcs=60,n.odgenes=1e3,odgenes = Selected_genes)
r$makeKnnGraph(k=30,type='PCA',distance = "cosine")
r$getKnnClusters(method=multilevel.community,type='PCA')
r$getDifferentialGenes(type = "PCA",clusterType = "community",verbose = T,z.threshold = 3,upregulated.only = T)

umap_plot = umap(r$reductions$PCA,n_neighbors = 30,spread = 4,
                 n_components = 2,metric = "cosine",verbose = T,)

plot(umap_plot,pch=21,bg=string.to.colors(r$clusters$PCA$community),
     xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2")
Mog_cluster_expression = aggregate(r$counts[,"Mog"],by=list(r$clusters$PCA$community),FUN=mean)


Mean_expression_cluster = aggregate(as.matrix(r$counts[,Selected_genes]),FUN=mean,by=list(r$clusters$PCA$community))
rownames(Mean_expression_cluster) = Mean_expression_cluster$Group.1
Mean_expression_cluster = t(Mean_expression_cluster[,-1])

Order_cluster =pheatmap(cor(Mean_expression_cluster,method = "spearman"),clustering_method = "ward.D")
Order_cluster = Order_cluster$tree_row$order
boxplot(r$counts[,"Plp1"]~factor(r$clusters$PCA$community,Order_cluster),outline=F)

Oligo_clusters = c(13,1,8)
Oligo_cells = colnames(data_count)[r$clusters$PCA$community%in%Oligo_clusters]
plot(umap_plot,pch=21,bg=string.to.colors(r$clusters$PCA$community%in%Oligo_clusters),
     xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2")


#III)Analysis of Oligo only

data_oligo = Merged_data_set[,Oligo_cells]

gene_size_data_oligo = rowSums(data_oligo)
hist(log10(gene_size_data_oligo),100)
abline(v=log10(10))

data_oligo = data_oligo[gene_size_data_oligo>50,]

Zero_proportion = rowSums(data_oligo==0)/ncol(data_oligo)
Mean_expression = rowMeans(data_oligo)

par(las=1)
plot(log10(Mean_expression),Zero_proportion,pch=21,bg="orange",xaxs='i',yaxs='i')
Zeros_excess_proportion = loess(Zero_proportion~log10(Mean_expression),degree = 2)
Zeros_excess_proportion = Zeros_excess_proportion$residuals
Zeros_excess_proportion = Zeros_excess_proportion[order(Zeros_excess_proportion,decreasing = T)]
Selected_genes = names(Zeros_excess_proportion[1:300])

Condition_oligos  = Condition[Oligo_cells]
Plate_oligos = Plate[Oligo_cells]
  
data_oligo = data_oligo[rownames(data_oligo)!="Smim20",]
r_oligo <- Pagoda2$new(data_oligo,log.scale=F)
r_oligo$adjustVariance(plot=T,gam.k=10)
r_oligo$calculatePcaReduction(nPcs=20,n.odgenes=1e3,odgenes = Selected_genes)
r_oligo$makeKnnGraph(k=30,type='PCA',distance = "cosine")
r_oligo$getKnnClusters(method=multilevel.community,type='PCA')
r_oligo$getDifferentialGenes(type = "PCA",clusterType = "community",verbose = T,z.threshold = 3,upregulated.only = T)


umap_plot_oligo = umap(r_oligo$reductions$PCA,n_neighbors = 30,
                           n_components = 2,metric = "cosine",verbose = T)
plot(umap_plot_oligo,pch=21,bg=string.to.colors(Condition_oligos),
     xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2")
plot(umap_plot_oligo,pch=21,bg=string.to.colors(r_oligo$clusters$PCA$community==2),
     xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2")

plot(umap_plot_oligo[Condition_oligos=="You",],pch=21,bg="grey",cex=0.7,
     xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2",xlim=c(-3,4))
plot(umap_plot_oligo[Condition_oligos=="Old",],pch=21,bg="orange",cex=0.7,
     xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2",xlim=c(-3,4))

plot(r_oligo$reductions$PCA[,c()])
plot(umap_plot_oligo,pch=21,bg=color_convertion(r_oligo$counts[,"Serpina3n"]),
     xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2")

#III)Meta-analysis part

DOL_genes = rownames(read.delim("Desktop/Meta_analysis_DOL_signature/Data_Log2FC/5XFAD_DOL.txt"))
DOL_genes = intersect(DOL_genes,rownames(data_count))

#A)Computing Log2FC

TPM_data = data_count[,Oligo_cells]
TPM_data = log2(1+ TPM_data/colSums(TPM_data)*10^6)
Condition_oligos = factor(Condition_oligos,levels = c("You","Old"))

Mean_oligo_expression = aggregate(as.matrix(t(TPM_data)),FUN = mean,by=list(Condition_oligos))
rownames(Mean_oligo_expression) = Mean_oligo_expression[,1]
Mean_oligo_expression = t(Mean_oligo_expression[,-1])

pdf("Desktop/Aging_mouse_SVZ//Figures//Log2FC_Oligo.pdf",width = 8,height = 8)
par(las=1,bty="l")
plot(Mean_oligo_expression,pch=21,xaxs="i",yaxs="i",xlim=c(0,max(Mean_oligo_expression)*1.1),
     ylim=c(0,max(Mean_oligo_expression)*1.1),xlab="Young mice",ylab="Old mice",cex=1.5,cex.lab=1.4,cex.axis=1.4,
     bg=string.to.colors(rownames(Mean_oligo_expression)%in%DOL_genes,colors = c("grey80","firebrick2")))
abline(0,1,lwd=2,lty=2,col="black")
plot(Mean_oligo_expression,pch=21,xaxs="i",yaxs="i",xlim=c(0,max(Mean_oligo_expression)*1.1),
     ylim=c(0,max(Mean_oligo_expression)*1.1),
     bg=string.to.colors(rownames(Mean_oligo_expression)%in%DOL_genes,colors = c("grey","orange")))
text(Mean_oligo_expression[,1],Mean_oligo_expression[,2],labels = rownames(Mean_oligo_expression))
dev.off()

Log2_FC = lm(Mean_oligo_expression[,2]~Mean_oligo_expression[,1])
Log2_FC = Log2_FC$residuals



Enrichment_DOL_signature = bulk.gsea(Log2_FC,set.list = list("DOL"=DOL_genes))
pdf("/home/data/Pierre/Untitled/home/pbost//Desktop/Aging_mouse_SVZ/Figures//Enrichment_DOL.pdf",width = 6,height = 6) 
par(las=1)
gsea(Log2_FC,DOL_genes)
dev.off()

##C)Cloglog regression 


cloglog_regression = function(x,l,condition) {
  x_bin = x>0
  model_temp = glm(x_bin~1+offset(log(l))+condition,family = binomial(link="cloglog"))
  model_temp = summary(model_temp)
  return(model_temp$coefficients[,4])
}

lib_size_oligo = lib_size[Oligo_cells]
Cloglog_result = apply(as.matrix(data_oligo),MARGIN = 1,FUN = function(x) {cloglog_regression(x,lib_size_oligo,Condition_oligos)})
Cloglog_result = t(Cloglog_result)
Cloglog_result = as.data.frame(Cloglog_result)
x= p.adjust(Cloglog_result[,2],method = "fdr")
x[x==0] = 10^-262
Cloglog_result$fdr = x

DOL_genes = read.delim("Documents/DOL_genes.txt")
DOL_genes = DOL_genes$Genes

DE_table = data.frame(log2FC = Log2_FC[rownames(Cloglog_result)],
                      Log_fdr = -log10(Cloglog_result$fdr),row.names =rownames(Cloglog_result) )
DE_table = DE_table[!grepl(rownames(DE_table),pattern = "Hbb-"),]
DE_table = DE_table[!grepl(rownames(DE_table),pattern = "Hba-"),]


pdf("Desktop/Aging_mouse_SVZ/Figures/Volcanoplot.pdf",width = 7,height = 7,useDingbats = F)
par(las=1,bty="l")
plot(DE_table,xaxs='i',yaxs='i',xlim=c(-7,7),ylim=c(0,280),pch=21,
     xlab="Log2FC (Old versus Young)",ylab="-Lgo10(FDR)",cex.lab=1.3,cex.axis=1.3,
     bg=string.to.colors(rownames(DE_table)%in%DOL_genes,colors = c("grey","orange")),cex=1.3)
plot(DE_table,xaxs='i',yaxs='i',xlim=c(-7,7),ylim=c(0,280),pch=21,
     xlab="Log2FC (Old versus Young)",ylab="-Lgo10(FDR)",cex.lab=1.3,cex.axis=1.3,
     bg=string.to.colors(rownames(DE_table)%in%DOL_genes,colors = c("grey","orange")),cex=1.3)
text(DE_table[DOL_genes,],labels =DOL_genes )
dev.off()


