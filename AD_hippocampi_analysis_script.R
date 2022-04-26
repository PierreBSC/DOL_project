library(RColorBrewer)
library(pagoda2)
library(uwot)
library(Matrix)
library(MASS)
library(igraph)
library(pheatmap)
library(fifer)
library(Matrix)
library(biomaRt)
library(ggplot2)


#0)Function definition

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


#I)Data loading and filtering

setwd("GSE153895_RAW/")
Annotation = read.delim("Annotation_data.txt")

list_files = list.files("GSE153895_RAW/",full.names = T)
x = list.files("GSE153895_RAW/",full.names = F)
x = base::strsplit(x,split = "_")
x = unlist(lapply(x,FUN = function(x) x[1]))

data_cout = read.delim(list_files[1])
data_cout = as(as.matrix(data_cout),"dgCMatrix")
Condition =  rep(x[1],ncol(data_cout))
   
for (k in 2:length(list_files)) {
  print(k)
  u = read.delim(list_files[k])
  u = as(as.matrix(u),"dgCMatrix")
  Condition =  c(Condition,rep(x[k],ncol(u)))
  data_cout = cbind(data_cout,u)
}


Condition_full = Annotation[Condition,"genotype group:ch1"]

Lib_size = colSums(data_cout)
Gene_size = rowSums(data_cout)

hist(log10(1+Gene_size),40,main="Total gene UMIs distribution",xlab="Total UMIs (Log10)",
     xaxs="i",col="grey",yaxs='i')
abline(v=log10(50),lty=2,lwd=2,col="red")
hist(log10(1+Lib_size),40,main="Total cellular UMIs distribution",xlab="Total UMIs (Log10)",
     xaxs="i",col="grey",yaxs='i',xlim=c(2,6),ylim=c(0,8000))

data_cout_bis = data_cout[Gene_size>100,]

list_genes_ENSEMBL = rownames(data_cout_bis)

ensembl = useEnsembl(biomart = "ensembl", 
                     dataset = "mmusculus_gene_ensembl")
gene_sequence=getBM(mart = ensembl,attributes = c("ensembl_gene_id","mgi_symbol"),
                    filters = "ensembl_gene_id",values = list_genes_ENSEMBL) 
gene_sequence_bis = gene_sequence$mgi_symbol
names(gene_sequence_bis) = gene_sequence$ensembl_gene_id
gene_sequence_bis = gene_sequence_bis[rownames(data_cout_bis)]
rownames(data_cout_bis) = gene_sequence_bis
x = table(gene_sequence_bis)
data_cout_bis = data_cout_bis[names(which(x==1)),]

#Removing dead/dying cells

Mitochondrial_genes = as.character(rownames(data_cout_bis))
Mitochondrial_genes = Mitochondrial_genes[grepl(Mitochondrial_genes,pattern = "mt-")]

Proportion_mitochondrial_read = colSums(data_cout_bis[Mitochondrial_genes,])/Lib_size

#II)Clustering and cell type assignment


#A)Clustering

r <- Pagoda2$new(data_cout_bis,log.scale=F)
r$adjustVariance(plot=T,gam.k=5)
r$calculatePcaReduction(nPcs=100,odgenes = Selected_genes)
r$makeKnnGraph(k=30,type='PCA',distance = "cosine")
r$getKnnClusters(method=multilevel.community,type='PCA')

r$getDifferentialGenes(type = "PCA",clusterType = "community",verbose = T,z.threshold = 3)
umap_plot = umap(r$reductions$PCA,n_neighbors = 30,spread = 3,
                 n_components = 2,metric = "cosine",verbose = T)

plot(umap_plot,pch=21,bg=string.to.colors(r$clusters$PCA$community),main="Cluster distribution",
     xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2")
plot(umap_plot,pch=21,bg=color_convertion(r$counts[,"C4b"]),main="Cluster distribution",
     xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2")

#B)Meta-clustering

Mean_expression_cluster = aggregate(as.matrix(r$counts[,Selected_genes]),FUN = mean,by=list(r$clusters$PCA$community))
rownames(Mean_expression_cluster) =Mean_expression_cluster$Group.1
Mean_expression_cluster = t(Mean_expression_cluster[,-1])

pdf("Desktop/Genentech_data/Figure_paper/Heatmap_all_clusters.pdf",width = 8,height = 8)
Meta_clustering = pheatmap(cor(Mean_expression_cluster,method = "spearman"),breaks = seq(0,1,length.out=100),
                           clustering_method = "ward.D2",show_colnames = F,show_rownames = T)
dev.off()
Order_cluster = Meta_clustering$tree_col$order

par(las=1,bty="l")
boxplot(r$counts[,"Acta2"]~factor(r$clusters$PCA$community,levels = Order_cluster),outline=F,)
boxplot(r$counts[,"Mbp"]~factor(r$clusters$PCA$community,levels = Order_cluster),outline=F)

Condition_full = Annotation[Condition,"characteristics_ch1"]

#C)Plotting marker gene


Selected_genes_to_plot = c("Hexb","Cx3cr1","Plp1","Mbp","Slc1a2","Gdf10","Acta2","Vtn","Rgs5","Cldn5","Ly6e",
                   "Tmem212","Ccdc153","Ttr","Enpp2","Npy","Nrgn","Col1a2","Cd3d","Pdgfra","Kcnj8","Reln","Pcp4","Cst7","Trem2","Dkk2","Apoe","Ctsb")
data_expression=log2(1+t(data_cout_bis[Selected_genes_to_plot,])/colSums(data_cout_bis)*10^6)

cluster_expression=factor(r$clusters$PCA$community[r$clusters$PCA$community%in%Order_cluster],Order_cluster)

violin_plot_gene=function(gene,title_input) {
  data_gene=data.frame(Expression=as.numeric(data_expression[,gene]),
                       Condition=cluster_expression)
  ggplot(data_gene, aes(x=Condition, y=Expression))  +  geom_violin(trim=TRUE,scale = "width",fill='grey', color="black",bw=1.5)+ 
    scale_y_continuous(name = "Expression log2(TPM)",limits = c(0,max(data_expression[,gene]))) +
    theme_classic() + ggtitle(title_input) + theme(plot.title = element_text(size=22),axis.text=element_text(size = 15),axis.title = element_text(size = 12)) + scale_x_discrete(labels = levels(r$clusters$PCA$community),name=" ")  
}

pdf("Desktop/Genentech_data//Figure_paper//Barplot_marker.pdf",width = 6,height = 4.5)
print(violin_plot_gene('Hexb','Microglia : Hexb'))
print(violin_plot_gene('Cx3cr1','Microglia : Cx3cr1'))
print(violin_plot_gene('Plp1','Oligodendrocytes : Plp1'))
print(violin_plot_gene('Mbp','Oligodendrocytes : Mbp'))
print(violin_plot_gene('Slc1a2','Astrocytes : Slc1a2'))
print(violin_plot_gene('Gdf10','Radial glia : Gdf10'))
print(violin_plot_gene('Acta2','Myocytes : Acta2'))
print(violin_plot_gene('Kcnj8','Pericytes : Kcnj8'))
print(violin_plot_gene('Cldn5','Endothelial cells : Cldn5'))
print(violin_plot_gene('Tmem212','Ependymal cells : Tmem212'))
print(violin_plot_gene('Ccdc153','Ependymal cells : Ccdc153'))
print(violin_plot_gene('Ttr','CP epithelial cells : Ttr'))
print(violin_plot_gene('Enpp2','CP epithelial cells : Enpp2'))
print(violin_plot_gene('Pdgfra','OPCs : Pdgfra'))
print(violin_plot_gene('Col1a2','Fibroblasts : Col1a2'))
print(violin_plot_gene('Cd3d','T-cells : Cd3d'))
print(violin_plot_gene('Nrgn','Neurons : Nrgn'))
print(violin_plot_gene('Reln','CR cells : Reln'))
print(violin_plot_gene('Pcp4','Purkinje : Pcp4'))
dev.off()




###III)Oligo analysis 

Oligo_cluster = c(17,19,23,3,11,16,13,26)
Names_oligo_cells = rownames(r$counts)[which(r$clusters$PCA$community%in%Oligo_cluster)]

Genetic_oligo = Condition_full[which(r$clusters$PCA$community%in%Oligo_cluster)]
Sample_oligo = Condition[which(r$clusters$PCA$community%in%Oligo_cluster)]


data_oligo = data_cout_bis[,Names_oligo_cells]

Zero_proporion_oligo = Matrix::rowSums((data_oligo)==0) / ncol(data_oligo)
Total_expression_oligo = Matrix::rowSums(data_oligo)
plot(log10(Total_expression_oligo),Zero_proporion_oligo[names(Total_expression_oligo)])
Regression_zero_oligo = loess(Zero_proporion_oligo~log10(Total_expression_oligo),subset = Total_expression_oligo>0)
Regression_zero_oligo = Regression_zero_oligo$residuals
Regression_zero_oligo = Regression_zero_oligo[order(Regression_zero_oligo,decreasing = T)]
Selected_oligo_genes = names(Regression_zero_oligo[1:200])
plot(Regression_zero_oligo[1:2000],typle="l")

r_oligo <- Pagoda2$new(data_oligo,log.scale=F)
r_oligo$adjustVariance(plot=T,gam.k=10)
r_oligo$calculatePcaReduction(nPcs=15,odgenes = Selected_oligo_genes )
r_oligo$makeKnnGraph(k=30,type='PCA',distance = "cosine")
r_oligo$getKnnClusters(method=multilevel.community,type='PCA')

x = leiden(r_oligo$graphs)
names(x) = rownames(r_oligo$counts)
r_oligo$clusters$PCA$leiden = x


r_oligo$getDifferentialGenes(type = "PCA",clusterType = "community",verbose = T,z.threshold = 3)
r_oligo$getDifferentialGenes(type = "PCA",clusterType = "leiden",verbose = T,z.threshold = 3)


N_cluster = length(unique(r_oligo$clusters$PCA$leiden))
optimal_palette =colorRamp(brewer.pal(11, "Spectral"))
optimal_palette = optimal_palette((1:N_cluster)/N_cluster)
optimal_palette = optimal_palette / 255
optimal_palette = rgb(optimal_palette)

umap_oligo_plot = umap(r_oligo$reductions$PCA,n_neighbors = 15,pca = NULL,
                 n_components = 2,metric = "cosine",verbose = T)


DOL_cluster = 11
pdf("Desktop/Genentech_data/Figure_paper/UMAP_Oligo_cluster.pdf",useDingbats = F,width = 8,height = )

x = Genetic_oligo =="genotype group: Control"
plot(umap_oligo_plot[x,],pch=21,bg=string.to.colors(r_oligo$clusters$PCA$leiden==DOL_cluster,colors = c("grey70",'orange'))[x],
     xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2",xlim=c(-6.5,6.5),ylim=c(-5,5),main="WT")

x = Genetic_oligo =="genotype group: PS2/APP/P301L"
plot(umap_oligo_plot[x,],pch=21,bg=string.to.colors(r_oligo$clusters$PCA$leiden==DOL_cluster,colors = c("grey70",'orange'))[x],
     xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2",xlim=c(-6.5,6.5),ylim=c(-5,5),main="PS2/APP/P301L")

x = Genetic_oligo =="genotype group: P301L"
plot(umap_oligo_plot[x,],pch=21,bg=string.to.colors(r_oligo$clusters$PCA$leiden==DOL_cluster,colors = c("grey70",'orange'))[x],
     xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2",xlim=c(-6.5,6.5),ylim=c(-5,5),main="P301L")


dev.off()

#


Mean_expression_cluster_oligo = aggregate(as.matrix(r_oligo$counts[,Selected_oligo_genes]),FUN = mean,by=list(r_oligo$clusters$PCA$leiden))
rownames(Mean_expression_cluster_oligo) =Mean_expression_cluster_oligo$Group.1
Mean_expression_cluster_oligo = t(Mean_expression_cluster_oligo[,-1])

Meta_clustering_oligo = pheatmap(cor(Mean_expression_cluster_oligo,method = "spearman"),
                           clustering_method = "ward.D2",show_colnames = F,show_rownames = T)
Order_cluster_oligo = Meta_clustering_oligo$tree_col$order

boxplot(r_oligo$counts[,"Serpina3n"]~factor(r_oligo$clusters$PCA$leiden,Order_cluster_oligo),outline=F)
#

Table_cluster_mice_oligo = table(Sample_oligo,r_oligo$clusters$PCA$leiden)
Table_cluster_mice_oligo = Table_cluster_mice_oligo/rowSums(Table_cluster_mice_oligo)

Table_to_plot_oligo = data.frame(DOL_proportion = Table_cluster_mice_oligo[,"11"]*100,
                           Alternative_DOL = Table_cluster_mice_oligo[,"10"]*100,
                           Inflammed_Oligo = Table_cluster_mice_oligo[,"12"]*100,
                           Condition=Annotation[rownames(Table_cluster_mice_oligo),"genotype group:ch1"])



#Exporting the results of the analysis

pdf("Desktop/Genentech_data/Figure_paper//DOL_proportion.pdf",width = 9,height = 9,useDingbats = F)
par(las=1,bty="l")
prism.plots(DOL_proportion~Condition,Table_to_plot_oligo,cex=3,pch=21,bg="red3",col="black",
            yaxs="i",xaxs="i")
dev.off()

pdf("Desktop/Genentech_data/Result_090421/Alternative_DOL_proportion.pdf",width = 9,height = 9)
par(las=1,bty="l")
prism.plots(Alternative_DOL~Condition,Table_to_plot_oligo,cex=3,pch=21,bg="red3",col="black",
            yaxs="i",xaxs="i")
dev.off()


Normalized_transcript_Serpina3n = log2(1+10^3*r_oligo$misc$rawCounts[,"Serpina3n"]/Lib_size[rownames(r_oligo$misc$rawCounts)])

pdf("Desktop/Genentech_data/Result_090421/Serpina3n_expression_oligo.pdf",width = 12,height = 7)
par(las=1,bty="l")
boxplot(Normalized_transcript_Serpina3n~factor(r_oligo$clusters$PCA$leiden,Order_cluster_oligo),outline=F,
        xlab="Oligodendrocyte cluster",ylab="Normalised Serpina3n expression",yaxs="i",cex.lab=1.2)
dev.off()


write.table(r_oligo$diffgenes$PCA$leiden$`11`,quote=F,
            file = "Desktop/Genentech_data/Result_090421/DE_table_DOL_cluster.txt",sep="\t")
write.table(r_oligo$diffgenes$PCA$leiden$`10`,quote=F,
            file = "Desktop/Genentech_data/Result_090421/DE_table_alternative_DOL_cluster.txt",sep="\t")


##Computing Log2FC and p-value for the DOL signature 

#Log2FC
TPM_data_oligo=log2(1+t(data_oligo[,])/colSums(data_oligo)*10^6)
Mean_log2F_oligo = aggregate(as.matrix(TPM_data_oligo),FUN = mean, by=list(r_oligo$clusters$PCA$leiden==11))
rownames(Mean_log2F_oligo) = Mean_log2F_oligo$Group.1
Mean_log2F_oligo = t(Mean_log2F_oligo[,-1])
Log2_FC_oligo = Mean_log2F_oligo[,2]-Mean_log2F_oligo[,1]

#DE analysis

cloglog_regression = function(x,l,condition) {
  x_bin = x>0
  model_temp = glm(x_bin~1+offset(log(l))+condition,family = binomial(link="cloglog"))
  model_temp = summary(model_temp)
  return(model_temp$coefficients[,4])
}

condition =r_oligo$clusters$PCA$leiden==11
data_oligo_bis = cbind(data_oligo[,condition],data_oligo[,sample(which(!condition),size = sum(condition))])
condition_bis = c(rep("DOL",sum(condition)),rep("Not_DOL",sum(condition)))
Cloglog_result = apply(as.matrix(data_oligo_bis),MARGIN = 1,FUN = function(x) {cloglog_regression(x,colSums(data_oligo_bis),condition_bis)})
Cloglog_result = t(Cloglog_result)


pdf("Desktop/Genentech_data/Figure/Volcanoplot.pdf",width = 9,height = 9,useDingbats = F)
par(las=1,bty="l")
plot(Log2_FC_oligo,-log10(p.adjust(Cloglog_result[,2],method = "BH")),xlim=c(-7,7),ylim=c(0,130),xaxs="i",yaxs="i",
     pch=21,bg=string.to.colors(rownames(Cloglog_result)%in%DOL_genes,colors = c("grey","orange")),
     xlab="Log2FC (DOLs vs other oligos)",ylab="-Log10(FDR)",cex=1.4,cex.lab=1.3,cex.axis=1.3)
plot(Log2_FC_oligo,-log10(p.adjust(Cloglog_result[,2],method = "BH")),xlim=c(-7,7),ylim=c(0,130),xaxs="i",yaxs="i",
     pch=21,bg=string.to.colors(rownames(Cloglog_result)%in%DOL_genes,colors = c("grey","orange")),
     xlab="Log2FC (DOL vs other oligos)",ylab="-Log10(FDR)",cex=1.4,cex.lab=1.3,cex.axis=1.3)
text(Log2_FC_oligo,-log10(p.adjust(Cloglog_result[,2],method = "BH")),labels = names(Log2_FC_oligo))
dev.off()

  View(t(Cloglog_result))
  
pdf("Desktop/Genentech_data/Figure_paper/Log2FC_DOL.pdf",width = 8,height = 8,useDingbats = F)
par(las=1,bty="l")
plot(Mean_log2F_oligo,pch=21,bg="red3",xaxs="i",yaxs='i',
     xlab="Mean expression in homeostatic Oligos (Log2(TPM)",cex.lab=1.2,
     ylab="Mean expression in DOLs (Log2(TPM)",xlim=c(0,17),ylim=c(0,17))
abline(0,1,lwd=2,lty=2)
plot(Mean_log2F_oligo,pch=21,bg="red3",xaxs="i",yaxs='i',
     xlab="Mean expression in homeostatic Oligos (Log2(TPM)",cex.lab=1.2,
     ylab="Mean expression in DOLs (Log2(TPM)",xlim=c(0,17),ylim=c(0,17))
text(Mean_log2F_oligo[,1],Mean_log2F_oligo[,2],labels = rownames(Mean_log2F_oligo))
dev.off()

library(liger)
DOL_genes = read.delim("Documents/DOL_genes.txt")
DOL_genes = DOL_genes$Genes

x = bulk.gsea(values = Log2_FC_oligo,set.list  = list(DOL =DOL_genes),  n.rand = 10000,)
pdf("Desktop/Genentech_data/Figure_paper/Gsea_DOL.pdf",width = 8,height = 8,useDingbats = F)
par(las=1,bty="l")
gsea(values = Log2_FC_oligo,DOL_genes,  n.rand = 10000)
dev.off()


#P-value


condition_vector = r_oligo$clusters$PCA$leiden==11
L_vector = colSums(data_oligo)

Cloglog_regression = function(x,condition_vector,L_vector) {
  x_binarized = factor(as.numeric(x>0))
  model_temp = glm(x_binarized ~ 1 + offset(log(L_vector))+condition_vector,family = binomial(link = "cloglog"))
  model_temp = summary(model_temp)
  return(model_temp$coefficients[2,3])
}
List_z_score = apply(TPM_data_oligo[,],MARGIN = 2,FUN = function(x) {Cloglog_regression(x,condition_vector,L_vector)})
List_z_score = abs(List_z_score)
List_p_value = (1-pnorm(List_z_score,mean = 0,sd = 1))*2

###IV)Microglia

Microglia_cluster = c(7,1,8,12)
Microglia_cells = r$clusters$PCA$community%in%Microglia_cluster

data_microglia = data_cout_bis[,Microglia_cells]

Zero_proporion_microglia = Matrix::rowSums((data_microglia)==0) / ncol(data_microglia)
Total_expression_microglia = Matrix::rowSums(data_microglia)
plot(log10(Total_expression_microglia),Zero_proporion_microglia[names(Total_expression_microglia)])
Regression_zero_microglia = loess(Zero_proporion_microglia~log10(Total_expression_microglia),subset = Total_expression_microglia>0)
Regression_zero_microglia = Regression_zero_microglia$residuals
Regression_zero_microglia = Regression_zero_microglia[order(Regression_zero_microglia,decreasing = T)]

plot(Regression_zero_microglia[1:2000],typle="l")
Selected_microglia_genes = names(Regression_zero_microglia[1:500])

r_microglia <- Pagoda2$new(data_microglia,log.scale=F)
r_microglia$adjustVariance(plot=T,gam.k=10)
r_microglia$calculatePcaReduction(nPcs=20,n.odgenes = Selected_microglia_genes )
r_microglia$makeKnnGraph(k=30,type='PCA',distance = "cosine")
r_microglia$getKnnClusters(method=multilevel.community,type='PCA')
r_microglia$getDifferentialGenes(type = "PCA",clusterType = "community",verbose = T,z.threshold = 3)

x = leiden(r_microglia$graphs)
names(x) = rownames(r_microglia$counts)
r_microglia$clusters$PCA$leiden = x
r_microglia$getDifferentialGenes(type = "PCA",clusterType = "leiden",verbose = T,z.threshold = 3)


umap_microglia_plot = umap(r_microglia$reductions$PCA,n_neighbors = 30,spread = 3,pca = NULL,
                       n_components = 2,metric = "correlation",verbose = T)

plot(umap_microglia_plot,pch=21,bg=string.to.colors(r_microglia$clusters$PCA$leiden%in%c(11)),main="Cluster distribution",
     xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2")


Mean_expression_cluster_microglia = aggregate(as.matrix(r_microglia$counts[,Selected_microglia_genes]),FUN = mean,by=list(r_microglia$clusters$PCA$leiden))
rownames(Mean_expression_cluster_microglia) =Mean_expression_cluster_microglia$Group.1
Mean_expression_cluster_microglia = t(Mean_expression_cluster_microglia[,-1])


Cluster_to_remove_microglia = c(13,15,16,7,6,10)

Meta_clustering_microglia = pheatmap(cor(Mean_expression_cluster_microglia[,-Cluster_to_remove_microglia],method = "spearman"),
                                     clustering_method = "ward.D2",show_colnames = F,show_rownames = T)
Order_cluster_microglia = Meta_clustering_microglia$tree_col$order

#IFN_cluster = 12 
#Microglia with oligodendrocytes  = 6
Selected_microglia_cells = !r_microglia$clusters$PCA$leiden%in%Cluster_to_remove_microglia

umap_microglia_plot = umap(r_microglia$reductions$PCA[Selected_microglia_cells,],n_neighbors = 30,pca = NULL,
                           n_components = 2,metric = "correlation",verbose = T)
plot(umap_microglia_plot,pch=21,bg=string.to.colors(r_microglia$clusters$PCA$leiden[Selected_microglia_cells]),main="Cluster distribution",
          xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2")

##Plot gene expression

Selected_genes_to_plot_microglia = c("Spp1","Cst7","Ctsd","Dkk2")
data_expression=log2(1+t(data_microglia[Selected_genes_to_plot_microglia,Selected_microglia_cells])/colSums(data_microglia[,Selected_microglia_cells])*10^6)

cluster_expression=factor(r_microglia$clusters$PCA$leiden[Selected_microglia_cells],levels = Order_cluster_microglia)

pdf("Desktop/Genentech_data/Figure_paper//DAM_marker.pdf",width = 9,height = 9,useDingbats = F)
print(violin_plot_gene('Cst7','DAM : Cst7'))
print(violin_plot_gene('Dkk2','DAM : Dkk2'))
dev.off()


Sample_microglia = Condition[which(r$clusters$PCA$community%in%Microglia_cluster)]
Sample_microglia = Sample_microglia[Selected_microglia_cells]


Table_cluster_mice_microglia = table(Sample_microglia,r_microglia$clusters$PCA$leiden[Selected_microglia_cells])
Table_cluster_mice_microglia = Table_cluster_mice_microglia/rowSums(Table_cluster_mice_microglia)

Table_to_plot_microglia = data.frame(DAM_proportion = rowSums(Table_cluster_mice_microglia[,c("4","11")])*100,
                                     Alternative_DAM = (Table_cluster_mice_microglia[,c("1")])*100,
                           Condition=Annotation[rownames(Table_cluster_mice_microglia),"genotype group:ch1"])


plot(Table_to_plot_microglia[Table_to_plot_microglia$Condition!="PS2/APP/P301L/TREM2KO","Alternative_DAM"],
     Table_to_plot_oligo[Table_to_plot_oligo$Condition!="PS2/APP/P301L/TREM2KO","DOL_proportion"])


#Exporting the results of the analysis

pdf("Desktop/Genentech_data/Figure_paper//DAM_proportion.pdf",width = 9,height = 9,useDingbats = F)
par(las=1,bty="l")
prism.plots(DAM_proportion~Condition,Table_to_plot_microglia,cex=3,pch=21,bg="red3",col="black",
            yaxs="i",xaxs="i")
dev.off()


###
length(Order_cluster)
for (k in 1:length(Order_cluster)) {
  Selected_cluster = Order_cluster[k]
  Table_temp = r$diffgenes$PCA$community[[Selected_cluster]]
  write.table(Table_temp,file = paste("Desktop/Genentech_data/Table_cell_type/Table_cluster_",k,".txt",sep = ""),sep = "\t",quote=F)
}

