library(pagoda2)
library(fifer)
library(uwot)
library(pheatmap)
library(liger)
library(CountClust)


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

data_S1 = read.10x.matrices("Downloads/GSE165098_RAW/S1///")
data_S2 = read.10x.matrices("Downloads/GSE165098_RAW/S2///")
data_S3 = read.10x.matrices("Downloads/GSE165098_RAW/S3///")
data_S4 = read.10x.matrices("Downloads/GSE165098_RAW/S4///")
data_S6 = read.10x.matrices("Downloads/GSE165098_RAW/S6///")
data_S7 = read.10x.matrices("Downloads/GSE165098_RAW/S7///")

l = colnames(data_S1)
l = substr(l,start = 5,stop = 600)
colnames(data_S1) = l

l = colnames(data_S2)
l = substr(l,start = 5,stop = 600)
colnames(data_S2) = l

l = colnames(data_S3)
l = substr(l,start = 5,stop = 600)
colnames(data_S3) = l

l = colnames(data_S4)
l = substr(l,start = 5,stop = 600)
colnames(data_S4) = l

l = colnames(data_S6)
l = substr(l,start = 5,stop = 600)
colnames(data_S6) = l

l = colnames(data_S7)
l = substr(l,start = 5,stop = 600)
colnames(data_S7) = l

data_location = read.delim("tissue_positions_list.csv",sep=",",header = F,row.names = 1)
data_location = data_location[,c(2,3)]
colnames(data_location) = c("Location_X","Location_Y")
data_location$Location_X = data_location$Location_X*6400/78
data_location$Location_Y = data_location$Location_Y*6400/128

data_location_S1 = data_location[colnames(data_S1),]
data_location_S2 = data_location[colnames(data_S2),]
data_location_S3 = data_location[colnames(data_S3),]
data_location_S4 = data_location[colnames(data_S4),]
data_location_S6 = data_location[colnames(data_S6),]
data_location_S7 = data_location[colnames(data_S7),]

data_raw = cbind(data_S1,data_S2,data_S3,data_S4,data_S6,data_S7)
Sample_vector = c(rep("S1",ncol(data_S1)),
           rep("S2",ncol(data_S2)),
           rep("S3",ncol(data_S3)),
           rep("S4",ncol(data_S4)),
           rep("S6",ncol(data_S6)),
           rep("S7",ncol(data_S7)))
data_location_merged = rbind(data_location_S1,data_location_S2,data_location_S3,data_location_S4,data_location_S6,data_location_S7)
data_location_merged$Sample = Sample_vector


Annotation_sample = getGEO("GSE165098")
Annotation_sample = Annotation_sample$GSE165098_series_matrix.txt.gz
View(pData(Annotation_sample))

#B)Filtering

Lib_size = colSums(data_raw)
hist(log10(Lib_size),100)
Gene_size = rowSums(data_raw)
hist(log10(Gene_size+1),100)

data_count = data_raw[Gene_size>100,Lib_size>1000]
rm(data_raw,data_S1,data_S2,data_S3,data_S4,data_S6,data_S7)

colnames(data_count) = paste("Cell",1:ncol(data_count),sep = "_")
x = table(rownames(data_count))

x = names(which(x>1))
data_count = data_count[!rownames(data_count)%in%x,]
data_location_count = data_location_merged[Lib_size>1000,]

#II) Data analysis

#A)Analysis per se

r <- Pagoda2$new(data_count,log.scale=F)
r$adjustVariance(plot=T,gam.k=5)
r$calculatePcaReduction(nPcs=50,n.odgenes=3000)
r$makeKnnGraph(k=30,type='PCA',distance = "cosine")
r$getKnnClusters(method=multilevel.community,type='PCA')

r$getDifferentialGenes(type = "PCA",clusterType = "community",verbose = T,z.threshold = 3)

umap_plot = umap(r$reductions$PCA,n_neighbors = 15,spread = 3,
                 n_components = 2,metric = "cosine",verbose = T)

par(las=1,bty="l")
plot(umap_plot,xlab="UMAP 1",ylab="UMAP 2",pch=21,bg=string.to.colors(r$clusters$PCA$community))


Plot_group_visium = function(sample_number="S1",group=NULL) {
  
  if (is.null(group)) {
    group = r$clusters$PCA$community[data_location_count$Sample==sample_number]
  }
  
  par(las=1,bty="l")
  plot(data_location_count$Location_X[data_location_count$Sample==sample_number],
       data_location_count$Location_Y[data_location_count$Sample==sample_number],pch=21,
       bg=string.to.colors(group),xlab="X position",ylab="Y position")
}


Plot_gene_visium = function(sample_number="S1",gene="Serpina3n") {
  par(las=1,bty="l")
  plot(data_location_count$Location_X[data_location_count$Sample==sample_number],main=gene,
       data_location_count$Location_Y[data_location_count$Sample==sample_number],pch=21,
       bg=color_convertion(r$counts[data_location_count$Sample==sample_number,gene]),xlab="X position",ylab="Y position")
}
boxplot(r$counts[,"Serpina3n"]~data_location_count$Sample,outline=F)

Plot_gene_visium(sample_number = "S7")

#B)LDA based analysis 

Selected_genes = r$getOdGenes(2000)
Selected_genes = Selected_genes[!grepl(Selected_genes,pattern = "mt-")]

Model_LDA_merged = c()
for (k in c(10,15,20,25,30))  {
  print(k)
  Model_LDA_temp = FitGoM(as.matrix(t(data_count[Selected_genes,])),K = k,tol=100,options="BIC")
  Model_LDA_merged[[k]] = Model_LDA_temp
}

List_BIC = unlist(lapply(Model_LDA_merged,FUN = function(x) {x$BIC}))
plot(List_BIC)


Selected_models = Model_LDA_merged[[30]]
Mixing = Selected_models$fit$omega
Contribution = Selected_models$fit$theta
Marginal_distribution = rowSums(data_count[Selected_genes,])/sum(data_count[Selected_genes,])

Get_relevance_table = function(Contribution,Marginal_distribution,lambda=0.5) {
  Relevance_table = apply(Contribution,MARGIN = 2,FUN = function(x) {lambda*log(x) + (1-lambda)*log(x/Marginal_distribution)})
}

Compute_KL = function(p,q) {
  KL_terms = p*log(p/q)
  KL_sums = sum(KL_terms,na.rm = T)
  return(KL_sums)
  
}

Get_specificity = function(Contribution,Marginal_distribution) {
  KL_terms = apply(Contribution,MARGIN = 2,FUN = function(x) {Compute_KL(x,Marginal_distribution)})
  return(KL_terms)
}

Contribution_normalized = Get_relevance_table(Contribution,Marginal_distribution,lambda=0.5)
Topic_specificity = Get_specificity(Contribution,Marginal_distribution)

boxplot(Mixing[,12]~data_location_count$Sample,outline=F,col=c("orange","orange","grey","grey","orange","grey"),
        xlab="Sample",ylab="Inflammatory score",cex.lab=1.3)
x = apply(Contribution_normalized,MARGIN = 2,FUN = rank)

Plot_topic_visium = function(sample_number="S1",topic=1,max_scale=0.08) {
  par(las=1)
  plot(data_location_count$Location_X[data_location_count$Sample==sample_number],main=paste("Sample",sample_number),cex=1,cex.lab=1.5,
       data_location_count$Location_Y[data_location_count$Sample==sample_number],pch=21,col="grey40",xlim=c(0,6400),ylim=c(0,6400),
       bg=color_convertion(Mixing[data_location_count$Sample==sample_number,topic],max_scale =max_scale ),xlab="X position (µm)",ylab="Y position (µm)",xaxs='i',yaxs='i')
}


pdf("Downloads/Figure_Visium_LPS_mice_DOL_topic.pdf",width = 3.5,height = 12)
par(mfrow=c(3,1))
Plot_topic_visium("S1",topic = 21,max_scale=0.08)
Plot_topic_visium("S2",topic = 21,max_scale=0.08)
Plot_topic_visium("S6",topic = 21,max_scale=0.08)
dev.off()


pdf("Downloads/Figure_Visium_LPS_mice_Oligo_topic.pdf",width = 3.5,height = 12)
par(mfrow=c(3,1))
Plot_topic_visium("S1",topic = 2,max_scale=0.5)
Plot_topic_visium("S2",topic = 2,max_scale=0.5)
Plot_topic_visium("S6",topic = 2,max_scale=0.5)
dev.off()
plot(Mixing[,2],Mixing[,21])


x = Contribution_normalized[,12]
x = x[order(x,decreasing = T)]
plot(x[1:50],pch=21,bg="red3",xlab="Rank",ylab="Normalized contribution",cex.lab=1.3,cex=1.3)
text(1:50,x[1:50],labels = names(x[1:50]),pos = 4)

Condition_vector = factor(data_location_count$Sample,levels = c("S1","S2","S6","S3","S4","S7"))

pdf("Downloads/Boxplot_topic_DOL.pdf",width = 6,height = 6,useDingbats = F)
par(las=1,bty="l")
boxplot(Mixing[,21]*100~Condition_vector,outline=F,col=rep(c("orange","grey"),each=3),
        xlab="Mice",ylab="DOL-like score",cex.lab=1.3,ylim=c(0,5),yaxs='i')
legend("topright",fill = c("orange","grey"),bty='n',legend = c("LPS","Control"),cex=1.5)
dev.off()

pdf("Downloads/Boxplot_Serpina3n.pdf",width = 6,height = 6,useDingbats = F)
par(las=1,bty="l")
boxplot(r$counts[,"Serpina3n"]~Condition_vector,outline=F,col=rep(c("orange","grey"),each=3),
        xlab="Mice",ylab="Serpina3n expression",cex.lab=1.3,ylim=c(0,1.5),yaxs='i')
dev.off()


DOL_signature = read.delim("Documents/DOL_genes.txt")
DOL_signature = DOL_signature$Genes
x = Contribution_normalized[,21]
x = (x-mean(x))/sd(x)
m = liger::bulk.gsea(values = x,set.list   =list(DOL_signature) )
m = liger::bulk.gsea(values = x,set.list   =list(DOL_signature) )

pdf("Downloads/GSEA_plot.pdf",width = 6,height = 6,useDingbats = F)
par(las=1)
liger::gsea(values = x,(DOL_signature) )
dev.off()
