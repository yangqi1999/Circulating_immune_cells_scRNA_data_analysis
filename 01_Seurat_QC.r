library(Seurat)
library(dplyr)
library(cowplot)
library(reshape2)
library(DoubletFinder)
library(data.table)
library(ggplot2)
args<-commandArgs(T)
if(length(args)!=6){
print ("Usage:")
print ("cd work_dir ;Rcript anno.R count_matrix sample max_nFeature_RNA max_mt doublet_rate")
print ("Note:")
print ("Aim to run single-sample umap")
q(save = "no", status = 0, runLast = TRUE)
}

setwd(args[1])
object.data <- Read10X(data.dir=args[2],gene.column = 1)
sample <- args[3]
colnames(object.data)<-paste(colnames(object.data),sample,sep="_")
object_name <- CreateSeuratObject(counts = object.data, project =sample, min.cells = 3, min.features = 200)

object_name[["percent.mt"]] <- PercentageFeatureSet(object = object_name, pattern = "ATP6|COX1|COX2|COX3|CYTB|ND1|ND2|ND3|ND4|ND4L|ND5|ND6")
## draw the count of various mt.percent

if (sum(object_name@meta.data$percent.mt)== 0){
  object_name@meta.data$percent.mt[1]<-0.000001
}
object_name <- subset(x = object_name, subset = nFeature_RNA >= 200 & nFeature_RNA <= as.numeric(args[4]) & nCount_RNA >= 500 & percent.mt <= as.numeric(args[5]))
object_name <- SCTransform(object_name, verbose = FALSE)
object_name <- RunPCA(object_name, features = VariableFeatures(object = object_name))
object_name <- FindNeighbors(object_name, dims = 1:30)
object_name <- FindClusters(object_name, resolution = 0.5)
object_name <- RunUMAP(object_name, dims = 1:30)
sweep.res.list_SM <- paramSweep_v3(object_name, PCs = 1:30,sct=TRUE)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
homotypic.prop <- modelHomotypic(object_name@meta.data$SCT_snn_res.0.5)
nExp_poi <- round(as.numeric(args[6])*ncol(object_name))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
object_name <- doubletFinder_v3(object_name, PCs = 1:30, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = FALSE,sct=TRUE)
DFTitle <- grep(pattern="DF.classifications",colnames(object_name@meta.data),value=TRUE)
object_name$doublet_finder <- object_name@meta.data[,DFTitle]
pdf("doublet.out.umap.pdf",width=12,height=10)
print(DimPlot(object_name, reduction = "umap", group.by = "doublet_finder"))
print(DimPlot(object_name, reduction = "umap", group.by = "SCT_snn_res.0.5"))
dev.off()
object_Singlet <- subset(object_name,doublet_finder == "Singlet")
object_Singlet<-SCTransform(object_Singlet, verbose = FALSE)

pdf(sprintf("%s_Vlnplot.pdf",sample),width=15,height=10)
VlnPlot(object_Singlet,features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
dev.off()
object_Singlet <- RunPCA(object = object_Singlet, verbose = FALSE)
pdf(sprintf("%s_DimHeatmap_30pc.pdf",sample),width=15,height=10)
DimHeatmap(object = object_Singlet, dims = 1:30, cells= 200, balanced = TRUE)
dev.off()
object_Singlet <- RunUMAP(object = object_Singlet, reduction = "pca", dims = 1:30)
object_Singlet <- FindNeighbors(object = object_Singlet, reduction = "pca", dims = 1:30)
object_Singlet <- FindClusters(object_Singlet, resolution = 0.5)
pdf(sprintf("%s_sample-umap.pdf",sample),width=15,height=10)
DimPlot(object =object_Singlet, reduction = "umap")
dev.off()
saveRDS(object_Singlet, file = "combined_analysis.rds")
