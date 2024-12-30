library(Seurat)
library(dplyr)
library(cowplot)
library(reshape2)
args = commandArgs(T)
ReadFile <- function(x){
        x <- as.vector(x)
        x1 <- x[1]
        x2 <- x[2]
        FileTmp <- readRDS(x2)
        FileTmp$stim <- x1
        return(FileTmp)
}

FileName <- read.table(args[1]) ###Each library path file
FileList <- apply(FileName,1,ReadFile)

pancreas.features <- SelectIntegrationFeatures(object.list = hm_list, nfeatures = 3000)
pancreas.list <- PrepSCTIntegration(object.list = hm_list,anchor.features = pancreas.features)
pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, normalization.method = "SCT", anchor.features = pancreas.features)
combined <- IntegrateData(anchorset = pancreas.anchors, normalization.method = "SCT")

pdf("Vlnplot.pdf",width=15,height=10)
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, group.by = "stim")
dev.off()
combined <- RunPCA(object = combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(object = combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(object = combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)
pdf("align_Clusters.pdf",width=20,height=10)
p1 <- DimPlot(object = combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(object = combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
dev.off()
pdf("dimplot_sample.pdf",width=30,height=5)
DimPlot(object = combined, reduction = "umap", split.by = "stim")
dev.off()
saveRDS(combined, file = "combined_analysis.rds")

DefaultAssay(combined)<-"RNA"
combined <- NormalizeData(combined)
all.genes <- rownames(combined)
combined <- ScaleData(combined,features=all.genes)
combined.markers <- FindAllMarkers(object = combined, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5, test.use ="wilcox")
topn <- combined.markers %>% group_by(cluster)
topn <- select(topn, gene, everything())
topn$gene=gsub("-","_",topn$gene)
write.table(topn[topn$p_val_adj <= 0.05,],file="RNA_GeneDiffExpFilter.xls",sep="\t",col.names = TRUE,row.names = F,quote=F)

