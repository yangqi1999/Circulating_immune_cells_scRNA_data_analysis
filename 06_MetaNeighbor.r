options(Seurat.object.assay.version = 'v3')
library(Seurat)
library(MetaNeighbor)
library(SummarizedExperiment)
library(pheatmap)
args <- commandArgs(T)

if(length(args)!=7){
print ("Usage:")
print ("Rscript 06_MetaNeighbor.r sp1_pseudocell_10_.Rds sp1_pseudocell_10_.pheno.csv sp2_pseudocell_10_.Rds sp2_pseudocell_10_.pheno.csv one2one.genelist sp1_name sp2_name")
q(save = "no", status = 0, runLast = TRUE)
}

sp1 <- readRDS(args[1])
sp1_anno <- read.table(args[2],sep=",",header=T)
sp2 <- readRDS(args[3])
sp2_anno <- read.table(args[4],sep=",",header=T)
sp1_sp2 <- read.table( args[5],row.names=1, sep="\t")

comm <- intersect(rownames(sp2),sp1_sp2$V2)
gene_mtx <- sp2[comm,]
sp2 <- CreateSeuratObject(gene_mtx)
sp2$Sample_ID <- sp2_anno$Sample_ID
sp2$Study_ID <- sp2_anno$Study_ID
sp2$Celltype <- sp2_anno$Celltype

comm <- intersect(rownames(sp1),rownames(sp1_sp2))
gene_mtx <- sp1[comm,]
rownames(gene_mtx)<-sp1_sp2[comm,]
sp1 <- CreateSeuratObject(gene_mtx)
sp1$Sample_ID <- sp1_anno$Sample_ID
sp1$Study_ID <- sp1_anno$Study_ID
sp1$Celltype <- sp1_anno$Celltype

all <- merge(sp1,sp2)
data <- SummarizedExperiment(all)
assay(data) <- all$RNA$counts
data$Sample_ID <- all$Sample_ID
data$Study_ID <- all$Study_ID
data$Celltype <- all$Celltype
var_genes = variableGenes(dat = data, exp_labels = data$Study_ID)
write.table(as.data.frame(var_genes),"var_genes.xls",sep="\t",quote=F,row.names=F,col.names=F )
celltype_NV = MetaNeighborUS(var_genes = var_genes,
                             dat = data,
                             study_id = data$Study_ID,
                             cell_type = data$Celltype,
                             fast_version = TRUE)
write.table(celltype_NV,"celltype_NV.txt",sep="\t",quote=F,row.names=T,col.names=T)

top_hits = topHits(cell_NV = celltype_NV,
                   dat = data,
                   study_id = data$Study_ID,
                   cell_type = data$Celltype,
                   threshold = 0.8)
write.table(top_hits,"top_hits_0.8.txt",sep="\t",quote=F,row.names=F,col.names=T)

pdf_name<-paste0(args[6],"_",args[7],"_meta.pdf")
pdf(pdf_name)
cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)
gplots::heatmap.2(celltype_NV,
                  margins=c(8,8),
                  keysize=1,
                  key.xlab="AUROC",
                  key.title="NULL",
                  trace = "none",
                  density.info = "none",
                  col = cols,
                  breaks = breaks,
                  offsetRow=0.1,
                  offsetCol=0.1,
                  cexRow = 0.7,
                  cexCol = 0.7)
dev.off()




