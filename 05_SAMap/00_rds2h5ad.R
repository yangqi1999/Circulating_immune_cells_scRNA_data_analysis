library(optparse)
option_list <- list(
        make_option(c("-n", "--name"), type="character", default=TRUE, help="Input the species name.", metavar="character"),
        make_option(c("-r", "--input1"), type="character", default=TRUE, help="Input the path of the rds file.", metavar="character"),
	make_option(c("-a", "--input2"), type="character", default=TRUE, help="Input the path of the anno file.", metavar="character"),
        make_option(c("-o", "--output"), type="character", default=TRUE, help="Input the path of the first result file -- h5Seurat.", metavar="character")
        );
opt = parse_args(OptionParser(option_list = option_list, usage = "This Script is written to convert from rds to pr_rds and pr_h5ad!"))

library(Seurat)
library(SeuratData)
library(SeuratDisk)

object <- readRDS(opt$input1)
object_refine <- CreateSeuratObject(counts = object, project = opt$name)
anno <- read.table(opt$input2,sep=",",header=T)
anno$Study_ID=NULL
rownames(anno) <- anno$Sample_ID
anno$Sample_ID <- NULL
object_refine$celltype <- anno
SaveH5Seurat(object_refine, filename = opt$output)
Convert(opt$output, dest = "h5ad")


