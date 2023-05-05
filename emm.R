library(Seurat)
pbmc.emm <- Read10X("./single-cell data/EMM_eas/")
pbmc.emm <- CreateSeuratObject(counts = pbmc.emm)
dim(pbmc.emm)
library(readxl)
sample.info1 <- read_xls("./single-cell data/EMM_eas/12276_2022_866_MOESM2_ESM.xls",
                         sheet = "Supplementary Table 4.")
colnames(sample.info1) <- c(t(sample.info1[1,]))
sample.info1 <- sample.info1[-1,] 
sample.info2 <- read_xls("./single-cell data/EMM_eas/12276_2022_866_MOESM2_ESM.xls",
                         sheet = "Supplementary Table 4. (cont'd)")
sample.info <- merge(sample.info1, sample.info2, all = T)
dim(sample.info)
pbmc.emm.barcode <- sample.info[which(sample.info$Disease == "CHIP- Mild COVID-19"),]$BarcodeID
pbmc.emm.cell <- sample.info[which(sample.info$Disease == "CHIP- Mild COVID-19"),]$Celltype
pbmc.emm <- pbmc.emm[,pbmc.emm.barcode]
View(pbmc.emm@meta.data)
pbmc.emm <- AddMetaData(pbmc.emm, pbmc.emm.cell, col.name = "Celltype")
pbmc.emm <- AddMetaData(pbmc.emm, sample.info[which(sample.info$Disease == "CHIP- Mild COVID-19"),]$Sample, col.name = "sample")
#QC and selecting cells
pbmc.emm[["percent.mt"]] <- PercentageFeatureSet(pbmc.emm, pattern = "^MT-")
VlnPlot(pbmc.emm, feature = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
pbmc.emm <- subset(pbmc.emm, subset = nFeature_RNA > 1000 & nFeature_RNA < 15000 & percent.mt < 15)
dim(pbmc.emm)
#normalization
pbmc.emm <- NormalizeData(pbmc.emm, normalization.method = "LogNormalize", scale.factor = 10000)

#feature selection
pbmc.emm <- FindVariableFeatures(pbmc.emm, selection.method = "vst", nfeatures = 2000)
plot1 <- VariableFeaturePlot(pbmc.emm)
top10 <- head(VariableFeatures(pbmc.emm), 10)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = T)
plot1 + plot2

#scaling
pbmc.emm <- ScaleData(pbmc.emm, features = VariableFeatures(pbmc.emm))
#pbmc.emm <- ScaleData(pbmc.emm, vars.to.regress = "percent.mt") 

#PCA
pbmc.emm <- RunPCA(pbmc.emm, features = VariableFeatures(pbmc.emm))
DimPlot(pbmc.emm, reduction = "pca")
DimHeatmap(pbmc.emm, dims = 1:10, cells = 500, balanced = T)
ElbowPlot(pbmc.emm, ndims = 40)

#cell clustering
pbmc.emm <- FindNeighbors(pbmc.emm, dims = 1:20)
pbmc.emm <- FindClusters(pbmc.emm, resolution = 0.3)#0.4-1.2
head(Idents(pbmc.emm), 5)

#umap or tsne
pbmc.emm <- RunUMAP(pbmc.emm, dims = 1:20)
DimPlot(pbmc.emm, reduction = "umap")
Idents(pbmc.emm) <- pbmc.emm@meta.data$labels
DimPlot(pbmc.emm, reduction = "umap", label = T, group.by = c("seurat_clusters", "labels", "Celltype"))

#cell types annotation
library(SingleR)
load("./single-cell data/references/hpca.RData")
annotation <- SingleR(test = pbmc.emm[["RNA"]]@data, ref = hpca.se,
                      assay.type.test=1, labels = hpca.se$label.main)
save(annotation, file = "./single-cell data/EMM_eas/cell.annotation.RData")
pbmc.emm <- AddMetaData(pbmc.emm, annotation$labels, col.name = "labels")
Idents(pbmc.emm) <- pbmc.emm@meta.data$labels
DimPlot(pbmc.emm, group.by = c("seurat_clusters", "labels", "Celltype"), reduction = "umap", label = T, pt.size = 0.5)

#adjustment
pbmc.emm <- subset(pbmc.emm, labels %in% c("B_cell", "Monocyte", "NK_cell", "T_cells"))

table(neutrophil$labels)
pbmc.emm@meta.data$labels[which(pbmc.emm@meta.data$Celltype == "NK cell")] <- "NK_cell"
pbmc.emm@meta.data$labels[which(pbmc.emm@meta.data$Celltype == "CD8 T cell, non-EM-like")] <- "T_cells"
pbmc.emm@meta.data$labels[which(pbmc.emm@meta.data$Celltype == "Monocyte, nonclassical")] <- "Monocyte"

pbmc.emm@meta.data$labels[which(pbmc.emm@meta.data$seurat_clusters == "0")] <- "T_cells"

cluster0 <- subset(pbmc.emm@meta.data, seurat_clusters == "0")
table(cluster0$Celltype)


pbmc.emm.filt <- subset(pbmc.emm, subset = labels == "B_cell" | labels == "Monocyte" |
                          labels == "T_cells" | labels == "NK_cell")
pbmc.emm.filt@meta.data$labels[which(pbmc.emm.filt@meta.data$seurat_clusters == "1")] <- "T_cells"

pbmc.emm.filt.meta <- pbmc.emm.filt@meta.data
save(pbmc.emm.filt.meta, file = "./single-cell data/EMM_eas/pbmc.emm.filt.meta.RData")

pbmc.emm.filt.counts <- pbmc.emm.filt[["RNA"]]@counts
save(pbmc.emm.filt.counts, file = "./single-cell data/EMM_eas/pbmc.emm.filt.counts.RData")
png(filename = "./single-cell data/EMM_eas/cluster.label.ref.png", 
    width = 5000, height = 2000, res = 300)
DimPlot(pbmc.emm, group.by = c("seurat_clusters", "labels"), reduction = "umap", label = T, pt.size = 0.5)

dev.off()
