cell.cohort2 <- readRDS("./single-cell data/cell/cohort2.annote.rds")
pbmc.cell.2 <- subset(cell.cohort2, Condition == "Mild")
dim(pbmc.cell.2)
#QC and selecting cells
#pbmc.cell.2[["percent.mt"]] <- PercentageFeatureSet(pbmc.cell.2, pattern = "^MT-")
VlnPlot(cov, feature = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
#pbmc.cell.2.qc <- subset(pbmc.cell.2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#dim(pbmc.cell.2.qc)
#normalization
#pbmc.cell.2 <- NormalizeData(pbmc.cell.2, normalization.method = "LogNormalize", scale.factor = 10000)

#feature selection
pbmc.cell.2 <- FindVariableFeatures(pbmc.cell.2, selection.method = "vst", nfeatures = 2000)

#scaling
pbmc.cell.2 <- ScaleData(pbmc.cell.2, features = VariableFeatures(pbmc.cell.2))
#pbmc.cell.2 <- ScaleData(pbmc.cell.2, vars.to.regress = "percent.mt") 

#PCA
pbmc.cell.2 <- RunPCA(pbmc.cell.2, features = VariableFeatures(object = pbmc.cell.2))
DimPlot(pbmc.cell.2, reduction = "pca")
DimHeatmap(pbmc.cell.2, dims = 1:10, cells = 500, balanced = T)
ElbowPlot(pbmc.cell.2, ndims = 40)

#cell clustering
pbmc.cell.2 <- FindNeighbors(pbmc.cell.2, dims = 1:20)
pbmc.cell.2 <- FindClusters(pbmc.cell.2, resolution = 0.5)#0.4-1.2
head(Idents(pbmc.cell.2), 5)

#umap
#pbmc.cell.2 <- RunUMAP(pbmc.cell.2, dims = 1:20)
DimPlot(pbmc.cell.2, reduction = "umap")

#cell types annotation
library(SingleR)
load("./single-cell data/references/hpca.RData")
pbmc.cell.2.data <- GetAssayData(pbmc.cell.2, slot = "data") #normalized data for annotation
annotation <- SingleR(test = pbmc.cell.2[["RNA"]]@data, ref = hpca.se,
                      assay.type.test=1, labels = hpca.se$label.main) 
load("./single-cell data/cell/cell.annotation.RData")
save(annotation, file = "./single-cell data/cell/cell.annotation.RData")
table(annotation$labels)
pbmc.cell.2 <- AddMetaData(pbmc.cell.2, annotation$labels, col.name = "labels")
Idents(pbmc.cell.2) <- pbmc.cell.2@meta.data$labels
DimPlot(pbmc.cell.2, group.by = c("seurat_clusters", "labels", "celltype"), reduction = "umap", label = T, pt.size = 0.5)
#adjustment
pbmc.cell.2 <- subset(pbmc.cell.2, label %in% c("B_cell", "Monocyte", "NK_cell", "T_cells"))

cell.cluster <- data.frame(pbmc.cell.2$seurat_clusters, stringsAsFactors = F)
colnames(cell.cluster) <- "cluster"
cell.cluster <- filter(cell.cluster, !grepl(c("4|17|20|19"), cluster))
pbmc.cell.2 <- subset(pbmc.cell.2, seurat_clusters %in% unique(cell.cluster))

table(pbmc.cell.2@meta.data$labels)
macrophage <- subset(pbmc.cell.2@meta.data, labels == "Macrophage")
dim(macrophage)

pbmc.cell.2$label[which(pbmc.cell.2$celltype == "Non-classical Monocytes")] <- "Monocyte"
pbmc.cell.2$label[which(pbmc.cell.2$celltype == "B")] <- "B_cell"

cluster10 <- subset(pbmc.cell.2@meta.data, seurat_clusters == "10")
pbmc.cell.2@meta.data$label[which(pbmc.cell.2@meta.data$seurat_clusters == "18")] <- "B_cell"
table(cluster10$labels)
pbmc.cell.2.filt <- subset(pbmc.cell.2, subset = labels == "B_cell" | labels == "Monocyte" | 
                           labels == "NK_cell" | labels == "T_cells")

DimPlot(pbmc.cell.2.filt, group.by = c("seurat_clusters", "labels"), reduction = "umap", label = T, pt.size = 0.5)

png(filename = "./single-cell data/cell/cluster.label.ref.2.png", 
    width = 5500, height = 2000, res = 300)
p <- DimPlot(pbmc.cell.2, group.by = c("seurat_clusters", "label"), reduction = "umap", label = T, pt.size = 0.5)
plot(p)
dev.off()



