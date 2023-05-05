cell.cohort1 <- readRDS("./single-cell data/cell/cohort1.annote.rds")
pbmc.cell <- subset(cell.cohort1, Condition == "Mild")
#QC and selecting cells
#pbmc.cell[["percent.mt"]] <- PercentageFeatureSet(pbmc.cell, pattern = "^MT-")
VlnPlot(cov, feature = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
#pbmc.cell.qc <- subset(pbmc.cell, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#dim(pbmc.cell.qc)
#normalization
#pbmc.cell <- NormalizeData(pbmc.cell, normalization.method = "LogNormalize", scale.factor = 10000)

#feature selection
pbmc.cell <- FindVariableFeatures(pbmc.cell, selection.method = "vst", nfeatures = 2000)

#scaling
pbmc.cell <- ScaleData(pbmc.cell, features = VariableFeatures(pbmc.cell))
#pbmc.cell <- ScaleData(pbmc.cell, vars.to.regress = "percent.mt") 

#PCA
pbmc.cell <- RunPCA(pbmc.cell, features = VariableFeatures(object = pbmc.cell))
DimPlot(pbmc.cell, reduction = "pca")
DimHeatmap(pbmc.cell, dims = 1:10, cells = 500, balanced = T)
ElbowPlot(pbmc.cell, ndims = 40)

#cell clustering
pbmc.cell <- FindNeighbors(pbmc.cell, dims = 1:20)
pbmc.cell <- FindClusters(pbmc.cell, resolution = 0.5)#0.4-1.2
head(Idents(pbmc.cell), 5)

#umap
#pbmc.cell <- RunUMAP(pbmc.cell, dims = 1:20)
DimPlot(pbmc.cell, reduction = "umap")

#cell types annotation
library(SingleR)
load("./single-cell data/references/hpca.RData")
pbmc.cell.data <- GetAssayData(pbmc.cell, slot = "data") #normalized data for annotation
annotation <- SingleR(test = pbmc.cell.data, ref = hpca.se,
                      assay.type.test=1, labels = hpca.se$label.main) 
load("./single-cell data/cell/cell.annotation.RData")
save(annotation, file = "./single-cell data/cell/cell.annotation.RData")
table(annotation$labels)
pbmc.cell <- AddMetaData(pbmc.cell, annotation$labels, col.name = "labels")
Idents(pbmc.cell) <- pbmc.cell@meta.data$label
DimPlot(pbmc.cell, group.by = c("seurat_clusters", "labels", "celltype"), reduction = "umap", label = T, pt.size = 0.5)
#adjustment
pbmc.cell <- subset(pbmc.cell, label %in% c("B_cell", "Monocyte", "NK_cell", "T_cells"))

cell.cluster <- data.frame(pbmc.cell$seurat_clusters, stringsAsFactors = F)
colnames(cell.cluster) <- "cluster"
cell.cluster <- filter(cell.cluster, !grepl(c("4|17|20|19"), cluster))
pbmc.cell <- subset(pbmc.cell, seurat_clusters %in% unique(cell.cluster))

table(pbmc.cell@meta.data$labels)
macrophage <- subset(pbmc.cell@meta.data, labels == "Macrophage")
dim(macrophage)

pbmc.cell@meta.data$label[which(pbmc.cell@meta.data$celltype == "Classical Monocytes")] <- "Monocyte"
pbmc.cell@meta.data$label[which(pbmc.cell@meta.data$celltype == "HLA-DR+ CD83+ Monocytes")] <- "Monocyte"
pbmc.cell@meta.data$label[which(pbmc.cell@meta.data$celltype == "CD163+ Monocytes (1_d7)")] <- "Monocyte"
pbmc.cell@meta.data$label[which(pbmc.cell@meta.data$celltype == "HLA-DR- S100A+ Monocytes")] <- "Monocyte"
pbmc.cell@meta.data$label[which(pbmc.cell@meta.data$celltype == "Non-classical Monocytes")] <- "Monocyte"
pbmc.cell@meta.data$label[which(pbmc.cell@meta.data$celltype == "CD4+ T")] <- "T_cells"
pbmc.cell@meta.data$label[which(pbmc.cell@meta.data$celltype == "CD8+ T")] <- "T_cells"
pbmc.cell@meta.data$label[which(pbmc.cell@meta.data$celltype == "NK")] <- "NK_cell"
pbmc.cell@meta.data$label[which(pbmc.cell@meta.data$celltype == "B")] <- "B_cell"

cluster10 <- subset(pbmc.cell@meta.data, seurat_clusters == "10")
pbmc.cell@meta.data$labels[which(pbmc.cell@meta.data$seurat_clusters == "0")] <- "Monocyte"
table(cluster10$labels)
pbmc.cell.filt <- subset(pbmc.cell, subset = labels == "B_cell" | labels == "Monocyte" | 
                           labels == "NK_cell" | labels == "T_cells")

DimPlot(pbmc.cell.filt, group.by = c("seurat_clusters", "labels"), reduction = "umap", label = T, pt.size = 0.5)

png(filename = "./single-cell data/cell/cluster.label.ref.png", 
    width = 5000, height = 2000, res = 300)
p <- DimPlot(pbmc.cell, group.by = c("seurat_clusters", "label"), reduction = "umap", label = T, pt.size = 0.5)
plot(p)
dev.off()



