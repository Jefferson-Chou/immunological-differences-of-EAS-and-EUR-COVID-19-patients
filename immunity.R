immunity.rds <- readRDS("./single-cell data/immunity/Final_nCoV_0716_upload.rds")
library(stringr)
batches <- stringr::str_subset(unique(immunity.rds@meta.data$batch),"COV")
cov <- subset(immunity.rds, batch %in% batches)
rm(list = "immunity.rds")
#QC and selecting cells
#cov[["percent.mt"]] <- PercentageFeatureSet(cov, pattern = "^MT-")
VlnPlot(cov, feature = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
cov <- subset(cov, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
dim(cov)
cov <- AddMetaData(cov, rep("immunity", dim(cov)[2]), col.name = "source")
##screen empty droplets (usually can be skipped)
#library(DropletUtils)
#ed.test <- emptyDrops(GetAssayData(cov, slot = "counts", assay = "RNA"))

##screen doublets
#library(scran)
#db.test <- doubletCluster(GetAssayData(cov, slot = "counts", assay = "RNA"),
#                          clusters = cov@meta.data$seurat_clusters)
#library(scater)
#doublets <- rownames(db.test)[isOutlier(db.test$N, type = "lower", log = T)]

##check cell cycle genes
#g2m <- CaseMatch(search = cc.genes$g2m.genes, match = rownames(cov)) 
#s <- CaseMatch(search = cc.genes$s.genes, match = rownames(cov)) 
#cov <- CellCycleScoring(object = cov, g2m.features = g2m, 
#                       s.features = s)
#cov <- RunPCA(cov, features = c(g2m, s))
#DimPlot(cov, reduction = "pca", group.by = "Phase")

#normalization
#cov <- NormalizeData(cov, normalization.method = "LogNormalize", scale.factor = 10000)

#feature selection
cov <- FindVariableFeatures(cov, selection.method = "vst", nfeatures = 2000)
#plot1 <- VariableFeaturePlot(cov)
#top10 <- head(VariableFeatures(cov), 10)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = T)
#plot1 + plot2

#scaling
cov <- ScaleData(cov, features = VariableFeatures(cov))
#cov <- ScaleData(cov, vars.to.regress = "percent.mt") 

#PCA
cov <- RunPCA(cov, features = VariableFeatures(object = cov))
#DimPlot(cov, reduction = "pca")
#DimHeatmap(cov, dims = 1:10, cells = 500, balanced = T)
#ElbowPlot(cov, ndims = 40)

#cell clustering
cov <- FindNeighbors(cov, dims = 1:20)
cov <- FindClusters(cov, resolution = 0.5)#0.4-1.2
head(Idents(cov), 5)

#umap or tsne
DimPlot(cov, reduction = "umap", group.by = "seurat_clusters")

#cell types annotation
library(celldex)
hpca.se <- HumanPrimaryCellAtlasData()

library(SingleR)
load("./single-cell data/references/hpca.RData")
load("./single-cell data/immunity/immunity.annotation.RData")
#anno.cov <- SingleR(test = cov[["RNA"]]@data, ref = list(DICE = dice.se, BPE = bpe.se), 
                  # assay.type.test=1, labels = list(dice.se$label.main, bpe.se$label.main))
immunity.annotation <- SingleR(test = cov[["RNA"]]@data, ref = hpca.se,
                      assay.type.test=1, labels = hpca.se$label.main) 
table(immunity.annotation$labels)
cov <- AddMetaData(cov, immunity.annotation$labels, col.name = "labels")
Idents(cov) <- cov@meta.data$labels
DimPlot(cov, group.by = c("seurat_clusters", "labels", "cell_type"), 
        reduction = "umap", label = T, pt.size = 0.5)

#adjustment
table(cov@meta.data$labels)
cov <- subset(cov, labels %in% c("B_cell", "Monocyte", "NK_cell", "T_cells"))
DimPlot(cov.filt, reduction = "umap", label = T, group.by = c("seurat_clusters", "labels"))
#adjust via cell type
neutrophil <- subset(cov@meta.data, labels == "Neutrophils")
table(neutrophil$labels)

cov@meta.data$labels[which(cov@meta.data$labels == "Pro-B_cell_CD34+")] <- "B_cell"
cov@meta.data$labels[which(cov@meta.data$labels == "Pre-B_cell_CD34-")] <- "B_cell"


#adjust via cluster id
cluster1 <- subset(cov@meta.data, seurat_clusters == "1")
table(cluster1$labels)
table(cluster1$cell_type)

cov@meta.data$labels[which(cov@meta.data$seurat_clusters == "1")] <- "T_cells"

cov.filt <- subset(cov, seurat_clusters %in% cov.cluster$cluster)
table(cov.filt@meta.data$labels)
table(cov.filt@meta.data$seurat_clusters)

DimPlot(cov.filt, group.by = c("seurat_clusters", "labels"), 
        reduction = "umap", label = T, pt.size = 0.5)

load("./SNP/gene.to.do.RData")
DotPlot(cov, features = gene.to.do[1:15])

#remove specific clusters
library(dplyr)
cov.cluster <- data.frame(unique(cov$seurat_clusters))
colnames(cov.cluster) <- "cluster"
cov.cluster <- dplyr::filter(cov.cluster, !grepl(c("6|8|10"), cluster)) #select clusters WITHOUT including cluster 6 9 10 11

png(filename = "./single-cell data/immunity/cluster.label.ref.png", 
    width = 5000, height = 2000, res = 300)
p <- DimPlot(cov, group.by = c("seurat_clusters", "labels", "cell_type"), reduction = "umap", label = T, pt.size = 0.5)
plot(p)
dev.off()
