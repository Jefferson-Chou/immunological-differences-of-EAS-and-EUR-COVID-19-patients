# data reading
library(Seurat)
for(i in 1:length(list.files(path = "./single-cell data/Nat. Med/", pattern = ".h5$"))){
  balf.nm <- Read10X_h5(paste0("./single-cell data/Nat. Med/", list.files(path = "./single-cell data/Nat. Med/", pattern = ".h5$")[i]),
                        use.names = T, unique.features = T) %>% CreateSeuratObject(counts = , min.cells = 3, min.features = 200) %>% 
    AddMetaData(balf.nm, strsplit(list.files(path = "./single-cell data/Nat. Med/", pattern = ".h5$")[i], "_")[[1]][2], col.name = "sample")
  assign(paste0("balf.nm.", strsplit(list.files(path = "./single-cell data/Nat. Med/", pattern = ".h5$")[i], "_")[[1]][2]), balf.nm)
  rm(list = "balf.nm")
}

balf.nm.GSM3660650 <- Read10X("./single-cell data/Nat. Med/GSM3660650/") %>% 
  CreateSeuratObject(counts = , min.cells = 3, min.features = 200) %>% 
  AddMetaData(object = , "GSM3660650", col.name = "sample")

balf.nm.C141$group <- "mild"
balf.nm.C142$group <- "mild"
balf.nm.C143$group <- "severe"
balf.nm.C144$group <- "mild"
balf.nm.C145$group <- "severe"
balf.nm.C146$group <- "severe"
balf.nm.C148$group <- "severe"
balf.nm.C149$group <- "severe"
balf.nm.C152$group <- "severe"
balf.nm.C51$group <- "healthy"
balf.nm.C52$group <- "healthy"
balf.nm.C100$group <- "healthy"
balf.nm.GSM3660650$group <- "healthy"