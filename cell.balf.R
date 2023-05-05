balf.cell <- readRDS("./single-cell data/cell_balf/BAL.Rds")
dim(balf.cell)
library(org.Hs.eg.db)
balf.cell.ensembl <- rownames(balf.cell)
balf.cell.symbol <- bitr(balf.cell.ensembl, fromType = "ENSEMBL", toType = "SYMBOL", 
                         OrgDb = org.Hs.eg.db, drop = F)
dim(balf.cell.symbol)
symbol.na <- balf.cell.symbol[which(is.na(balf.cell.symbol$SYMBOL) == T),]$ENSEMBL
ensembl.symbol <- read.table("./single-cell data/references/ensembl_symbol.txt",
                             header = T, sep = ",")
symbol.compensate <- subset(ensembl.symbol, Gene.stable.ID %in% balf.cell.ensembl)
#ensembl.symbol[which(ensembl.symbol$Gene.stable.ID %in% balf.cell.ensembl),]$Gene.name
match(balf.cell.ensembl, symbol.compensate$Gene.stable.ID)
match(c(1,2,5,7),c(2,7,9,0))
