sce <- readRDS('~/Projects/hcrn/data/sce.rds')

x <- readRDS('~/Projects/hcrn/data/edgeR_immune.rds')

gene.lists <- lapply(x, function(tt){
    IDs <- rownames(tt$table)[tt$table$logFC > 0 & tt$table$FDR < .05]
    unique(rowData(sce)$Symbol[which(rownames(sce) %in% IDs)])
})

# 1 = CD8 T
# 3 = Treg
# 5 = CD4 T

write.table(gene.lists[[1]], file='~/Desktop/CD8.csv', quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(gene.lists[[3]], file='~/Desktop/Treg.csv', quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(gene.lists[[5]], file='~/Desktop/CD4.csv', quote=FALSE, row.names = FALSE, col.names = FALSE)
