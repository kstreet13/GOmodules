
require(rjson)
mods <- fromJSON(file='modules/clean_bpmodules.json')

# section
# category (row)
# module (box)
# node (module components) (gene list)
# gene(s)

# extract the modules, nodes, and genes
# some modules have 0 genes
df <- lapply(mods, function(mc){
    if(length(mc$leaf_genes) > 0){
        return(
            data.frame(
                # can't rely on labels, because some are ""
                module = mc$module_id,
                node = mc$node_id,
                gene = sapply(mc$leaf_genes, function(g){ g$gene_symbol })
            )
        )
    }else{
        return(NULL)
    }
})
df <- do.call(rbind, df)

# read in single-cell data
require(SingleCellExperiment)
sce <- readRDS('~/Projects/hcrn/data/sce.rds')
rownames(sce) <- rowData(sce)$Symbol

###
# use counts matrix to build module % activation vector for each cell
###
# subset to relevant genes
counts <- assay(sce,'counts')[rownames(sce) %in% df$gene, ]

# logical matrix indicating if node is active in given cell
node.active <- sapply(unique(df$node), function(nd){
    colSums(counts[which(rownames(counts) %in% df$gene[which(df$node == nd)]), , drop=FALSE]) > 0
})
# module % activation per cell
mod.scores <- sapply(unique(df$module), function(md){
    rowMeans(node.active[, which(colnames(node.active) %in% df$node[which(df$module == md)]), drop = FALSE])
})

geneNorm <- mod.scores / rowSums(mod.scores) # somehow much worse?
logCountNorm <- mod.scores / log(sce$total_counts) # not bad

{
    pca <- BiocSingular::runPCA(mod.scores, rank = nrow(mod.scores))
    resids <- lm(pca$x ~ log(sce$total_counts))$residuals
    pca <- BiocSingular::runPCA(resids, rank = nrow(mod.scores))
    umap1 <- uwot::umap(pca$x)
    plot(umap1, asp=1, col = colorby(factor(sce$seurat_clusters_no2s)))
    plot(umap1, asp=1, col = colorby(log1p(sce$total_genes)))
    
    umap2 <- uwot::umap(pca$x, metric = "cosine")
    plot(umap2, asp=1, col = colorby(factor(sce$seurat_clusters_no2s)))
    plot(umap2, asp=1, col = colorby(log1p(sce$total_genes)))
    
    
} # regress out total count on PCs



# PCA
pca <- BiocSingular::runPCA(logCountNorm, rank = 5)
# strongly impacted by sequencing depth
plot(pca$x,asp=1, col=colorby(log1p(sce$total_counts)))
plot(pca$x,asp=1, col=colorby(log1p(sce$total_genes)))
# but still some biological signal
plot(pca$x,asp=1, col=colorby(factor(sce$seurat_clusters_no2s)))
plot(pca$x[,2:3],asp=1, col=colorby(factor(sce$seurat_clusters_no2s)))

pairs(pca$x[,2:5],asp=1, col=colorby(factor(sce$seurat_clusters_no2s)))

# need a way to normalize for sequencing depth/number of genes
# regress it out?

# may be a good case for NMF

# straight to UMAP?
umap1 <- uwot::umap(logCountNorm)
plot(umap1, asp=1, col = colorby(factor(sce$seurat_clusters_no2s)))
plot(umap1, asp=1, col = colorby(log1p(sce$total_genes)))

umap2 <- uwot::umap(logCountNorm, metric = "cosine")
plot(umap2, asp=1, col = colorby(factor(sce$seurat_clusters_no2s)))
plot(umap2, asp=1, col = colorby(log1p(sce$total_genes)))

