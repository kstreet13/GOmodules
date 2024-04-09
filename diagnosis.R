# figure out what the important factors are in the UMAP plot
# because we have structure there, but can't see it when we put two cells into the online viewer

require(mclust)
mc <- Mclust(umap, G = 10:20)
clus <- factor(mc$classification)
clus[clus %in% c(17,18)] <- 8
clus <- factor(clus)

plot(umap, asp=1, col=colorby(clus))
pal <- colorby(factor(1:lenu(clus)))
centers <- t(sapply(levels(clus), function(clID){
    colMeans(umap[which(clus==clID),])
}))
points(centers,pch=1,cex=2.5)
points(centers,pch=16,cex=2.5, col=1)
text(centers, labels = levels(clus), col = pal, font=2)



# what distinguishes 5 and 15 from everything else?
ind1 <- which(clus %in% c(5,15))
ind2 <- which(!clus %in% c(5,15))
x <- t(sapply(1:ncol(mod.scores), function(ii){
    t <- t.test(mod.scores[ind1,ii], mod.scores[ind2,ii])
    w <- wilcox.test(mod.scores[ind1,ii], mod.scores[ind2,ii])
    return(c(t$statistic, w$p.value))
}))
x <- as.data.frame(x)
names(x) <- c('stat','pval')
x$pval[x$pval==0] <- min(x$pval[x$pval>0])
x$padj <- p.adjust(x$pval, method = 'fdr')


plot(x$stat, -log10(x$pval))
# three clear outliers
hits <- which(x$stat < -190)
hits <- hits[order(x$pval[hits], decreasing = FALSE)]
# hits <- order(x$stat,decreasing=TRUE)[1:4]


sapply(hits, function(h){
    term <- colnames(mod.scores)[h]
    mi <- which.max(sapply(mods,function(m){m$module_id})==term)
    mods[[mi]]$module_label
})

layout(matrix(1:4,2,2))
plot(umap, asp=1, col=colorby(log1p(mod.scores[,hits[1]]), alpha=.5))
plot(umap, asp=1, col=colorby(log1p(mod.scores[,hits[2]]), alpha=.5))
plot(umap, asp=1, col=colorby(log1p(mod.scores[,hits[3]]), alpha=.5))
plot(umap, asp=1, col=colorby(log1p(mod.scores[,hits[4]]), alpha=.5))
layout(1)
# no perfect markers






hits <- NULL
for(clID in unique(clus)){
    ind1 <- which(clus == clID)
    ind2 <- which(!clus == clID)
    x <- sapply(1:ncol(mod.scores), function(ii){
        t.test(mod.scores[ind1,ii], mod.scores[ind2,ii])$statistic
    })
    h.i <- order(abs(x),decreasing=TRUE)[1:4]
    hits.i <- data.frame(module = h.i, cluster = rep(clID, 4), dir = sign(x[h.i]))
    hits <- rbind(hits, hits.i)
}
hits$mod_label <- sapply(hits$module, function(h){
    term <- colnames(mod.scores)[h]
    mi <- which.max(sapply(mods,function(m){m$module_id})==term)
    mods[[mi]]$module_label
})
# the three biggest drivers are 293, 590, and 712

