require(phangorn)
require(ggplot2)
require(ggtree)

trees <- read.nexus('sinotibetan-sagart-9k-proxy.nex')
tree <- maxCladeCred(trees)

ST.data <- read.csv('numevo_st_temp.csv')
new.data <- ST.data[1,]
for (i in 2:nrow(ST.data)) {
  if (ST.data[i,]$glottocode != ST.data[i-1,]$glottocode) {
    new.data <- rbind(new.data,ST.data[i,])
  }
}
ST.data <- new.data
rownames(ST.data) <- ST.data$glottocode
ST.data <- ST.data[,c('phrules', 'fusion', 'exponence','N_interrupt', 'adjacent', 'fixed_posit', 'slot')]

ST.data <- ST.data[which((rownames(ST.data) %in% tree$tip.label)),]

tree <- keep.tip(tree,which(tree$tip.label %in% rownames(ST.data)))

g <- ggtree(tree) + geom_tiplab(size=3,align=TRUE)

pdf('numevo_tree')
gheatmap(g,apply(ST.data,2,as.factor), offset=1500, width=0.5, font.size=3, 
         colnames_angle=45, hjust=0) + 
  scale_fill_manual(breaks=c("0","1"), 
                    values=c("#00BFC4","#F8766D"),name="classifier")
dev.off()