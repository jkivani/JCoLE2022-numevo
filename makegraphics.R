require(rstan)
require(bayesplot)
require(tikzDevice)

x <- NULL

for (i in 1:7) {
  load(paste('rates_ratio_',i,'.Rdata'))
  x <- cbind(x,extract(fit.full)$ratio)
}

props <- apply(x,2,function(x){length(which(x>1))/length(x)})

colnames(x) <- colnames(ST.data)
colnames(x) <- c('Mult. exp.','Phon. rules','Fixed pos.','Not inter.','First slot','Concat','Adj')
colnames(x) <-paste(colnames(x),' (',props,')',sep='')

tikz('numevo_posterior.tex',width=4*2,height=1.5*2)
mcmc_hist(x,facet_args=list(ncol=4)) + vline_at(1)
dev.off()