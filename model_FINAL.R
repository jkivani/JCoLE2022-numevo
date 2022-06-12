require(phytools)
require(rstan)
#require(cmdstanr)

set.seed(1234)

k = as.integer(commandArgs(trailingOnly = TRUE)[1])

model.code <- "functions {
  matrix evprob(real z, real alpha, real beta) {
    matrix[2,2] P;
    P[1,1] = (beta/(alpha+beta)) + (alpha/(alpha+beta)*exp(-(alpha+beta)*z));
    P[1,2] = (alpha/(alpha+beta)) - (alpha/(alpha+beta)*exp(-(alpha+beta)*z));
    P[2,1] = (beta/(alpha+beta)) - (beta/(alpha+beta)*exp(-(alpha+beta)*z));
    P[2,2] = (alpha/(alpha+beta)) + (beta/(alpha+beta)*exp(-(alpha+beta)*z));
    return P;
  }
  //compute likelihood via Felsenstein's Pruning Algorithm
  real pruning(int N, int B, int[] child, int[] parent, real[] brlen, matrix tiplik, real alpha, real beta) {
    matrix[N,2] lambda;                  //likelihoods at tips+nodes
    vector[2] pi;                         //stationary probability
    lambda = log(tiplik);
    for (b in 1:B) {
      matrix[2,2] P = evprob(brlen[b], alpha, beta); //via matrix exponentiation
      for (d in 1:2) {
        lambda[parent[b],d] += log(dot_product(P[d],exp(lambda[child[b]])));
      }
    }
    pi[1] = log(beta) - log(alpha+beta) + lambda[parent[B],1];
    pi[2] = log(alpha) - log(alpha+beta) + lambda[parent[B],2];
    return(log_sum_exp(pi));
  }
}
data {
  int<lower=1> N; //number of tips+internal nodes+root
  int<lower=1> T; //number of tips
  int<lower=1> B; //number of branches
  int<lower=1> child[B];                //child of each branch
  int<lower=1> parent[B];               //parent of each branch
  real<lower=0> brlen[B];               //length of each branch
  matrix[N,2] tiplik;     //likelihoods for data at tips in tree
  }
parameters {
  real<lower=0,upper=1> p_diff;
  real<lower=0> s[3];
}
transformed parameters {
  vector[2] llik;
  llik[1] = log(1-p_diff) + pruning(N,B,child,parent,brlen,tiplik,s[1],s[1]);
  llik[2] = log(p_diff) + pruning(N,B,child,parent,brlen,tiplik,s[2],s[3]);
}
model {
  p_diff ~ beta(1,1);
  s ~ gamma(1,1);
  target += log_sum_exp(llik);
}
generated quantities {
  int z;
  real ratio;
  z = bernoulli_rng(exp(llik[2]-log_sum_exp(llik)));
  if (z == 0) {
    ratio = 1;
  }
  if (z == 1) {
    ratio = s[2]/s[3];
  }
}"

trees <- read.nexus('sinotibetan-sagart-9k-proxy-numevo.nex')

ST.data <- read.csv('numevo_sinotibetan.csv')
ST.data <- ST.data[,c(1,4:10)]
ST.data <- aggregate(.~glottocode,ST.data,sum)
rownames(ST.data) <- ST.data$glottocode

ST.data <- ST.data[,2:8]
ST.data[ST.data>1] = 1
ST.data <- ST.data[,order(colMeans(ST.data))]

inds <- sample(1:length(trees),50)

bin.states <- cbind(to.matrix(as.character(ST.data[,k]),seq=c('0','1')))
rownames(bin.states) <- rownames(ST.data)

fit.list <- list()
for (t in inds) {
  tree <- trees[[t]]
  tree <- keep.tip(tree,which(tree$tip.label %in% rownames(bin.states)))
  tree <- reorder.phylo(tree,'pruningwise')
  bin.states.t <- bin.states[tree$tip.label,]
  D <- ncol(bin.states.t)
  bin.states.t <- rbind(bin.states.t,matrix(1,nrow=tree$Nnode,ncol=ncol(bin.states.t))) #likelihoods for tips + nodes
  parent <- tree$edge[,1]
  child <- tree$edge[,2]
  tree$edge.length[tree$edge.length==0] <- 1
  b.lens <- tree$edge.length/1000
  N <- length(unique(c(parent,child)))
  T <- length(child[which(!child %in% parent)])
  tip.lik <- bin.states.t
  data.list <- list(N=N,
                    T=T,
                    B=length(parent),
                    brlen=b.lens,
                    child=child,
                    parent=parent,
                    tiplik=tip.lik)
  fit <- stan(model_code=model.code,data=data.list,thin=10)
  fit.list <- append(fit.list,fit)
  fit.full <- sflist2stanfit(fit.list)
  save.image(paste('rates_ratio_',k,'.Rdata'))
}
