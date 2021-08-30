library(igraph)

df <- read.csv("eth-dataset.csv",colClasses=c("from_address"="character", "to_address"="character"))

nodes <- structure(list(source = df$from_address, destination = df$to_address))

graph <- graph_from_data_frame(nodes)

vertices <- sort(unique(V(graph)))
vertices

## How many components?
no.clusters(graph)

## How big are these?
table(clusters(graph)$csize)

simple_graph <- simplify(graph)

ego_size(graph, order = 1, nodes = V(graph), mode = c("all", "out","in"), mindist = 0)
ego_graph <- make_ego_graph(graph, order = 1, nodes = V(graph), mode = c("all","out", "in"), mindist = 0)
cluster_135 <- ego_graph[[640]]
cluster_55 <- ego_graph[[611]]
cluster_22 <- ego_graph[[573]]
cluster_39 <- ego_graph[[351]]
cluster_21 <- ego_graph[[628]]
cluster_12 <- ego_graph[[606]]
simple_12 <- simplify(cluster_12)
simple_21 <- simplify(cluster_21)
simple_39 <- simplify(cluster_39)
simple_135 <- simplify(cluster_135)
simple_55 <- simplify(cluster_55)
simple_22 <- simplify(cluster_22)

#clustering
kc <- cluster_fast_greedy(as.undirected(simple_graph))
length(kc)
sizes(kc)
membership(kc)
plot(kc,graph,vertex.label=NA,vertex.size = 3,edge.arrow.size=.3)
dendPlot(kc, mode="dendrogram",leaflab = "none")
#dendPlot(kc, mode="dendrogram",leaflab = "none",xlim = c(1, 300))

#drinkChain
kc1 <- cluster_fast_greedy(as.undirected(simple_135))
plot(kc1,cluster_135,vertex.size = 3,edge.arrow.size=.5,vertex.label =  NA)
dendPlot(kc1, mode="phylo")

#Maximine coin
kc2 <- cluster_fast_greedy(as.undirected(simple_55))
plot(kc2,cluster_55,vertex.size = 3,edge.arrow.size=.5,vertex.label =  NA)
dendPlot(kc2, mode="phylo")

#0xaa090b9143c61b08259e185088119c971fa6707c normal address
#kc3 <- cluster_fast_greedy(as.undirected(simple_22))
#plot(kc3,cluster_22,vertex.size = 3,edge.arrow.size=.5,vertex.label =  NA)
#dendPlot(kc3, mode="phylo")

#pointing out from upbit2
kc4 <- cluster_fast_greedy(as.undirected(simple_39))
plot(kc4,cluster_39,vertex.size = 3,edge.arrow.size=.5,vertex.label =  NA)
dendPlot(kc4, mode="phylo")

#IDEX
kc5 <- cluster_fast_greedy(as.undirected(simple_21))
plot(kc5,cluster_21,vertex.size = 3,edge.arrow.size=.5,vertex.label =  NA)
dendPlot(kc5, mode="phylo")

#Crypto.com: MCO Token
kc6 <- cluster_fast_greedy(as.undirected(simple_12))
plot(kc6,cluster_12,vertex.size = 3,edge.arrow.size=.5,vertex.label =  NA)
dendPlot(kc6, mode="phylo")

#mathematical modelling
nv <- vcount(graph)
ne <- ecount(graph)
degs <- igraph::degree(graph)
ntrials <- 12000

num.comm.rg <- numeric(ntrials)
for ( i in ( 1:ntrials)){
  g.rg <- sample_gnm(nv,ne)
  c.rg <- cluster_fast_greedy(g.rg)
  num.comm.rg[i] <-length(c.rg)
}

num.comm.grg <- numeric(ntrials)
for ( i in ( 1:ntrials)){
  g.grg <- sample_degseq(degs,method = "vl")
  c.grg <- cluster_fast_greedy(g.grg)
  num.comm.grg[i] <-length(c.grg)
}

rslts <- c(num.comm.rg,num.comm.grg)
indx <- c(rep(0,ntrials), rep(1,ntrials))
counts <- table(indx,rslts)/ntrials

barplot(counts, beside=TRUE, col=c("blue","red"),
        xlab = "Number of Communities",
        ylab = "Relative Frequency",
        legend=c("Fixed Size", "Fixed Degree Sequence"))


#statistical modelling: networkblock model
library(blockmodels)
set.seed(42)
A.fblog <- as.matrix(as_adjacency_matrix(as.undirected(graph)))
fblog.sbm <- BM_bernoulli("SBM_sym", A.fblog)
fblog.sbm$estimate()
ICLs <- fblog.sbm$ICL
Q <- which.max(ICLs)

Z <- fblog.sbm$memberships[[Q]]$Z
cl.labs <- apply(Z,1,which.max)
nv <- vcount(graph)
summary(Z[cbind(1:nv,cl.labs)])
cl.cnts <- as.vector(table(cl.labs))
alpha <- cl.cnts/nv
alpha
Pi.mat <- fblog.sbm$model_parameters[[Q]]$pi
Pi.mat[3,]

ntrials <- 1000
Pi.mat <- (t(Pi.mat)+Pi.mat)/2 
deg.summ <- list(ntrials)
for(i in (1:ntrials)){ 
  blk.sz <- rmultinom(1,nv,alpha)
  g.sbm <- sample_sbm(nv,pref.matrix=Pi.mat,
                      block.sizes=blk.sz,
                      directed=FALSE)
  deg.summ[[i]] <- summary(igraph::degree(g.sbm))
}
Reduce('+',deg.summ)/ntrials
summary(igraph::degree(graph))
plot(fblog.sbm$ICL,xlab="Q",ylab="ICL",type="b")
lines(c(Q,Q),c(min(ICLs),max(ICLs)),col="red",lty=2)

edges <- as_edgelist(graph,names=FALSE)
neworder<-order(cl.labs)
m<-t(matrix(order(neworder)[as.numeric(edges)],2))
plot(1, 1, xlim = c(0, nv + 1), ylim = c(nv + 1, 0),
     type = "n", axes= FALSE, xlab="Classes",
     ylab="Classes",main="Reorganized Adjacency matrix")
rect(m[,2]-0.5,m[,1]-0.5,m[,2]+0.5,m[,1]+0.5,col=1)
rect(m[,1]-0.5,m[,2]-0.5,m[,1]+0.5,m[,2]+0.5,col=1)
cl.lim <- cl.cnts
cl.lim <- cumsum(cl.lim)[1:(length(cl.lim)-1)]+0.5
clip(0,nv+1,nv+1,0)
abline(v=c(0.5,cl.lim,nv+0.5), 
       h=c(0.5,cl.lim,nv+0.5),col="red")

g.cl <- graph_from_adjacency_matrix(Pi.mat,
                                    mode="undirected",
                                    weighted=TRUE)
vsize <- 100*sqrt(alpha)
ewidth <- 10*E(g.cl)$weight
name <- V(graph)$name
class.by.name <- as.matrix(table(cl.labs,name))
pie.vals <- lapply(1:Q, function(i) 
  as.vector(class.by.name[i,]))
my.cols <- topo.colors(length(unique(name)))
plot(g.cl, edge.width=ewidth,
     vertex.shape="pie", vertex.pie=pie.vals,
     vertex.pie.color=list(my.cols),
     vertex.size=vsize, vertex.label.dist=0.1*vsize,
     vertex.label.degree=pi)

