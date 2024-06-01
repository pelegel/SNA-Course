library(sna)
library(network)
library(igraph)
library(intergraph)


filePath=choose.files()
net<-read.csv(filePath,header=TRUE, row.names = 1)


#network
myNet<-as.network.matrix(net, loops=TRUE, multiple=FALSE, ignore.eval=FALSE, names.eval='weight', matrix.type = "adjacency")
myNet
myNet1<-as.network.matrix(net, loops=TRUE, multiple=FALSE,  directed=FALSE, ignore.eval=FALSE, names.eval='weight', matrix.type = "adjacency")
myNet1

#filePath=choose.files()
myUDNet<-as.network.matrix(net, directed=FALSE, loops=TRUE, multiple=FALSE, ignore.eval=FALSE, names.eval='weight', matrix.type = "adjacency")
myUDNet



#igraph
myIgraph<-asIgraph(myNet)
plot(myIgraph)
myIgraph
myUDIgraph<-asIgraph(myUDNet)
plot(myUDIgraph)
#myUDIgraph
#myUDIgraph1<-asIgraph(myNet1)


#----------------------------------- 1 - community detection ------------------------------------#
#edge betweeness
eb <- cluster_edge_betweenness(myUDIgraph)
eb

eb$membership
eb$modularity
modularity(eb)
sizes(eb)
plot_dendrogram(eb)


#louvain
#run louvain 100 times
bestcl<-cluster_louvain(myUDIgraph, resolution=1)
bestmod<-modularity(bestcl)
for (i in 1:100){
  cl<-cluster_louvain(myUDIgraph, resolution=1)
  mod<-modularity(cl)
  if (mod > bestmod){
    bestmod<-mod
    bestcl<-cl
  }
}


cl<-cluster_louvain(myUDIgraph, resolution=1)
modularity(cl)

sizes(bestcl)
cl$membership
cl$modularity

c1 <- c()
c2 <- c()
c3 <- c()
c4 <- c()
c5 <- c()
c6 <- c()
c7 <- c()
c8 <- c()
c9 <- c()
for (i in 1:467){
  if (bestcl$membership[i] == 1){
    c1 <- append(c1, vertex_attr(myUDIgraph, 'vertex.names', index = i))
  }
  if (bestcl$membership[i] == 2){
    c2 <- append(c2, vertex_attr(myUDIgraph, 'vertex.names', index = i))
  }
  if (bestcl$membership[i] == 3){
    c3 <- append(c3, vertex_attr(myUDIgraph, 'vertex.names', index = i))
  }
  if (bestcl$membership[i] == 4){
    c4 <- append(c4, vertex_attr(myUDIgraph, 'vertex.names', index = i))
  }
  if (bestcl$membership[i] == 5){
    c5 <- append(c5, vertex_attr(myUDIgraph, 'vertex.names', index = i))
  }
  if (bestcl$membership[i] == 6){
    c6 <- append(c6, vertex_attr(myUDIgraph, 'vertex.names', index = i))
  }
  if (bestcl$membership[i] == 7){
    c7 <- append(c7, vertex_attr(myUDIgraph, 'vertex.names', index = i))
  }
  if (bestcl$membership[i] == 8){
    c8 <- append(c8, vertex_attr(myUDIgraph, 'vertex.names', index = i))
  }
  if (bestcl$membership[i] == 9){
    c9 <- append(c9, vertex_attr(myUDIgraph, 'vertex.names', index = i))
  }
}


layout_g<-layout_with_fr(myUDIgraph)
par(mfrow=c(1,2))

plot(myUDIgraph, edge.arrow.size=0,vertex.label=NA, 
     layout=layout_g,
     vertex.color=membership(eb),
     vertex.size=igraph::degree(myUDIgraph,mode = "in")+5)

plot(myUDIgraph, edge.arrow.size=0,vertex.label=NA, 
     layout=layout_g,
     vertex.color=membership(cl),
     vertex.size=igraph::degree(myUDIgraph,mode = "in")+5)

plot(eb,myUDIgraph,layout=layout_g)
plot(cl,myUDIgraph,layout=layout_g)
par(mfrow=c(1,1))
myUDIgraph = make_graph(myUDIgraph)
coords = layout_with_fr(myUDIgraph)
# plot the graph
plot(myUDIgraph, layout=coords, vertex.label=NA, vertex.size=10)
c1 = cluster_fast_greedy(myUDIgraph)
modularity(c1)
crossing(c1, myUDIgraph)
plot(c1, myUDIgraph, layout=coords)
plot_dendrogram(c1)



#--------------------------------------- 2 - triad census ---------------------------------------#

tc<- sna::triad.census(myNet, mode="digraph")
dc<- sna::dyad.census(myNet)

#random nets with same dyads
RandNets<-rguman(100,network.size(myNet),network.edgecount(myNet),
                 mut=dc[1], asym = dc[2], null=dc[3],method="exact")
#triad census of the random nets
tcRand<- sna::triad.census(RandNets, mode="digraph")

t.test(tcRand[,1], mu=tc[1])
hist(tcRand[,1],xlim=c(16335649 ,16355258))
abline(v=tc[1],col='blue', lwd=3)

t.test(tcRand[,2], mu=tc[2])
hist(tcRand[,2],xlim=c(51347 ,56911.68))
abline(v=tc[2] ,col='blue', lwd=3)

t.test(tcRand[,3], mu=tc[3])
hist(tcRand[,3],xlim=c(435346,467544.7))
abline(v=tc[3] ,col='blue', lwd=3)

t.test(tcRand[,4], mu=tc[4])
hist(tcRand[,4], xlim=c(15,54))
abline(v=tc[4] ,col='blue', lwd=3)

t.test(tcRand[,5], mu=tc[5])
hist(tcRand[,5], xlim=c(16,341))
abline(v=tc[5],col='blue', lwd=3)

t.test(tcRand[,6], mu=tc[6])
hist(tcRand[,6], xlim=c(31,154))
abline(v=tc[6],col='blue', lwd=3)

t.test(tcRand[,7], mu=tc[7])
hist(tcRand[,7], xlim=c(537.1827,3401))
abline(v=tc[7],col='blue', lwd=3)

t.test(tcRand[,8], mu=tc[8])
hist(tcRand[,8], xlim=c(535.8713,1747))
abline(v=tc[8],col='blue', lwd=3)

t.test(tcRand[,9], mu=tc[9])
hist(tcRand[,9], xlim=c(-0.004018727,3))
abline(v=tc[9],col='blue', lwd=3)

t.test(tcRand[,10], mu=tc[10])
hist(tcRand[,10], xlim=c(-0.00984217,1))
abline(v=tc[10],col='blue', lwd=3)

t.test(tcRand[,11], mu=tc[11])
  hist(tcRand[,11], xlim=c(4437.715,16578))
abline(v=tc[11],col='blue', lwd=3)

t.test(tcRand[,12], mu=tc[12])
hist(tcRand[,12], xlim=c(0.06527686,22))
abline(v=tc[12],col='blue', lwd=3)

t.test(tcRand[13], mu=tc[13])
hist(tcRand[,13], xlim=c(0,58))
abline(v=tc[13],col='blue', lwd=3)

t.test(tcRand[14], mu=tc[14])
hist(tcRand[,14], xlim=c(0,18))
abline(v=tc[14],col='blue', lwd=3)

t.test(tcRand[,15], mu=tc[15])
hist(tcRand[,15], xlim=c(4.563691,324))
abline(v=tc[15],col='blue', lwd=3)

t.test(tcRand[,16], mu=tc[16])
hist(tcRand[,16], xlim=c(13.99579 ,1053))
abline(v=tc[16],col='blue', lwd=3)




#--------------------------------------- 3 - small world ----------------------------------------#

######### is it a smallworld net ###########3

g2Rand<-erdos.renyi.game(gorder(myIgraph), gsize(myIgraph), type ="gnm", directed = FALSE,loops = FALSE)
trg2<-transitivity(myIgraph)

trrand<-transitivity(g2Rand)
tr_ratio<-trg2/trrand ## has to bee greater then 10

CPLg2<-average.path.length(myIgraph)
CPL_rand<-average.path.length(g2Rand)
CPL_ratio<-CPLg2/CPL_rand ## has to be around 1



#--------------------------------------- 4 - scale free -----------------------------------------#

fit_power_law(degree(myIgraph, mode="in")+1,implementation = "plfit")
fit_power_law(degree(myIgraph, mode="out")+1,implementation = "plfit")

plot(degree.distribution(myIgraph, mode="in"), log="xy", xlab = "log(Probability)", ylab="log(degree in)")
plot(degree.distribution(myIgraph, mode="out"), log="xy", xlab = "log(Probability)", ylab="log(degree out)")

hist(degree(myIgraph, mode="in"), col="blue")
hist(degree(myIgraph, mode="out"), col="blue")



g <- barabasi.game(gorder(myIgraph))
d <- degree(g, mode="in")
fit1 <- fit_power_law(d+1, 10)
fit2 <- fit_power_law(d+1, 10, implementation="R.mle")

fit1$alpha
stats4::coef(fit2)
fit1$logLik
stats4::logLik(fit2)


library(ggplot2)
library(dplyr)


