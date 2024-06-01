install.packages('sna')
install.packages('network')
install.packages('igraph')
#install.packages('intergraph')

library(sna)
library(network)
library(igraph)
library(intergraph)


filePath=choose.files()
net<-read.csv(filePath,header=TRUE, row.names = 1)

#network
myNet<-as.network.matrix(net, loops=TRUE, multiple=FALSE, ignore.eval=FALSE, names.eval='weight', matrix.type = "adjacency")
myNet


#igraph
myIgraph<-asIgraph(myNet)
plot(myIgraph)
myIgraph
sum(get.edge.value(myNet, "weight"))


#metrics
vertices <- network.size(myNet)
vertices

edges <- network.edgecount(myNet)
edges


density1<-gden(myNet)       #SNA
density1

density2<-edge_density(myIgraph)    #Igraph
density2

diameter <- diameter(myIgraph)
diameter


far <- farthest_vertices(myIgraph) 
far


betweenness <- which.max(betweenness(myIgraph))
betweenness


betweennessall <- (betweenness(myIgraph))
betweennessall
sort(betweennessall)


closeness <- which.max(closeness(myIgraph))
closeness

maxdegree <- which.max(degree(myIgraph))
maxdegree



#geodesic distribution
distances <- distances(myIgraph, mode = "in")
distances
vec <- as.vector(distances)
vec
geodesic_hist<-hist(vec,breaks=30, xlab = 'Distance', ylab = 'Count', main = 'Geodesic Distribution', col = 'lightblue', freq = FALSE )
#geodist <- geodist(myNet)
#geodist

average_path_length <- average.path.length(myIgraph)
average_path_length




#components
weak_components <- components(myIgraph)$no
weak_components

strong_components <-component.dist(myNet, connected="strong")
strong_components
weak_components <-component.dist(myNet, connected="weak")
weak_components

cutpoints <- cutpoints(myNet)
cutpoints


centralization.degree(myIgraph)

#metrics correlations
degree_betweenness_corr <- cor(degree(myIgraph),betweenness(myIgraph))
degree_betweenness_corr

closeness_betweenness_corr <- cor(closeness(myIgraph),betweenness(myIgraph))
closeness_betweenness_corr

closeness_degree_corr <- cor(closeness(myIgraph),(myIgraph))
closeness_degree_corr




#page rank
pagerank <- page_rank(myIgraph)$vector
pagerank

maxpagerank <- which.max(pagerank)
maxpagerank

centralization.pagerank(myIgraph)


#Eccentricity 
eccentricity <- eccentricity(myIgraph, vids = V(myIgraph), mode = c("in"))
maxeccentricity <- which(eccentricity %in% max(eccentricity))
maxeccentricity








#-----------Pajek-------------
library(sna)
library(network)
library(igraph)


net1<-read.paj(choose.files())
net1
plot(net1)



#igraph
library("igraph")
myIgraph1<-asIgraph(net1)
plot(myIgraph1)
myIgraph1
sum(get.edge.value(net1, "weight"))


#metrics
vertices <- network.size(net1)
vertices

edges <- network.edgecount(net1)
edges


density1<-gden(net1)       #SNA
density1

density2<-edge_density(myIgraph1)    #Igraph
density2

diameter <- diameter(myIgraph1)
diameter


far <- farthest_vertices(myIgraph1) 
far


#geodesic distribution
distances <- distances(myIgraph1, mode = "in")
distances
vec <- as.vector(distances)
vec

distances_vec <- c()
for (x in 1:length(vec)) {
  if (vec[x] != Inf){
      distances_vec<-append(distances_vec, vec[x])    
  }
}

mean_geodesic<-mean(distances_vec)
sd_geodesic<-sd(distances_vec)


mean_distance<- sd(vec)
mean_distance
geodesic_hist<-hist(vec, xlab = 'Distance', ylab = 'Count', main = 'Geodesic Distribution', col = 'lightblue', freq = FALSE )
#geodist <- geodist(myNet)
#geodist


average_path_length <- average.path.length(myIgraph1)
average_path_length



#components
weak_components <- components(myIgraph1, mode="weak")$no
weak_components

strong_components <- components(myIgraph1, mode="strong")$no
strong_components



#isolates
isolates<-isolates(net1)
isolates

#cutpoints
cutpoints <- cutpoints(net1)
cutpoints


#central nodes
betweenness <- which.max(betweenness(myIgraph1))
betweenness

betweennessall <- (betweenness(myIgraph1))
betweennessall
sort(betweennessall)


maxdegree <- which.max(degree(myIgraph1, mode="total"))
maxdegree


closeness<-closeness(myIgraph1)
closeness_vex<-c()
for (x in 1:length(closeness)){
  if (closeness[x] != max(closeness)){
    closeness_vex<-append(closeness_vex,closeness[x])
  }
}

closeness <- which.max(closeness_vex)
closeness


centralization.degree(myIgraph)

#metrics correlations
degree_betweenness_corr <- cor(degree(myIgraph1),betweenness(myIgraph1))
degree_betweenness_corr

closeness_betweenness_corr <- cor(closeness(myIgraph1),betweenness(myIgraph1))
closeness_betweenness_corr

closeness_degree_corr <- cor(closeness(myIgraph1),degree(myIgraph1))
closeness_degree_corr




#page rank
pagerank <- as.numeric(page_rank(myIgraph1)$vector)
pagerank
maxpagerank <- which.max(pagerank)
maxpagerank

centralize(pagerank)
centralize(eccentricity)


degree_centralization<-centralization.degree(myIgraph1, mode="in")$centralization
degree_centralization
betweenness_centralization<-centralization.betweenness(myIgraph1)$centralization
betweenness_centralization


#Eccentricity 
eccentricity <- eccentricity(myIgraph1, vids = V(myIgraph1), mode = c("in"))
maxeccentricity <- which(eccentricity %in% max(eccentricity))
maxeccentricity
mineccentricity <- which(eccentricity %in% min(eccentricity))
mineccentricity

centralize(eccentricity)
centralize(eccentricity)


rgraph_degrees<-c(0)
rgraph_betweennesses<-c(0)

#random graph
for (x in 1:100) {
  rgnm<-rgnm(1,467,2179, mode="digraph")
  rgraph<-graph_from_adjacency_matrix(rgnm, mode = "directed")
  rgraph_degrees<-append(rgraph_degrees,centralization.degree(rgraph, mode="in")$centralization,x)
  rgraph_betweennesses<-append(rgraph_betweennesses,centralization.betweenness(rgraph)$centralization,x)
}

mean_rgraph_degree<-mean(rgraph_degrees)
mean_rgraph_degree
mean_rgraph_betweenness<-mean(rgraph_betweennesses)
mean_rgraph_betweenness


degree_centralization/mean_rgraph_degree
betweenness_centralization/mean_rgraph_betweenness
