setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(bnlearn)
library(igraph)
library(pcalg)
library(miic)
library(qgraph)

data("cosmicCancer")
data <- cosmicCancer
dim(data)

# ------- Remove samples with NA values
data_cleaned <- data[complete.cases(data), ]
dim(data_cleaned)
# We have removed 8 samples with NA values
colSums(is.na(data))[colSums(is.na(data)) != 0]
# The variable with NA is "Ploidy"

# ------- Remove variables with constant value
data_cleaned <- data_cleaned[, apply(data_cleaned, MARGIN=2, FUN=function(x){length(unique(x))}) != 1]
dim(data_cleaned)
# We have removed 14 variables with constant values

# ------- As factor
col <- names(data_cleaned)
data_cleaned[col] <- lapply(data_cleaned[col], FUN=as.factor)


########################################################
# Functions
########################################################

# Plot with igraph
plot_graph <- function(bnObject, area_ratio=8, size_offset=7, size_ratio=2){
  graph <- graph_from_adjacency_matrix(amat(bnObject))
  
  # Customize colors
  V(graph)$color <- ifelse(names(V(graph)) == tolower(names(V(graph))), "orange",
                           ifelse(names(V(graph)) == toupper(names(V(graph))), "cyan",
                                  "green"))
  
  # Customize size
  deg <- degree(graph, mode = "all")  
  V(graph)$size <- log(deg)*size_ratio + size_offset
  
  # Remove node with 0 degree
  graph <- delete.vertices(graph, which(deg==0))    
  
  # Create layout
  e <- get.edgelist(graph, names=F)
  layout <- qgraph.layout.fruchtermanreingold(e, vcount=vcount(graph),
                                              area=area_ratio*(vcount(graph)^2), 
                                              repulse.rad=(vcount(graph)^3.1))
  
  # Plot
  plot(graph, edge.arrow.size=0.07, layout=layout)
}

# Plot with graphviz 
plot_graph2 <- function(bnObject, fontsize=0, layoutFun="dot"){
  g <- Rgraphviz::layoutGraph(bnlearn::as.graphNEL(bnObject), layoutType=layoutFun)
  graph::nodeRenderInfo(g) <- list(fontsize=fontsize, shape="circle")
  
  # Customize colors
  node_names <- names(graph::nodeRenderInfo(g)$col)
  graph::nodeRenderInfo(g)$col[node_names == tolower(node_names)] <- "red"
  graph::nodeRenderInfo(g)$col[node_names == toupper(node_names)] <- "blue"
  
  # Remove node with 0 degree
  deg <- degree(graph_from_adjacency_matrix(amat(hc_res)))  
  graph::nodeRenderInfo(g)$col[deg == 0] <- "white"
  graph::nodeRenderInfo(g)$fontsize[deg == 0] <- 1
  
  # Plot
  Rgraphviz::renderGraph(g)
}


########################################################
# Network reconstruction with the hill-climbing approach
########################################################
# ------- 2
hc_res <- hc(data)
# We need to remove NaN/NA values

# ------- 3
hc_res <- hc(data_cleaned)

# ------- 4
plot_graph(hc_res, size_offset=7, size_ratio=2, area_ratio=10)
plot_graph2(hc_res)



########################################################
# Network reconstruction with the PC approach
########################################################



########################################################
# Network reconstruction with the MIIC approach
########################################################