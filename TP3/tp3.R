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

as_factor_matrix <- function(matrix){
  col <- names(matrix)
  matrix[col] <- lapply(matrix[col], FUN=as.factor)
  return(matrix)
}

#####################
# Plot with igraph
#####################
plot_graph <- function(bnObject, bnObject_adj=NULL, area_ratio=8, size_offset=7, size_ratio=2, highlight_betweenness=FALSE,
                       highlight_link_muta_expr=FALSE){
  if(is.null(bnObject_adj))
    graph <- graph_from_adjacency_matrix(amat(bnObject), mode="directed")
  else
    graph <- graph_from_adjacency_matrix(bnObject_adj, mode="directed")
  
  # Customize colors
  V(graph)$color <- ifelse(names(V(graph)) == tolower(names(V(graph))), "orange",
                           ifelse(names(V(graph)) == toupper(names(V(graph))), "cyan",
                                  "green"))
  
  # Highlist top nodes and arcs according to betweenness centrality
  if(highlight_betweenness){
    # --- nodes
    nodes_condition <- V(graph) %in% V(graph)[rev(order(betweenness(graph)))[1:10]]
    # node shape
    V(graph)$shape <- ifelse(nodes_condition, "square", "circle")
    # node contour color
    V(graph)$frame.color <- ifelse(nodes_condition, "black", "gray")
    # --- edges
    betweenness_edges_condition <- E(graph) %in% E(graph)[rev(order(edge_betweenness(graph)))[1:10]]
    # edge color
    E(graph)$color <- ifelse(betweenness_edges_condition, "black", "gray")
    # edge width
    E(graph)$width <- ifelse(betweenness_edges_condition, 4, 3)
  }
  
  # Highlight link between mutated and expressed genes
  if(highlight_link_muta_expr){
    edges <- t( apply(as.matrix(attr(E(graph), "vnames")), MARGIN=1, FUN=function(x){ unlist(strsplit(x, "[|]")) }) )
    link_muta_expr_edges_condition <- (edges[,1]==toupper(edges[,1]) & edges[,2]==tolower(edges)[,2]) |
                        (edges[,2]==toupper(edges[,2]) & edges[,1]==tolower(edges)[,1])
    # edge color
    #E(graph)$color <- ifelse(link_muta_expr_edges_condition, "red", "gray")
    # node symbol font
    nodes <- unique(c(edges[link_muta_expr_edges_condition]))
    nodes <- nodes[nodes == tolower(nodes)]
    V(graph)$label.font <- ifelse(V(graph) %in% V(graph)[nodes], 2, 1)
  }
  
  # Manage if both highlight are requested
  # if(highlight_betweenness & highlight_link_muta_expr){
  #   E(graph)$color <- ifelse(betweenness_edges_condition & link_muta_expr_edges_condition, "darkred", 
  #                            ifelse(betweenness_edges_condition, "black", 
  #                                   ifelse(link_muta_expr_edges_condition, "red", "gray")))
  # }
  
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
  plot(graph, edge.arrow.size=0.12, layout=layout)
}


#####################
# Plot with graphviz 
#####################
plot_graph2 <- function(bnObject, fontsize=0, layoutFun="dot"){
  g <- Rgraphviz::layoutGraph(bnlearn::as.graphNEL(bnObject), layoutType=layoutFun)
  graph::nodeRenderInfo(g) <- list(fontsize=fontsize, shape="circle")
  
  # Customize colors
  node_names <- names(graph::nodeRenderInfo(g)$col)
  graph::nodeRenderInfo(g)$col[node_names == tolower(node_names)] <- "red"
  graph::nodeRenderInfo(g)$col[node_names == toupper(node_names)] <- "blue"
  
  # Remove node with 0 degree
  deg <- degree(graph_from_adjacency_matrix(amat(bnObject)))  
  graph::nodeRenderInfo(g)$col[deg == 0] <- "white"
  graph::nodeRenderInfo(g)$fontsize[deg == 0] <- 1
  
  # Plot
  Rgraphviz::renderGraph(g)
}


#####################
# Run
#####################
run_algorithm <- function(data, data_cleaned, algo="hc"){
  
}


########################################################
# Network reconstruction with the hill-climbing approach
########################################################
# ------- 2
hc_res <- hc(data)
# We need to remove NaN/NA values

# ------- 3
hc_res <- hc(data_cleaned)
hc_graph <- graph_from_adjacency_matrix(amat(hc_res))

# ------- 4
plot_graph(hc_res, size_offset=7, size_ratio=2, area_ratio=10)
plot_graph2(hc_res)

# ------- 5
# Get hubs
degree <- degree(hc_graph)
sort(degree, decreasing = T)[1:10]

# Get variables related to "Ploidy"
hc_res$nodes$Ploidy$nbr

# Get top 10 nodes and edges with betweenness centrality
# nodes
sort(betweenness(hc_graph), decreasing = T)[1:10]
# edges
E(hc_graph)[rev(order(edge_betweenness(hc_graph)))[1:10]]

plot_graph(hc_res, size_offset=7, size_ratio=2, area_ratio=10, highlight_betweenness=T)

# Get the mutated genes (lower case nodes) that are related to gene expression (upper case nodes)
edges <- hc_res$arcs[ (hc_res$arcs[,1]==toupper(hc_res$arcs[,1]) & hc_res$arcs[,2]==tolower(hc_res$arcs)[,2]) |
               (hc_res$arcs[,2]==toupper(hc_res$arcs[,2]) & hc_res$arcs[,1]==tolower(hc_res$arcs)[,1]), ]
nodes <- unique(c(edges))
nodes <- nodes[nodes == tolower(nodes)]

plot_graph(hc_res, size_offset=7, size_ratio=2, area_ratio=10, highlight_betweenness=T, highlight_link_muta_expr=T)



########################################################
# Network reconstruction with the PC approach
########################################################
# ------- 2
# pc_data <- data
# col <- names(pc_data)
# pc_data[col] <- lapply(pc_data[col], FUN=function(x){ as.integer(x)-1 } )
# 
# nlev=apply(pc_data, MARGIN=2, FUN=function(x){ length(attr(as.factor(x), "levels")) })
# pc_res <- pc(suffStat=list(dm = pc_data, nlev=nlev, adaptDF = FALSE), indepTest=disCItest, alpha=0.05, labels=colnames(pc_data))

# We need to remove NaN/NA values

# ------- 3
pc_data <- data_cleaned
col <- names(pc_data)
pc_data[col] <- lapply(pc_data[col], FUN=function(x){ as.integer(x)-1 } )

nlev=apply(pc_data, MARGIN=2, FUN=function(x){ length(attr(as.factor(x), "levels")) })
pc_res <- pc(suffStat=list(dm = pc_data, nlev=nlev, adaptDF = FALSE), indepTest=disCItest, alpha=0.05, labels=colnames(pc_data))

pc_adj <- as(pc_res, "amat")
pc_graph <- graph_from_adjacency_matrix(pc_adj)

# ------- 4
plot_graph(as.bn(pc_res), size_offset=7, size_ratio=2, area_ratio=10)
plot_graph2(pc_res)

# ------- 5
# Get hubs
degree <- degree(pc_graph)
sort(degree, decreasing = T)[1:10]
# GMPS  CCNE2   ECT2  CENPA   PRC1   RFC4 DIAPH3    RB1  AURKB  AURKA 
# 6      4      4      4      4      4      4      4      4      4 

# Get variables related to "Ploidy"
as.bn(pc_res)$nodes$Ploidy$nbr
# [1] "PPP2R2A"

# Get top 10 nodes and edges with betweenness centrality
# nodes
sort(betweenness(pc_graph), decreasing = T)[1:10]
# GMPS DIAPH3  CCNE2   ECT2   PRC1  CENPA   RFC4  CCNE1   tp53  AURKB 
# 22.0   13.5   12.5   12.5   11.0   10.0   10.0    9.0    8.0    7.0 
# edges
E(pc_graph)[rev(order(edge_betweenness(pc_graph)))[1:10]]
# + 10/99 edges from b7b3893 (vertex names):
#   [1] GMPS  ->RFC4   CCNE2 ->ECT2   DIAPH3->PRC1   ECT2  ->GMPS   tp53  ->GMPS   GPR180->DIAPH3 RB1   ->CCNE2  CCNE1 ->CENPA  GMPS  ->CHEK1  RFC4  ->HRASLS

plot_graph(pc_res, size_offset=7, size_ratio=2, area_ratio=10, highlight_betweenness=T)

# Get the mutated genes (lower case nodes) that are related to gene expression (upper case nodes)
edges <- as.bn(pc_res)$arcs[ (as.bn(pc_res)$arcs[,1]==toupper(as.bn(pc_res)$arcs[,1]) & as.bn(pc_res)$arcs[,2]==tolower(as.bn(pc_res)$arcs)[,2]) |
                        (as.bn(pc_res)$arcs[,2]==toupper(as.bn(pc_res)$arcs[,2]) & as.bn(pc_res)$arcs[,1]==tolower(as.bn(pc_res)$arcs)[,1]), ]
nodes <- unique(c(edges))
nodes[nodes == tolower(nodes)]
# [1] "flt1"   "lin9"   "erbb2"  "tp53"   "scube2"

plot_graph(as.bn(pc_res), size_offset=7, size_ratio=2, area_ratio=10, highlight_betweenness=T, highlight_link_muta_expr=T)


########################################################
# Network reconstruction with the MIIC approach
########################################################

data(cosmicCancer)
data(cosmicCancer_stateOrder)
# execute MIIC (reconstruct graph)
miic.res = miic(inputData = cosmicCancer, categoryOrder = cosmicCancer_stateOrder, latent = TRUE,
                confidenceShuffle = 100, confidenceThreshold = 0.001)

# plot graph
miic.plot(miic.res, igraphLayout=layout_with_drl)

# write graph to graphml format. Note that to correctly visualize
# the network we created the miic style for Cytoscape (http://www.cytoscape.org/).
#miic.write.network.cytoscape(g = miic.res, file = file.path(tempdir(),"/temp"))