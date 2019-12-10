setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(bnlearn)
library(igraph)
library(pcalg)

data(insurance)

compute_statistics <- function(real_graph, graph_2, do_print=TRUE){
  tp <- sum(attr(E(graph_2), "vnames") %in% attr(E(real_graph), "vnames"))
  fp <- sum(attr(E(graph_2), "vnames") %in% attr(E(real_graph), "vnames") == FALSE)
  fn <- sum(attr(E(real_graph), "vnames") %in% attr(E(graph_2), "vnames") == FALSE)
#compute_statistics <- function(adj_1, adj_2, do_print=TRUE){
  #tp <- sum(adj_1 == adj_2 & adj_1 != 0 & adj_2 != 0)
  #fp <- sum(adj_1 != adj_2 & adj_2 != 0)
  #fn <- sum(adj_1 != adj_2 & adj_1 != 0)
  precision <- tp / (tp + fp)
  recall <- tp / (tp + fn)
  fscore <- 2 * (precision * recall) / (precision + recall)
  
  if(do_print){
    cat("TP        = ", tp, "\n", sep=" ")
    cat("FP        = ", fp, "\n", sep=" ")
    cat("FN        = ", fn, "\n", sep=" ")
    cat("Precision = ", precision, "\n", sep=" ")
    cat("Recall    = ", recall, "\n", sep=" ")
    cat("Fscore    = ", fscore, "\n", sep=" ")
  }
  
  return(c("tp"=tp, "fp"=fp, "fn"=fn, "precision"=precision, "recall"=recall, 
           "fscore"=fscore))
}

plot_graph <- function(adj, real_graph=NULL, offset_size=15, layout=layout_with_dh){
  graph <- graph_from_adjacency_matrix(adj)
  deg <- degree(graph, mode = "all")
  V(graph)$size <- deg + offset_size
  if(! is.null(real_graph)){
    E(graph)$color <- ifelse(attr(E(graph), "vnames") %in% attr(E(real_graph), "vnames"), "grey","red") 
    E(graph)$width <- ifelse(attr(E(graph), "vnames") %in% attr(E(real_graph), "vnames"), 1, 2) 
  }
  plot(graph, edge.arrow.size=0.4, layout=layout)
}

# Plot with graphviz 
plot_graph2 <- function(bnObject, real_bnObject=NULL, fontsize=3, layoutFun="dot"){
  g <- Rgraphviz::layoutGraph(bnlearn::as.graphNEL(bnObject), layoutType=layoutFun)
  graph::nodeRenderInfo(g) <- list(fontsize=fontsize, shape="circle")
  
  if(! is.null(real_bnObject)){
    # Get skeleton
    bnObject2 <- bnlearn::skeleton(bnObject)
    g2 <- Rgraphviz::layoutGraph(bnlearn::as.graphNEL(bnObject2))
    g2_edge <- graph::edgeRenderInfo(g2)
    
    real_bnObject <- bnlearn::skeleton(real_bnObject)
    real_g <- Rgraphviz::layoutGraph(bnlearn::as.graphNEL(real_bnObject))    
    real_g_edge <- graph::edgeRenderInfo(real_g)
    
    new_col <- graph::edgeRenderInfo(g)$col
    g_names <- names(new_col)
    
    fp_edge <- bnlearn::compare(real_bnObject, bnObject2, arcs=T)$fp
    fp_edge_names <- apply(fp_edge, MARGIN=1, FUN=function(x){ paste(x[1],x[2],sep="~") })
    
    print(g_names[g_names %in% fp_edge_names])
    new_col[g_names %in% fp_edge_names] <- "red"
    
    graph::edgeRenderInfo(g) <- list(col=new_col)
        
    # real_g <- Rgraphviz::layoutGraph(bnlearn::as.graphNEL(real_bnObject))    
    # real_g_names <- names(graph::edgeRenderInfo(real_g)$col)
    # new_col <- graph::edgeRenderInfo(g)$col    
    # g_names <- names(new_col)
    # indexes <- ! g_names %in% real_g_names
    # new_col[indexes] <- "red"
    # graph::edgeRenderInfo(g) <- list(col=new_col)
  }
  
  Rgraphviz::renderGraph(g)
}



#######
# 1
#######
# b
modelstring = paste0("[Age][Mileage][SocioEcon|Age][GoodStudent|Age:SocioEcon]",
                     "[RiskAversion|Age:SocioEcon][OtherCar|SocioEcon][VehicleYear|SocioEcon:RiskAversion]",
                     "[MakeModel|SocioEcon:RiskAversion][SeniorTrain|Age:RiskAversion]",
                     "[HomeBase|SocioEcon:RiskAversion][AntiTheft|SocioEcon:RiskAversion]",
                     "[RuggedAuto|VehicleYear:MakeModel][Antilock|VehicleYear:MakeModel]",
                     "[DrivingSkill|Age:SeniorTrain][CarValue|VehicleYear:MakeModel:Mileage]",
                     "[Airbag|VehicleYear:MakeModel][DrivQuality|RiskAversion:DrivingSkill]",
                     "[Theft|CarValue:HomeBase:AntiTheft][Cushioning|RuggedAuto:Airbag]",
                     "[DrivHist|RiskAversion:DrivingSkill][Accident|DrivQuality:Mileage:Antilock]",
                     "[ThisCarDam|RuggedAuto:Accident][OtherCarCost|RuggedAuto:Accident]",
                     "[MedCost|Age:Accident:Cushioning][ILiCost|Accident]",
                     "[ThisCarCost|ThisCarDam:Theft:CarValue][PropCost|ThisCarCost:OtherCarCost]")
dag = model2network(modelstring)
# c
typeof(dag)
dag
# d
real_adj <- amat(dag)
#real_adj <- amat(cpdag(dag))
# e
real_graph <- graph_from_adjacency_matrix(real_adj)
plot_graph(real_adj)

plot_graph2(dag)


#######
# 2
#######
# b
hc_res <- hc(insurance)
# c
hc_adj <- amat(hc_res)
hc_graph <- graph_from_adjacency_matrix(hc_adj)
# d
plot_graph(hc_adj)
# e
hc_graph <- graph_from_adjacency_matrix(hc_adj)
hc_stats <- compute_statistics(real_graph, hc_graph)
hc_stats <- compute_statistics(graph_from_adjacency_matrix(amat(bnlearn::skeleton(dag))), graph_from_adjacency_matrix(amat(bnlearn::skeleton(hc_res))))
hc_stats <- compute_statistics(graph_from_adjacency_matrix(amat(cpdag(dag))), graph_from_adjacency_matrix(amat(cpdag(hc_res))))
bnlearn::compare(dag, hc_res)
# f
plot_graph(hc_adj, real_graph, offset_size=11)
plot_graph2(hc_res, fontsize=60)
plot_graph2(hc_res, dag, fontsize=60)
temp <- graphviz.compare(dag, hc_res)
# g


#######
# 3
#######
# a
data <- data.matrix(insurance) - 1
nlev=apply(data, MARGIN=2, FUN=function(x){ length(attr(as.factor(x), "levels")) })
pc_res <- pc(suffStat=list(dm = data, nlev=nlev, adaptDF = FALSE), indepTest=disCItest, alpha=0.05, labels=colnames(data))
# b
pc_adj <- as(pc_res, "amat")
# c
plot_graph(pc_adj)
# d
pc_graph <- graph_from_adjacency_matrix(pc_adj)
pc_stats <- compute_statistics(real_graph, pc_graph)
pc_stats <- compute_statistics(graph_from_adjacency_matrix(amat(bnlearn::skeleton(dag))), graph_from_adjacency_matrix(amat(bnlearn::skeleton(as.bn(pc_res)))))
pc_stats <- compute_statistics(graph_from_adjacency_matrix(amat(cpdag(dag))), graph_from_adjacency_matrix(amat(cpdag(as.bn(pc_res)))))
# e
plot_graph(pc_adj, real_graph, offset_size=11, layout=layout_with_dh)
plot_graph2(as.bn(pc_res), dag, fontsize=33)


#######
# 4
#######
# a
aracne_res <- aracne(insurance)
# b
aracne_adj <- amat(aracne_res)
# c
plot_graph(aracne_adj)
# d
aracne_graph <- graph_from_adjacency_matrix(aracne_adj)
aracne_stats <- compute_statistics(real_graph, aracne_graph)
aracne_stats <- compute_statistics(graph_from_adjacency_matrix(amat(bnlearn::skeleton(dag))), graph_from_adjacency_matrix(amat(bnlearn::skeleton(aracne_res))))
aracne_stats <- compute_statistics(graph_from_adjacency_matrix(amat(cpdag(dag))), graph_from_adjacency_matrix(amat(cpdag(aracne_res))))
# e
plot_graph(aracne_adj, real_graph, offset_size=11)
plot_graph2(aracne_res, dag, fontsize=37)
# dot, neato, twopi, circo and fdp
plot_graph2(aracne_res, dag, fontsize=40, layoutFun="fdp")

