setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(bnlearn)
library(igraph)
library(pcalg)

data(insurance)

compute_statistics <- function(real_graph, graph_2, do_print=TRUE){
  tp <- sum(E(graph_2) %in% E(real_graph))
  fp <- sum(E(graph_2) %in% E(real_graph) == FALSE)
  fn <- sum(E(real_graph) %in% E(graph_2) == FALSE)
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

plot_graph <- function(adj, real_graph=NULL, offset_size=15){
  graph <- graph_from_adjacency_matrix(adj)
  deg <- degree(graph, mode = "all")
  V(graph)$size <- deg + offset_size
  if(! is.null(real_graph)){
    E(graph)$color <- ifelse(E(graph) %in% E(real_graph), "grey","red") 
    E(graph)$width <- ifelse(E(graph) %in% E(real_graph), 1, 2) 
  }
  plot(graph, edge.arrow.size=0.4, layout=layout_with_dh)
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


#######
# 2
#######
# b
hc_res <- hc(insurance)
# c
hc_adj <- amat(hc_res)
# d
plot_graph(hc_adj)
# e
hc_stats <- compute_statistics(real_graph, hc_graph)
hc_stats <- compute_statistics(graph_from_adjacency_matrix(amat(cpdag(dag))), graph_from_adjacency_matrix(amat(cpdag(hc_res))))
hc_stats <- compute_statistics(graph_from_adjacency_matrix(amat(bnlearn::skeleton(dag))), graph_from_adjacency_matrix(amat(bnlearn::skeleton(hc_res))))
# f
plot_graph(hc_adj, real_graph, offset_size=11)
# g


#######
# 3
#######
# a
data <- data.matrix(insurance) - 1
pc_res <- pc(suffStat=list(C = cor(data), n=nrow(data)), indepTest=disCItest, alpha=0.05, labels=colnames(data))
# b
pc_adj <- amat(pc_res)
# c
pc_graph <- graph_from_adjacency_matrix(pc_adj)
plot(pc_graph, edge.arrow.size=0.1, layout=layout_with_lgl)
# d
hc_stats <- compute_statistics(real_graph, pc_graph)
# e



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
#aracne_stats <- compute_statistics(amat(cpdag(dag)), amat(cpdag(aracne_res)))
aracne_stats <- compute_statistics(real_graph, aracne_graph)
aracne_stats <- compute_statistics(graph_from_adjacency_matrix(amat(cpdag(dag))), graph_from_adjacency_matrix(amat(cpdag(aracne_res))))
aracne_stats <- compute_statistics(graph_from_adjacency_matrix(amat(bnlearn::skeleton(dag))), graph_from_adjacency_matrix(amat(bnlearn::skeleton(aracne_res))))
# e
plot_graph(aracne_adj, real_graph, offset_size=11)
