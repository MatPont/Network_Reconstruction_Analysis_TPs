library(bnlearn)
library(igraph)
library(pcalg)

data(insurance)

compute_statistics <- function(adj_1, adj_2, do_print=TRUE){
  tp <- sum(adj_1 == adj_2 & adj_1 != 0 & adj_2 != 0)
  fp <- sum(adj_1 != adj_2 & adj_2 != 0)
  fn <- sum(adj_1 != adj_2 & adj_1 != 0)
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
real_adj <- amat(cpdag(dag))
# e
real_graph <- graph_from_adjacency_matrix(real_adj)
plot(real_graph, edge.arrow.size=0.1)


#######
# 2
#######
# b
hc_res <- hc(insurance)
# c
hc_adj <- amat(hc_res)
# d
hc_graph <- graph_from_adjacency_matrix(hc_adj)
plot(hc_graph, edge.arrow.size=0.1, layout=layout_with_lgl)
# e
hc_stats <- compute_statistics(real_adj, hc_adj)
# f

# g
   


#######
# 3
#######
# a

# b
pc_adj <- amat(pc_res)
# c
pc_graph <- graph_from_adjacency_matrix(pc_adj)
plot(pc_graph, edge.arrow.size=0.1, layout=layout_with_lgl)
# d

# e



#######
# 4
#######
# a
aracne_res <- aracne(insurance)
# b
aracne_adj <- amat(aracne_res)
# c
aracne_graph <- graph_from_adjacency_matrix(aracne_adj)
plot(aracne_graph, edge.arrow.size=0.1, layout=layout_with_lgl)
# d

# e
