library(bnlearn)
library(igraph)
library(pcalg)

data(insurance)

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
real_adjacency_matrix <- amat(dag)
# e
real_graph <- graph_from_adjacency_matrix(real_adjacency_matrix)
plot(real_graph, edge.arrow.size=0.25)



#######
# 2
#######
# b
res_hc <- hc(insurance)
# c
adjacency_matrix <- amat(dag)
# d
graph <- graph_from_adjacency_matrix(adjacency_matrix)
plot(graph, edge.arrow.size=0.25, layout=layout_with_lgl)
# e
     