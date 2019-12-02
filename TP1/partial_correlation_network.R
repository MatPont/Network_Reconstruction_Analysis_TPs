library(miic)
library(ggraph)
library(igraph)
library(Hmisc)

data(hematoData)

data <- scale(as.matrix(hematoData))

hematopoietic <- c("Runx1", "Ikaros", "Myb", "Cbfa2t3h", "Gata1", "Mitf", "Nfe2", "Gfi1b", "Sfpi1", "Gfi1")
endothelial <- c("Erg", "Sox17", "Notch1", "Tbx3", "Tbx20", "Sox7", "HoxB4")
unclassified <- c("HoxB2", "HoxD8")

lambda = 0
lambda = diag(lambda,dim(data)[2],dim(data)[2])

data_cov=cov(data) + lambda
res_corr_filter = solve(data_cov) * (-1)

# Filter edges with a threshold
# thresh = 0.05
# res_corr_filter[abs(res_corr_filter) < thresh] <- 0

# Filter edges by keeping highest values
my_order <- t(apply(res_corr_filter, MARGIN=1, FUN=rank))


to_keep <- 4
res_corr_filter[my_order < (ncol(res_corr_filter) - to_keep) & my_order > to_keep-1 | (my_order < to_keep-1 & res_corr_filter > 0)] <- 0


# to_keep <- 5
# res_corr_filter[my_order < (ncol(res_corr_filter) - to_keep)] <- 0
# #res_corr_filter[my_order > to_keep] <- 0
# to_keep <- 3
# res_corr_filter[my_order < (ncol(res_corr_filter) - to_keep) & (my_order > to_keep & res_corr_filter > 0)] <- 0
# res_corr_filter[my_order < (ncol(res_corr_filter) - to_keep) & my_order > to_keep] <- 0

# Make graph
graph <- graph_from_adjacency_matrix(res_corr_filter, diag = FALSE, weight = TRUE, mode = "undirected")

# Parametrize
V(graph)$color <- ifelse(attr(V(graph), "names") %in% hematopoietic, "red", 
                         ifelse(attr(V(graph), "names") %in% endothelial , "violet", 
                                ifelse(attr(V(graph), "names") %in% unclassified , "gray", "deepskyblue")))
E(graph)$color <- ifelse(E(graph)$weight > 0, "blue","red")
E(graph)$width <- abs(E(graph)$weight) * 3
#deg <- degree(graph, mode = "all")
#V(graph)$size <- deg/std(deg)*20


# Layout
my_layout <- layout_with_dh
my_layout <- layout_with_gem
my_layout <- layout_with_graphopt
my_layout <- layout_with_lgl


# Plot
plot(graph, layout=my_layout, edge.arrow.size=0.25)
plot(graph)

# Choose layout
#layout_list_names <- c("layout_as_bipartite", "layout_as_star", "layout_as_tree", "layout_in_circle", "layout_nicely", "layout_on_grid", "layout_on_sphere", "layout_randomly", "layout_with_dh", "layout_with_fr", "layout_with_gem", "layout_with_graphopt", "layout_with_kk", "layout_with_lgl", "layout_with_mds", "layout_with_sugiyama")
#layout_list <- c(layout_as_bipartite, layout_as_star, layout_as_tree, layout_in_circle, layout_nicely, layout_on_grid, layout_on_sphere, layout_randomly, layout_with_dh, layout_with_fr, layout_with_gem, layout_with_graphopt, layout_with_kk, layout_with_lgl, layout_with_mds, layout_with_sugiyama)
# layout_list_names <- c("layout_with_dh", "layout_with_gem", "layout_with_graphopt", "layout_with_lgl")
# layout_list <- c(layout_with_dh, layout_with_gem, layout_with_graphopt, layout_with_lgl)
# i <- 1
# for(my_layout in layout_list){
#   print(layout_list_names[i])
#   try(plot(graph, layout=my_layout, edge.arrow.size=0.5))
#   readline("Wait")
#   i <- i + 1
# }

