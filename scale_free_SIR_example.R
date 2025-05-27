library(igraph)

set.seed(42)

n_nodes = 30
p = 0.2
T_max = 50

g = sample_pa(n = n_nodes, directed = FALSE)

deg = degree(g)
peripheral_nodes = which(deg == 1) # identify peripheral edges
peripheral_edges = E(g)[.from(peripheral_nodes) | .to(peripheral_nodes)]

# initialize percolation
E(g)$percolated = 0
initial_edge = sample(peripheral_edges, 1)
E(g)[initial_edge]$percolated = 1

percolation_history = matrix(0, nrow = ecount(g), ncol = T_max)
percolation_history[, 1] = E(g)$percolated

percolation_centrality = matrix(NA, nrow = ecount(g), ncol = T_max)
betweenness_centrality = matrix(NA, nrow = ecount(g), ncol = T_max)

percolated_edge_count = numeric(T_max)
percolated_edge_count[1] = sum(E(g)$percolated)

source("~/Downloads/Edge Percolation/edge_percolation_centrality.R")
percolation_centrality[, 1] = edge_percolation_centrality(g, E(g)$percolated)
betweenness_centrality[, 1] = edge_betweenness(g, directed = FALSE)

# Simulation
for (t in 2:T_max) {
  percolated_prev = percolation_history[, t - 1]
  percolated_new = percolated_prev
  percolated_edges = which(percolated_prev == 1)
  
  for (e_idx in percolated_edges) {
    nodes = ends(g, E(g)[e_idx])
    adjacent_edges = c(incident(g, nodes[1], mode = "all"), incident(g, nodes[2], mode = "all"))
    adjacent_edges = setdiff(adjacent_edges, percolated_edges)
    
    for (adj_e in adjacent_edges) {
      if (percolated_new[adj_e] == 0 && runif(1) < p) {
        percolated_new[adj_e] = 1
      }
    }
  }
  
  percolation_history[, t] = percolated_new
  E(g)$percolated = percolated_new
  percolated_edge_count[t] = sum(percolated_new)
  
  source("~/Downloads/Edge Percolation/edge_percolation_centrality.R")
  percolation_centrality[, t] = edge_percolation_centrality(g, percolated_new)
  betweenness_centrality[, t] = edge_betweenness(g, directed = FALSE)
}

betweenness_centrality = (2 / (n_nodes * (n_nodes - 1))) * betweenness_centrality

ratio = colMeans(percolation_centrality) / colMeans(betweenness_centrality)


layout_fixed = layout_with_fr(g)

par(mfrow = c(2, 5), mar = c(1, 1, 2, 1))
for (t in 1:T_max) {
  edge_colors = ifelse(percolation_history[, t] == 1, "red", "gray")
  plot(g,
       layout = layout_fixed,
       edge.color = edge_colors,
       edge.width = 2,
       vertex.size = 5,
       vertex.label = NA,
       main = paste("Time Step", t))
}


par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
plot(1:T_max, percolated_edge_count,
     type = "b", col = "firebrick", lwd = 2, pch = 19,
     xlab = "Time Step", ylab = "Number of Infected Edges",
     main = "",
     cex.axis = 1.5,      
     cex.lab = 1.5,       
     cex = 1.5)       
grid()

plot(1:T_max, ratio,
     type = "b", col = "darkblue", lwd = 2, pch = 19,
     xlab = "Time Step", ylab = "Mean EPC / Mean EBC",
     main = "",
     cex.axis = 1.5,      
     cex.lab = 1.5,       
     cex = 1.5)  
grid()


par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
plot(1:ecount(g), rowMeans(percolation_centrality),
     type = "b", col = "red", lwd = 2, pch = 19,
     xlab = "Edge", ylab = "Mean EPC",
     main = "",
     cex.axis = 1.5,      
     cex.lab = 1.5,       
     cex = 1.5)  
grid()

plot(1:ecount(g), rowMeans(betweenness_centrality),
     type = "b", col = "darkblue", lwd = 2, pch = 19,
     xlab = "Edge", ylab = "Mean EBC",
     main = "",
     cex.axis = 1.5,      
     cex.lab = 1.5,       
     cex = 1.5)  
grid()

par(mfrow = c(1, 1))
plot(g, layout = layout_fixed)
