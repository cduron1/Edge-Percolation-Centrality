library(igraph)

# Reconstructed edge list
edges = c(
  "m013", "f201",
  "f201", "m012",
  "m012", "f202",
  "m012", "f009",
  "f009", "m201",
  "m201", "f307",
  "f009", "m017",
  "f033", "m017",
  "m017", "f336",
  "m017", "f015",
  "m017", "f034",
  "f015", "m026",
  "m026", "f010",
  "f009", "m202",
  "m202", "f011",
  "f011", "m204",
  "f011", "m203",
  "m203", "f017",
  "f011", "m526",
  "f900", "m526",
  "m526", "f514",
  "m526", "f019",
  "m017", "f019",
  "f019", "m206",
  "m206", "f022",
  "f022", "m025",
  "f022", "m207",
  "m207", "f012",
  "f011", "m551",
  "m551", "f022",
  "m209", "f022",
  "m208", "f022",
  "f022", "m016",
  "m016", "f014",
  "m016", "f533",
  "m016", "f038",
  "m016", "f020",
  "f020", "m019",
  "f020", "m210",
  "m017", "f000"
)

g = graph(edges, directed = FALSE)

# Convert to character to preserve original vertex IDs as labels
edges_chr = as.character(edges)

p = 0.2
T_max = 40

deg = degree(g)
peripheral_nodes = which(deg == 1) # identify peripheral edges
peripheral_edges = E(g)[.from(peripheral_nodes) | .to(peripheral_nodes)]

# initialize percolation
E(g)$percolated = 0
initial_edge = peripheral_edges[1] # select edge as VPC did in their example; or use --> #initial_edge = sample(peripheral_edges, 1)
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

betweenness_centrality = (2 / (vcount(g) * (vcount(g) - 1))) * betweenness_centrality

ratio = colMeans(percolation_centrality) / colMeans(betweenness_centrality)


layout_fixed = matrix(c(
  -7.2474212, 8.8042300,
  -5.6834361, 7.7551449,
  -3.4344284, 6.3985421,
  -4.3569997, 7.7544673,
  -0.7977739, 4.3494752,
  -3.0354277, 4.0890278,
  -4.7824493, 3.8583603,
  1.9098215, 5.1863500,
  0.7128072, 6.2285130,
  3.4638951, 5.8051707,
  2.6240372, 8.2935224,
  1.5695060, 6.9575696,
  3.1685181, 10.7639744,
  3.5702581, 12.5694883,
  1.4239402, 1.7711740,
  3.9665346, -0.7092591,
  5.0211109, -1.5838874,
  6.5172158, -1.0695579,
  8.4324572, -1.2773334,
  4.9554542, 1.1109990,
  6.7699797, 1.4540620,
  6.1610545, 2.3925819,
  2.8655421, 1.2223070,
  2.8767778, -2.8194027,
  3.2280269, -6.7679139,
  1.7990118, -7.0459270,
  1.2878825, -8.4601822,
  -0.1467148, -9.6568255,
  3.6798992, -3.8722792,
  4.6098829, -7.4128939,
  3.2125426, -8.2240645,
  4.4611944, -10.7092451,
  3.6625194, -12.2234669,
  5.2490490, -12.1328280,
  6.1123317, -11.3312474,
  5.0691260, -13.7318091,
  4.8256769, -15.6149000,
  6.0575577, -15.3025644,
  2.8088197, 6.6586317
), ncol = 2, byrow = TRUE)


par(mfrow = c(1,1))
#for (t in c(1, 10, 20, 30, 40)) { #
for (t in 1:T_max ){
  edge_colors = ifelse(percolation_history[, t] == 1, "red", "gray")
  plot(g,
       layout = layout_fixed,
       edge.color = edge_colors,
       edge.width = 3,
       vertex.size = 5,
       vertex.label = NA)
}


par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
plot(1:T_max, percolated_edge_count,
     type = "b", col = "black", lwd = 2, pch = 19,
     xlab = "Time Step", ylab = "Number of Infected Edges",
     main = "",
     cex.axis = 1.5,      
     cex.lab = 1.5,       
     cex = 1.5)       
grid()

plot(1:T_max, ratio,
     type = "b", col = "black", lwd = 2, pch = 19,
     xlab = "Time Step", ylab = "Mean EPC / Mean EBC",
     main = "",
     cex.axis = 1.5,      
     cex.lab = 1.5,       
     cex = 1.5)  
grid()


par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
plot(1:ecount(g), rowMeans(percolation_centrality),
     type = "b", col = "black", lwd = 2, pch = 19,
     xlab = "Edge", ylab = "Mean EPC",
     main = "",
     cex.axis = 1.5,      
     cex.lab = 1.5,       
     cex = 1.5)  
grid()

plot(1:ecount(g), rowMeans(betweenness_centrality),
     type = "b", col = "black", lwd = 2, pch = 19,
     xlab = "Edge", ylab = "Mean EBC",
     main = "",
     cex.axis = 1.5,      
     cex.lab = 1.5,       
     cex = 1.5)  
grid()


par(mfrow = c(1, 1))
plot(g,
     edge.width = betweenness_centrality[1:ecount(g),1]*10,  # Scale edge width based on centrality
     vertex.label.cex = 0.5,
     vertex.size = 9,     
     edge.color = "black",  # Color the edges
     vertex.color = "orange",
     layout = layout_fixed)  # Color the vertices


par(mfrow = c(1, 1))
plot(g,
     edge.width = percolation_centrality[1:ecount(g),1]*10,  # Scale edge width based on centrality
     vertex.label.cex = 0.5,
     vertex.size = 9,     
     edge.color = "black",  # Color the edges
     vertex.color = "orange",
     layout = layout_fixed)  # Color the vertices

par(mfrow = c(1, 1))
plot(g,
     edge.width = percolation_centrality[1:ecount(g),11]*10,  # Scale edge width based on centrality
     vertex.label.cex = 0.5,
     vertex.size = 9,     
     edge.color = "black",  # Color the edges
     vertex.color = "orange",
     layout = layout_fixed)  # Color the vertices

par(mfrow = c(1, 1))
plot(g,
     edge.width = percolation_centrality[1:ecount(g),37]*10,  # Scale edge width based on centrality
     vertex.label.cex = 0.5,
     vertex.size = 9,     
     edge.color = "black",  # Color the edges
     vertex.color = "orange",
     layout = layout_fixed)  # Color the vertices


# Rankings based on percolation centrality
percolation_ranking = order(rowMeans(percolation_centrality), decreasing = TRUE)

# Rankings based on betweenness centrality
betweenness_ranking = order(rowMeans(betweenness_centrality), decreasing = TRUE)

# Create rank positions
percolation_rank = match(1:40, percolation_ranking)
betweenness_rank = match(1:40, betweenness_ranking)

# Compute Spearman correlation
cor(percolation_rank, betweenness_rank, method = "spearman")

