library(igraph)
library(microbenchmark)

# Define edges from Figure 1 (https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0053095&type=printable)
edges = matrix(c(
  1, 4,
  2, 4,
  3, 4,
  4, 5,
  5, 6,
  6, 7,
  6, 8
), byrow = TRUE, ncol = 2)

g = graph_from_edgelist(edges, directed = FALSE)
coords = matrix(c(-1.4657215, -1.8178889, 
           -2.8319484, -0.9782229, 
           -2.4833739,  0.5869833, 
           -1.4642752, -0.3993414, 
           0.3256479,  0.3607142, 
           2.0340507,  1.0822883, 
           3.3885511,  0.7233162, 
           2.7286221,  2.3082103), 
           byrow = TRUE, ncol = 2)


# Define percolation states for Fig a
xa = c(
  `1` = 0.2, # 1 - 4
  `2` = 0.1, # 2 - 4
  `3` = 0.2, # 3 - 4
  `4` = 0.1, # 4 - 5
  `5` = 0.1, # 5 - 6
  `6` = 0.5, # 6 - 7
  `7` = 0.5 # 6 - 8
)

# Define percolation states for Fig b
xb = c(
  `1` = 0.5, # 1 - 4
  `2` = 0.1, # 2 - 4
  `3` = 0.5, # 3 - 4
  `4` = 0.1, # 4 - 5
  `5` = 0.1, # 5 - 6
  `6` = 0.2, # 6 - 7
  `7` = 0.2 # 6 - 8
)

x = xb

source("~/Downloads/Edge Percolation/edge_percolation_centrality.R")
pc_edge = edge_percolation_centrality(g, x)

E(g)$width = pc_edge
E(g)$label = x

plot(g,
     layout = coords,
     edge.width = E(g)$width*20,
     edge.label = E(g)$label,
     edge.label.cex = 2,
     edge.label.color = "black",
     vertex.label.color = "black",
     vertex.color = "orange",
)

print(round(pc_edge, 3))