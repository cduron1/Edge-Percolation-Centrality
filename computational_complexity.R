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


# Define percolation states
x = c(
  `1` = 0.3,
  `2` = 0.5,
  `3` = 0.5,
  `4` = 0.2,
  `5` = 0.3,
  `6` = 0.2,
  `7` = 0.1
)

source("~/Downloads/Edge Percolation/edge_percolation_centrality.R")
pc_edge = edge_percolation_centrality(g, x)
print(round(pc_edge, 3))

# ## COMMENT: When all edges are percolate (i.e., in the same state), then
# ## pc_edge reduces to edge_betw
# edge_betw = edge_betweenness(g, directed = FALSE)
# edge_betw = ( 2 / (vcount(g)*(vcount(g) - 1) ) )*edge_betw
# print(round(edge_betw, 3))
# 
# ## COMMENT: To show numerical equivalency for directed graphs, change 
# ## directed = TRUE in lines 14 and 72, mode = "out" in line 42, and 2->1 in line 73


## Computational complexity numerical verification
sizes = c(10, 20, 30, 40, 50, 100)
times = numeric(length(sizes))

for (i in seq_along(sizes)) {
  N = sizes[i]
  cat("Running N =", N, "\n")
  
  # Generate scale-free network
  g = sample_pa(N, directed = FALSE)
  
  # Assign random percolation states to edges (between 0 and 1)
  x_t_edges = runif(ecount(g), min = 0, max = 1)
  
  # Time the percolation centrality computation
  times[i] = system.time({
    source("~/Downloads/Edge Percolation/edge_percolation_centrality.R")
    pc_edge = edge_percolation_centrality(g, x_t_edges)
  })[["elapsed"]]
}

# Plot log-log
plot(sizes, times, type = "b", log = "xy", 
     xlab = "Number of Vertices", 
     ylab = "Elapsed Time (seconds)", 
     main = "", 
     pch = 19, 
     cex.axis = 1.5,      
     cex.lab = 1.5,       
     cex = 1.5)           

# Fit and report slope of log-log line
log_model = lm(log(times) ~ log(sizes))
abline(log_model, col = "red", lty = 2)
cat("Estimated complexity O(n^", round(coef(log_model)[2], 2), ")\n")
