library(igraph)

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
  `1` = 1,
  `2` = 1,
  `3` = 1,
  `4` = 1,
  `5` = 1,
  `6` = 1,
  `7` = 1,
  `8` = 1
)

percolation_centrality = function(g, x_t) {
  N = vcount(g)
  pc = numeric(N) 
  
  total_x_t = sum(x_t)
  
  for (v in V(g)) {
    v_idx = as.integer(v)
    sum_term = 0
    
    for (s in V(g)) {
      s_idx = as.integer(s)
      if (s_idx == v_idx) next
      for (r in V(g)) {
        r_idx = as.integer(r)
        if (r_idx == v_idx || r_idx == s_idx) next
        
        # Calculate all shortest paths from s to r
        all_paths = all_shortest_paths(g, from = s, to = r, mode = "all")$res
        sigma_sr = length(all_paths)
        
        if (sigma_sr == 0) next
        
        # Count how many shortest paths from s to r pass through v
        sigma_sr_v = sum(sapply(all_paths, function(path) v_idx %in% path[-c(1, length(path))]))
        
        if ((total_x_t - x_t[v_idx]) != 0) {
          weight = (x_t[s_idx] / (total_x_t - x_t[v_idx]))
          term = (sigma_sr_v / sigma_sr) * weight
          sum_term = sum_term + term
        }
      }
    }
    pc[v_idx] = sum_term / (N - 2)
  }
  
  return(pc)
}

pc = percolation_centrality(g, x)
print(round(pc, 3))

betw = betweenness(g, directed = FALSE, normalized = TRUE)
print(round(betw, 3))


## COMMENT: To show numerical equivalency for directed graphs, change 
## directed = TRUE in lines 14 and 70 and mode = "out" in line 46