edge_percolation_centrality = function(g, x_t_edges) {
  N = vcount(g)
  E_count = ecount(g)
  pc_edge = numeric(E_count)
  
  edge_list = ends(g, E(g), names = FALSE)
  
  for (h in seq_len(E_count)) {
    sum_term = 0
    
    for (s in V(g)) {
      for (r in V(g)) {
        if (s == r) next
        
        all_paths = all_shortest_paths(g, from = s, to = r, mode = "all")$res # identify shortest path(s)
        sigma_sr = length(all_paths) # number of shortest path(s)
        if (sigma_sr == 0) next
        
        sigma_sr_h = 0
        source_edge_id = NA
        
        for (path in all_paths) {
          path = as.integer(path)
          if (length(path) < 2) next # if shortest path consists of 2 vertices (or 1 edge), skip
          
          edge_h_nodes = edge_list[h, ]
          has_h = any(sapply(seq_along(path)[-length(path)], function(i) {
            all(sort(c(path[i], path[i+1])) == sort(edge_h_nodes))
          }))
          
          if (has_h) {
            sigma_sr_h = sigma_sr_h + 1 # if shortest path contains edge h, increment counter
            
            if (is.na(source_edge_id)) {
              source_pair = sort(path[1:2])
              source_edge_id = which(apply(edge_list, 1, function(e) all(sort(e) == source_pair)))
            }
          }
        }
        
        if (sigma_sr_h > 0 && length(source_edge_id) == 1) {
          denom = sum(x_t_edges) - x_t_edges[h]
          if (denom != 0) {
            weight = x_t_edges[source_edge_id] / denom
            term = (sigma_sr_h / sigma_sr) * weight
            sum_term = sum_term + term
          }
        }
      }
    }
    
    pc_edge[h] = sum_term * ((E_count - 1) / (N * (N - 1)))
  }
  
  return(pc_edge)
}