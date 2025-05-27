library(igraph)
library(ggplot2)
library(reshape2)
library(viridis)
library(plotly)
library(scales)
library(ggrepel)

edge_list = c(1, 2,   # e1
               1, 7,   # e2
               2, 3,   # e3
               3, 7,   # e4
               3, 4,   # e5
               4, 6,   # e6
               4, 5,   # e7
               5, 6,   # e8
               6, 8,   # e9
               8, 9,  # e10
               2, 7,  # e11
               7, 8,   # e12
               7, 9,  # e13
               2, 10,  # e14
               2, 11, # e15
               4, 12, # e16
               4, 13,  # e17
               5, 14, # e18
               8, 16,  # e19
               1, 18,  # e20
               9, 18,  # e21
               1, 9,   # e22
               8, 15,  # e23
               8, 17)  # e24

# Create the graph
g = graph(edges = edge_list, directed = FALSE)

E_count = ecount(g)

x_t_edges_mat = t(read.csv("infection_v18.csv", header = FALSE))

T_max = dim(x_t_edges_mat)[2] 

pc_matrix = matrix(0, nrow = ecount(g), ncol = T_max)

# Compute edge percolation centrality over time
for (t in seq(1,T_max,10)) {
  x_t = x_t_edges_mat[, t]
  source("~/Downloads/Edge Percolation/edge_percolation_centrality.R")
  pc_matrix[, t] = edge_percolation_centrality(g, x_t)
}

time_steps = seq(1, T_max, 10)
non_zero_cols = colSums(pc_matrix != 0) > 0
nonzero_pc_matrix = pc_matrix[, non_zero_cols, drop = FALSE]
non_zero_time_steps = time_steps[non_zero_cols]

avg_epc = rowMeans(nonzero_pc_matrix)
rank_epc = order(avg_epc, decreasing = TRUE) # edge rank: 21 10  9 20  4  5 23 19 24 13 12  8  1  6  2 16 22 17 11  3 15 18 14  7

pc_df = melt(nonzero_pc_matrix)
colnames(pc_df) = c("Edge", "TimeIndex", "EPC")

# Replace TimeIndex with actual time values
pc_df$Time = non_zero_time_steps[pc_df$TimeIndex]

# Update Edge to factor with labels e1, e2, ...
pc_df$Edge = factor(pc_df$Edge, levels = unique(pc_df$Edge),
                    labels = paste0("e", seq_along(unique(pc_df$Edge))))

# Plot
ggplot(pc_df, aes(x = Time, y = Edge, fill = EPC)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma") +
  labs(title = "", x = "Time (in Steps of 10)", y = "Edge") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))




## use code below if calculating at every time step
# for (t in 1:T_max) {
#   x_t = x_t_edges_mat[, t]
#   source("~/Downloads/Edge Percolation/edge_percolation_centrality.R")
#   pc_matrix[, t] = edge_percolation_centrality(g, x_t)
# }
# 
# pc_df = melt(nonzero_pc_matrix)
# colnames(pc_df) = c("Edge", "Time", "EPC")
# 
# 
# # Update Edge to factor with labels e1, e2, ...
# pc_df$Edge = factor(pc_df$Edge, levels = unique(pc_df$Edge),
#                     labels = paste0("e", seq_along(unique(pc_df$Edge))))
# 
# ggplot(pc_df, aes(x = Time, y = Edge, fill = EPC)) +
#   geom_tile() +
#   scale_fill_viridis_c(option = "plasma") +
#   labs(title = "",
#        x = "Time Step", y = "Edge") +
#   theme_minimal() +
#   theme(axis.text.y = element_text(size = 8))