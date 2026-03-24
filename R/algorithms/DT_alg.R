# =============================================================================
# DT_alg.R
# Distance-to-Target (DT) navigation algorithm for the discretized SOS problem.
#
# Penalty function:
#   P_DT(c, p, dt) = c + (dt / (1-p))^(-log(1-p))
# =============================================================================

source("helpers.R")

# Function: Update_graph_DT
# Initializes edge weights on graph g by adding DT penalty terms for each
# obstacle that intersects each edge.
# -----------------------------------------------------------------------------
Update_graph_DT <- function(g, x, y, circle_info, r, dt = 50) {
  # INPUT:
  #   g           : igraph object (base lattice graph)
  #   x, y        : grid dimension parameters
  #   circle_info : data frame with columns x, y, cost, prob, status
  #   r           : obstacle radius
  #   dt          : distance-to-target parameter (default 50)
  # OUTPUT: 
  #   list with G_info (updated igraph) and Int_info (intersection matrix)
  n   <- nrow(circle_info)
  elg <- get.data.frame(g, what = "edges")
  colnames(elg) <- c("From", "To", "Cost")
  int_info <- matrix(0, ncol = n, nrow = nrow(elg))
  
  for (i in 1:n) {
    el <- Intersect_Obs(t(circle_info[i, 1:2]), r, x, y)
    n1 <- nrow(el)
    for (k in 1:n1) el[k, ] <- sort(el[k, ], decreasing = FALSE)
    
    p <- circle_info[i, 4]
    # DT penalty: c + (dt / (1-p))^(-log(1-p))
    penalty <- circle_info[i, 3] + (dt / (1 - p))^(-log(1 - p))
    
    for (j in 1:n1) {
      idx <- which(elg[, 1] == el[j, 1] & elg[, 2] == el[j, 2])
      elg[idx, 3] <- elg[idx, 3] + 0.5 * penalty
      int_info[idx, i] <- 1
    }
  }
  
  list(G_info   = graph.data.frame(elg, directed = FALSE),
       Int_info = int_info)
}

# Function: DT_Alg
# Main DT navigation algorithm. Iteratively finds shortest paths under current
# edge weights, traverse to the first ambiguous obstacle edge, disambiguates,
# and updates the graph accordingly.
# -----------------------------------------------------------------------------
DT_Alg <- function(obs_info, r, x, y, s, t, dt = 50) {
  # INPUT:
  #   obs_info : data frame with columns x, y, cost, prob, status
  #   r        : obstacle radius
  #   x, y     : grid dimension parameters
  #   s, t     : starting and target vertices
  #   dt       : distance-to-target parameter (default 50)
  # OUTPUT: 
  #  list with Optimal_found, Length_total, Cost_total,
  #   Optimal_path, Disambiguate_state

  vertice_list  <- Lattice_Vertices(x, y)
  G_original    <- Graph_Discretized(x, y)
  out_Ginfo     <- Update_graph_DT(G_original, x, y, obs_info, r, dt)
  G_ed          <- out_Ginfo$G_info
  Int_info      <- out_Ginfo$Int_info
  df_edge_ed    <- get.data.frame(G_ed, what = "edges")
  
  length_total <- 0
  cost_total   <- 0
  reach_t      <- FALSE
  path_record  <- s
  D_record     <- c()
  
  while (!reach_t) {
    output <- get.shortest.paths(
      G_ed,
      which(vertex.attributes(G_ed)$name == as.character(s)),
      which(vertex.attributes(G_ed)$name == as.character(t)),
      weights = df_edge_ed$Cost, output = "both", algorithm = "dijkstra"
    )
    V_list <- as.numeric(attributes(output$vpath[[1]])$names)
    
    D_state <- NULL
    for (i in 2:length(V_list)) {
      edge_idx    <- which(df_edge_ed$from == min(V_list[(i-1):i]) &
                             df_edge_ed$to == max(V_list[(i-1):i]))
      edge_length <- Dist_Euclidean(as.numeric(vertice_list[V_list[i-1], ]),
                                    as.numeric(vertice_list[V_list[i], ]))
      
      if (sum(Int_info[edge_idx, ]) != 0) {
        D_state  <- V_list[i - 1]
        D_record <- c(D_record, D_state)
        break
      } else {
        length_total <- length_total + edge_length
        path_record  <- c(path_record, V_list[i])
      }
    }
    
    if (is.null(D_state)) {
      reach_t <- TRUE
      output_final <- list(
        Optimal_found      = TRUE,
        Length_total       = length_total,
        Cost_total         = cost_total,
        Optimal_path       = path_record,
        Disambiguate_state = D_record
      )
    } else {
      s <- D_state
      obs_ind <- which(Int_info[edge_idx, ] == 1)
      
      if (length(obs_ind) > 1) {
        dists   <- sapply(obs_ind, function(oi)
          Dist_Euclidean(as.numeric(vertice_list[D_state, ]),
                         as.numeric(obs_info[oi, 1:2])))
        obs_ind <- obs_ind[which.min(dists)]
      }
      
      cost_total <- cost_total + obs_info[obs_ind, 3]
      p          <- obs_info[obs_ind, 4]
      penalty    <- obs_info[obs_ind, 3] + (dt / (1 - p))^(-log(1 - p))
      
      if (obs_info$status[obs_ind] == 1) {
        df_edge_ed[which(Int_info[, obs_ind] == 1), 3] <- Inf
      } else {
        df_edge_ed[which(Int_info[, obs_ind] == 1), 3] <-
          df_edge_ed[which(Int_info[, obs_ind] == 1), 3] - 0.5 * penalty
        Int_info[which(Int_info[, obs_ind] == 1), obs_ind] <- 0
      }
      G_ed <- graph.data.frame(df_edge_ed, directed = FALSE)
    }
  }
  return(output_final)
}
