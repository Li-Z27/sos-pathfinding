# =============================================================================
# helpers.R
# Shared helper functions for SOS pathfinding algorithms.
# =============================================================================

library(igraph)

# Function: Graph_Discretized
# Constructs an 8-adjacency weighted lattice graph on an (x+1) x (y+1) grid.
# -----------------------------------------------------------------------------
Graph_Discretized <- function(x, y) {
  # INPUT:
  #   x : grid width parameter
  #   y : grid height parameter
  # OUTPUT: 
  #   G : igraph weighted undirected graph
  n <- (x + 1) * (y + 1)
  temp <- graph.empty(n, directed = FALSE)
  G <- set.vertex.attribute(temp, "name", value = c(1:n))
  A <- get.adjacency(G, sparse = FALSE)

  # Interior vertices
  for (j in 1:(y - 1)) {
    for (i in 1:(x - 1)) {
      A[1 + i + j * (x + 1), i + j * (x + 1)] <- 1        # left
      A[1 + i + j * (x + 1), 2 + i + j * (x + 1)] <- 1        # right
      A[1 + i + j * (x + 1), 1 + i + (j - 1) * (x + 1)] <- 1    # down
      A[1 + i + j * (x + 1), 1 + i + (j + 1) * (x + 1)] <- 1    # up
      
      A[1 + i + j * (x + 1), i + (j + 1) * (x + 1)] <- sqrt(2) # left-up
      A[1 + i + j * (x + 1), 2 + i + (j + 1) * (x + 1)] <- sqrt(2) # right-up
      A[1 + i + j * (x + 1), i + (j - 1) * (x + 1)] <- sqrt(2) # left-down
      A[1 + i + j * (x + 1), 2 + i + (j - 1) * (x + 1)] <- sqrt(2) # right-down
    }
  }

  # Bottom boundary (j = 0)
  for (i in 1:(x - 1)) {
    A[1 + i, i] <- 1
    A[1 + i, 2 + i] <- 1
    A[1 + i, 1 + i + (x + 1)] <- 1
    A[1 + i, i + (x + 1)] <- sqrt(2)
    A[1 + i, 2 + i + (x + 1)] <- sqrt(2)
  }

  # Top boundary (j = y)
  for (i in 1:(x - 1)) {
    A[1 + i + y * (x + 1), i + y * (x + 1)] <- 1
    A[1 + i + y * (x + 1), 2 + i + y * (x + 1)] <- 1
    A[1 + i + y * (x + 1), 1 + i + (y - 1) * (x + 1)] <- 1
    A[1 + i + y * (x + 1), i + (y - 1) * (x + 1)] <- sqrt(2)
    A[1 + i + y * (x + 1), 2 + i + (y - 1) * (x + 1)] <- sqrt(2)
  }

  # Left boundary (i = 0)
  for (j in 1:(y - 1)) {
    A[1 + j * (x + 1), 1 + (j - 1) * (x + 1)] <- 1
    A[1 + j * (x + 1), 2 + j * (x + 1)] <- 1
    A[1 + j * (x + 1), 1 + (j + 1) * (x + 1)] <- 1
    A[1 + j * (x + 1), 2 + (j - 1) * (x + 1)] <- sqrt(2)
    A[1 + j * (x + 1), 2 + (j + 1) * (x + 1)] <- sqrt(2)
  }

  # Right boundary (i = x)
  for (j in 1:(y - 1)) {
    A[1 + x + j * (x + 1), 1 + x + (j - 1) * (x + 1)] <- 1
    A[1 + x + j * (x + 1), 1 + x + (j + 1) * (x + 1)] <- 1
    A[1 + x + j * (x + 1), x + j * (x + 1)] <- 1
    A[1 + x + j * (x + 1), x + (j - 1) * (x + 1)] <- sqrt(2)
    A[1 + x + j * (x + 1), x + (j + 1) * (x + 1)] <- sqrt(2)
  }

  # Four corners
  A[1, 2] <- 1; A[1, 2 + x] <- 1; A[1, 3 + x] <- sqrt(2)
  A[1 + y * (x + 1), 1 + (y - 1) * (x + 1)] <- 1
  A[1 + y * (x + 1), 2 + y * (x + 1)] <- 1
  A[1 + y * (x + 1), 2 + (y - 1) * (x + 1)] <- sqrt(2)
  A[1 + x, x] <- 1; A[1 + x, 2 * (1 + x)] <- 1; A[1 + x, 2 * (1 + x) - 1] <- sqrt(2)
  A[(1 + x) * (1 + y), (1 + x) * y] <- 1
  A[(1 + x) * (1 + y), (1 + x) * (1 + y) - 1]  <- 1
  A[(1 + x) * (1 + y), (1 + x) * y - 1] <- sqrt(2)

  G <- graph.adjacency(A, mode = "undirected", weighted = TRUE)
  return(G)
}

# Function: Intersect_Obs
# Finds all graph edges that intersect the boundary of an obstacle.
# -----------------------------------------------------------------------------
Intersect_Obs <- function(c, r, x, y) {
  # INPUT:
  #   c : vector (cx, cy) — obstacle center coordinates
  #   r : obstacle radius
  #   x : grid width parameter
  #   y : grid height parameter
  # OUTPUT: 
  #  el : 2-column matrix of edge endpoint vertex indices (edge list)
  cx <- c[1]; cy <- c[2]
  cx1 <- floor(cx); cy1 <- floor(cy)
  r1 <- ceiling(r) + 2
  coor_info <- Lattice_Vertices(x, y)

  X_temp <- matrix(
    c(rep(-r1:r1, times = length(seq(-r1, r1, by = 1))),
      rep(-r1:r1, each  = length(seq(-r1, r1, by = 1)))),
    ncol = 2
  )

  # Vertices just outside the obstacle boundary
  Intersect_temp1 <- function(vector_ij) {
    x_i <- vector_ij[1]; x_j <- vector_ij[2]
    d2 <- (cx - (cx1 + x_i))^2 + (cy - (cy1 + x_j))^2
    if (r^2 < d2 && d2 <= (r + sqrt(2))^2){
      return(c(cx1 + x_i, cy1 + x_j))
    }
  }
  temp  <- as.numeric(unlist(apply(X_temp, 1, Intersect_temp1)))
  case2 <- matrix(temp, ncol = 2, byrow = TRUE)

  # Vertices just inside and on the obstacle boundary
  Intersect_temp2 <- function(vector_ij) {
    x_i <- vector_ij[1]; x_j <- vector_ij[2]
    d2 <- (cx - (cx1 + x_i))^2 + (cy - (cy1 + x_j))^2
    if ((r - sqrt(2))^2 < d2 && d2 <= r^2){
      return(c(cx1 + x_i, cy1 + x_j))
    }
  }
  temp  <- as.numeric(unlist(apply(X_temp, 1, Intersect_temp2)))
  case3 <- matrix(temp, ncol = 2, byrow = TRUE)
  
  # Form edge list matrix
  el <- matrix(0, ncol = 2)
  n2 <- nrow(case2); n3 <- nrow(case3)
  if (n2 != 0 && n3 != 0) {
    X_temp2 <- matrix(
      c(rep(1:n3, times = n2),
        rep(1:n2, each  = n3)),
      ncol = 2
    )
    Intersect_temp3 <- function(vector_ij) {
      x_i <- vector_ij[1]; x_j <- vector_ij[2]
      d <- Dist_Euclidean(case2[x_j, ], case3[x_i, ])
      if (d == 1 || d == sqrt(2)) {
        e1 <- which(coor_info[, 1] == case2[x_j, 1] & coor_info[, 2] == case2[x_j, 2])
        e2 <- which(coor_info[, 1] == case3[x_i, 1] & coor_info[, 2] == case3[x_i, 2])
        return(sort(c(e1, e2), decreasing = FALSE))
      }
    }
    el_vector <- as.numeric(unlist(apply(X_temp2, 1, Intersect_temp3)))
    el <- matrix(el_vector, ncol = 2, byrow = TRUE)
  }
  return(el)
}

# Function: Dist_Euclidean
# Computes Euclidean distance between two points.
# -----------------------------------------------------------------------------
Dist_Euclidean <- function(point1, point2) {
  sqrt((point1[1] - point2[1])^2 + (point1[2] - point2[2])^2)
}


# Lattice_Vertices
# Generates (x,y) coordinates for all vertices in the (x+1)x(y+1) lattice.
# -----------------------------------------------------------------------------
Lattice_Vertices <- function(x, y) {
  # INPUT:
  #   x : grid width parameter
  #   y : grid height parameter
  # OUTPUT: 
  #   coords : matrix with columns (x_coord, y_coord), one row per vertex
  coords <- matrix(nrow = (x + 1) * (y + 1), ncol = 2)
  coords[, 1] <- rep(0:x, y + 1)
  coords[, 2] <- rep(0:y, each = x + 1)
  return(coords)
}

# Index_Coordinates
# Returns the vertex index for a given (x, y) grid coordinate.
# -----------------------------------------------------------------------------
Index_Coordinates <- function(m, x, y) {
  1 + m[1] + m[2] * (x + 1)
}
