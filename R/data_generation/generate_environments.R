# =============================================================================
# generate_environments.R
# Synthetic obstacle environment generation for SOS Monte Carlo simulations.
#
# Obstacle centers follow a Strauss spatial point process (spatial package).
# Sensor probability marks follow Beta distributions parameterized by lambda:
#   False obstacles: prob ~ Beta(4 - lambda, 4 + lambda)
#   True obstacles:  prob ~ Beta(4 + lambda, 4 - lambda)
# =============================================================================

library(spatial)

# Function: obs_gen_clutter
# Generates one environment containing false obstacles only.
# -----------------------------------------------------------------------------
obs_gen_clutter <- function(n_points, cost, lambda = 2,
                            x_start, x_end,
                            y_start, y_end,
                            inhibition_radius) {
  # INPUT:
  #   n_points          : number of obstacles
  #   cost              : disambiguation cost per obstacle
  #   lambda            : sensor accuracy parameter (default 2)
  #   x_start, x_end    : x-range of obstacle placement window
  #   y_start, y_end    : y-range of obstacle placement window
  #   inhibition_radius : minimum inter-obstacle distance for Strauss process
  # OUTPUT: data frame with columns x, y, cost, prob, status
  ppregion(xl = x_start, xu = x_end, yl = y_start, yu = y_end)
  pts  <- spatial::Strauss(n = n_points, c = 1, r = inhibition_radius)
  data.frame(
    x      = pts$x,
    y      = pts$y,
    cost   = rep(cost, n_points),
    prob   = rbeta(n_points, 4 - lambda, 4 + lambda),
    status = rep(0, n_points)
  )
}

# Function: obs_gen_obstacle
# Generates one environment containing true obstacles only.
# -----------------------------------------------------------------------------
obs_gen_obstacle <- function(n_points, cost, lambda = 2,
                             x_start, x_end,
                             y_start, y_end,
                             inhibition_radius) {
  # INPUT:
  #   n_points          : number of obstacles
  #   cost              : disambiguation cost per obstacle
  #   lambda            : sensor accuracy parameter (default 2)
  #   x_start, x_end    : x-range of obstacle placement window
  #   y_start, y_end    : y-range of obstacle placement window
  #   inhibition_radius : minimum inter-obstacle distance for Strauss process
  # OUTPUT: data frame with columns x, y, cost, prob, status
  ppregion(xl = x_start, xu = x_end, yl = y_start, yu = y_end)
  pts  <- spatial::Strauss(n = n_points, c = 1, r = inhibition_radius)
  data.frame(
    x      = pts$x,
    y      = pts$y,
    cost   = rep(cost, n_points),
    prob   = rbeta(n_points, 4 + lambda, 4 - lambda),
    status = rep(1, n_points)
  )
}

# Function: obs_gen_mixed
# Generates one environment containing false/true obstacles.
# -----------------------------------------------------------------------------
obs_gen_mixed <- function(n_false, n_true, cost, lambda = 2,
                          x_start, x_end,
                          y_start, y_end,
                          inhibition_radius) {
  clutter  <- obs_gen_clutter(n_false, cost, lambda,
                              x_start, x_end, y_start, y_end,
                              inhibition_radius)
  obstacle <- obs_gen_obstacle(n_true, cost, lambda,
                               x_start, x_end, y_start, y_end,
                               inhibition_radius)
  rbind(clutter, obstacle)
}
