# sos-pathfinding

R code for:

> Charyyev, P., Zhou, L., & Ceyhan, E. "Navigation algorithms for optimal pathfinding in stochastic obstacle scene problem."

## Overview

This repository contains implementations of six penalty-based navigation
algorithms for the discretized Stochastic Obstacle Scene (SOS) problem, along
with functions for generating synthetic obstacle environments used in the
Monte Carlo experiments reported in the paper.

## Repository Structure

```
sos-pathfinding/
├── README.md
└── R/
    ├── algorithms/
    │   ├── helpers.R       Shared helper functions
    │   ├── RD_alg.R        # Reset Disambiguation (RD)
    │   ├── AP_alg.R        # Additive Penalty (AP)
    │   ├── DT_alg.R        # Distance-to-Target (DT)
    │   ├── ACS_alg.R       # ACS(k)
    │   ├── TACS_alg.R      # Threshold-Adjusted ACS (TACS)
    │   └── MTACS_alg.R     # Multiple-Threshold ACS (MTACS)
    └── data_generation/
        └── generate_environments.R   # Synthetic environment generators
```

## Dependencies

```r
install.packages(c("igraph", "spatial"))
```

## Usage

Each algorithm file sources `helpers.R` automatically. Source the desired
algorithm file and call its main function with an obstacle data frame and
radius:

```r
source("R/algorithms/ACS_alg.R")   # also loads helpers.R

# obs_info: data frame with columns x, y, cost, prob, status
result <- ACS_Alg(obs_info, r = 4.5, acs_k = 3)

result$Length_total   # total Euclidean path length
result$Cost_total     # total disambiguation cost
result$Optimal_path   # vertex sequence of traversed path
```

## Generating Environments

```r
source("R/data_generation/generate_environments.R")

# Single environment (50 false obstacles, cost = 5, lambda = 2)
env <- obs_gen_clutter(n_points = 50, cost = 5, lambda = 2,
                       x_start = 10, x_end = 90,
                       y_start = 10, y_end = 90,
                       inhibition_radius = 9)
```

All obstacle environments are synthetically generated from specified
parametric distributions: obstacle centers follow a Strauss spatial point
process over the window, for example, [10,90] x [10,90], and sensor 
probability marks follow Beta distributions parameterized by the sensor 
accuracy parameter lambda.
