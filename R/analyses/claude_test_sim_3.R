# Gibert & Brassil 2014 Figure 3 Replication in R
# Predator-prey system with phenotypic variation

library(deSolve)
library(ggplot2)
library(dplyr)
library(pracma) # for numerical integration

# Define the probability density function for phenotypic distribution
pdist <- function(x, sigma, X) {
  return((1/sqrt(2*pi*sigma)) * exp(-(x-X)^2/(2*sigma)))
}

# Attack rate function
attack <- function(x, a, alpha, tau) {
  return(alpha * exp(-(x-a)^2/(2*tau^2)))
}

# Handling time function  
handling <- function(x, eta, a, nu) {
  return(eta - exp(-(x-a)^2/(2*nu^2)))
}

# Functional response integrand
functional_response_integrand <- function(x, a, alpha, tau, eta, nu, R, sigma, X) {
  att <- attack(x, a, alpha, tau)
  hand <- handling(x, eta, a, nu)
  pd <- pdist(x, sigma, X)
  
  return((att / (1 + att * hand * R)) * pd)
}

# Numerical integration of functional response
integrate_functional_response <- function(a, alpha, tau, eta, nu, R, sigma, X) {
  # Use numerical integration over a reasonable range
  x_range <- c(-10, 10)  # Should cover most of the distribution
  
  result <- integrate(function(x) {
    sapply(x, function(xi) {
      functional_response_integrand(xi, a, alpha, tau, eta, nu, R, sigma, X)
    })
  }, lower = x_range[1], upper = x_range[2])
  
  return(result$value)
}

# Define the ODE system
predator_prey_ode <- function(t, y, parms) {
  R <- y[1]  # Resource
  c <- y[2]  # Consumer
  
  with(parms, {
    # Calculate functional response
    func_resp <- integrate_functional_response(a, alpha, tau, eta, nu, R, sigma, X)
    
    # Resource equation
    dR_dt <- r * R * (1 - R/k) - R * c * func_resp
    
    # Consumer equation  
    dc_dt <- epsilon * R * c * func_resp - beta * c
    
    return(list(c(dR_dt, dc_dt)))
  })
}

# Parameters (from the Mathematica code)
params <- list(
  alpha = 2.0,
  eta = 2.0, 
  r = 0.3,
  beta = 0.1,
  k = 1.0,
  epsilon = 0.5,
  tau = 1.0,
  nu = 1.0,
  a = 0.5,  # This will be varied based on phenotypic distance
  X = 0.0
)

# Create parameter combinations
variance_vals <- seq(0.001, 2.001, by = 0.02)  # 101 values
phenotypic_dist_vals <- seq(0.0, 2.0, by = 0.02)  # 101 values

# Create all combinations
param_grid <- expand.grid(variance = variance_vals, 
                          phenotypic_dist = phenotypic_dist_vals)

# Function to simulate single parameter combination
simulate_single <- function(variance, phenotypic_dist) {
  # Update parameters
  current_params <- params
  current_params$sigma <- variance
  current_params$a <- phenotypic_dist
  
  # Initial conditions
  y0 <- c(R = 1.0, c = 0.1)
  
  # Time sequence
  times <- seq(0, 500, by = 0.5)
  
  # Solve ODE
  tryCatch({
    solution <- ode(y0, times, predator_prey_ode, current_params, method = "lsoda")
    
    # Extract final part of time series to classify dynamics
    final_times <- tail(solution, 200)  # Last 200 time points
    consumer_vals <- final_times[, "c"]
    
    # Classification criteria (adapted from Mathematica code)
    # Extinction: consumer values < 0.005
    extinction <- any(consumer_vals < 0.005)
    
    if (extinction) {
      return("extinction")
    }
    
    # For non-extinction cases, check for oscillations vs equilibrium
    # Look at variation in the final part
    consumer_var <- var(consumer_vals)
    consumer_range <- max(consumer_vals) - min(consumer_vals)
    
    # Equilibrium: low variation
    if (consumer_range < 0.1 && consumer_var < 0.01) {
      return("equilibrium")
    }
    
    # Limit cycles: intermediate oscillations
    if (consumer_range > 0.1 && max(consumer_vals) < 0.35) {
      return("limit_cycle")
    }
    
    # Large amplitude oscillations
    return("oscillatory")
    
  }, error = function(e) {
    return("failed")
  })
}

# Run simulations (this will take some time)
cat("Running simulations... This may take several minutes.\n")

# For faster computation, you might want to use parallel processing
library(parallel)
results <- mcmapply(simulate_single,
                   param_grid$variance,
                   param_grid$phenotypic_dist,
                   mc.cores = 4)

# Sequential version (safer but slower)
# results <- mapply(simulate_single, 
#                   param_grid$variance, 
#                   param_grid$phenotypic_dist)

# Create results data frame
results_df <- data.frame(
  variance = param_grid$variance,
  phenotypic_dist = param_grid$phenotypic_dist,
  outcome = results
)

# Create the bifurcation plot
p <- ggplot(results_df, aes(x = phenotypic_dist, y = variance, color = outcome)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c(
    "extinction" = "black",
    "equilibrium" = "white", 
    "limit_cycle" = "gray60",
    "oscillatory" = "gray80",
    "failed" = "red"
  )) +
  labs(
    x = "Phenotypic Distance",
    y = expression(sigma^2),
    color = "Dynamics",
    title = "Heat Map of Dynamical Outcomes (Gibert & Brassil 2014, Fig 3)"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    legend.position = "right",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5)
  ) +
  scale_x_continuous(breaks = c(0, 1, 2)) +
  scale_y_continuous(breaks = c(0, 1, 2)) +
  coord_fixed(ratio = 1)

# Display the plot
print(p)

# Optional: Save the plot
# ggsave("gibert_brassil_fig3_reproduction.png", p, width = 10, height = 8, dpi = 300)

# Print summary of results
cat("\nSimulation Summary:\n")
table(results_df$outcome)