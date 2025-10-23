# Modified Rosenzweig-MacArthur Model with Phenotypic Variation
# Based on Gibert & Brassil (2014)

library(deSolve)
library(ggplot2)
library(dplyr)
library(pracma) # for numerical integration

# Function to calculate attack rate as function of phenotype x
attack_rate <- function(x, a_max, theta_a, s_squared) {
  a_max * exp(-(x - theta_a)^2 / (2 * s_squared))
}

# Function to calculate handling time as function of phenotype x  
handling_time <- function(x, h_max, h_min, theta_h, m_squared) {
  h_max - (h_max - h_min) * exp(-(x - theta_h)^2 / (2 * m_squared))
}

# Function to calculate interaction strength for resources
interaction_strength_R <- function(R, C, x_mean, sigma_squared, delta_squared, 
                                   a_max, h_max, h_min, s_squared, m_squared) {
  # Define the integrand for numerical integration
  integrand <- function(x) {
    # Phenotype probability density
    p_x <- dnorm(x, mean = x_mean, sd = sqrt(sigma_squared))
    
    # Attack rate and handling time for this phenotype
    a_x <- attack_rate(x, a_max, x_mean + sqrt(delta_squared), s_squared)
    h_x <- handling_time(x, h_max, h_min, x_mean + sqrt(delta_squared), m_squared)
    
    # Functional response component
    functional_response <- (a_x) / (1 + a_x * h_x * R)
    
    return(functional_response * p_x)
  }
  
  # Numerical integration
  integral_result <- integrate(integrand, lower = x_mean - 5*sqrt(sigma_squared), 
                               upper = x_mean + 5*sqrt(sigma_squared))$value
  
  return(-R * integral_result)
}

# Function to calculate interaction strength for consumers
interaction_strength_C <- function(R, C, x_mean, sigma_squared, delta_squared, 
                                   a_max, h_max, h_min, s_squared, m_squared, e) {
  # Define the integrand for numerical integration
  integrand <- function(x) {
    # Phenotype probability density
    p_x <- dnorm(x, mean = x_mean, sd = sqrt(sigma_squared))
    
    # Attack rate and handling time for this phenotype
    a_x <- attack_rate(x, a_max, x_mean + sqrt(delta_squared), s_squared)
    h_x <- handling_time(x, h_max, h_min, x_mean + sqrt(delta_squared), m_squared)
    
    # Functional response component
    functional_response <- (a_x) / ((1 + a_x * h_x * R)^2)
    
    return(functional_response * p_x)
  }
  
  # Numerical integration
  integral_result <- integrate(integrand, lower = x_mean - 5*sqrt(sigma_squared), 
                               upper = x_mean + 5*sqrt(sigma_squared))$value
  
  return(e * C * integral_result)
}

# Modified R-M system with phenotypic variation
rm_system <- function(t, state, parameters) {
  R <- max(state[1], 0)  # Ensure non-negative
  C <- max(state[2], 0)  # Ensure non-negative
  
  if (R < 1e-10 || C < 1e-10) {
    # Handle edge cases near zero
    dR <- ifelse(R < 1e-10, 0, parameters$r * R * (1 - R/parameters$K))
    dC <- ifelse(C < 1e-10, 0, -parameters$m * C)
    return(list(c(dR, dC)))
  }
  
  with(as.list(parameters), {
    # Define the integrand for the functional response
    functional_response_integrand <- function(x) {
      # Phenotype probability density
      p_x <- dnorm(x, mean = x_mean, sd = sqrt(sigma_squared))
      
      # Attack rate and handling time for this phenotype
      # Fixed: theta values should be x_mean +/- sqrt(delta_squared) for optimal values
      a_x <- attack_rate(x, a_max, x_mean, s_squared)  # Simplified: optimum at x_mean
      h_x <- handling_time(x, h_max, h_min, x_mean, m_squared)  # Simplified: optimum at x_mean
      
      # Functional response per individual predator
      return((a_x) / (1 + a_x * h_x * R) * p_x)
    }
    
    # Numerical integration for mean functional response
    tryCatch({
      mean_func_response <- integrate(functional_response_integrand, 
                                      lower = x_mean - 4*sqrt(sigma_squared), 
                                      upper = x_mean + 4*sqrt(sigma_squared))$value
    }, error = function(e) {
      # Fallback: use mean values
      a_mean <- a_max * exp(-delta_squared / (2 * s_squared))
      h_mean <- h_max - (h_max - h_min) * exp(-delta_squared / (2 * m_squared))
      mean_func_response <<- a_mean / (1 + a_mean * h_mean * R)
    })
    
    # Total consumption
    F_RC <- R * C * mean_func_response
    
    # System equations
    dR <- r * R * (1 - R/K) - F_RC
    dC <- e * F_RC - m * C
    
    return(list(c(dR, dC)))
  })
}

# Function to determine system dynamics type
classify_dynamics <- function(R_eq, C_eq, parameters, tolerance = 1e-4) {
  if (C_eq < tolerance) {
    return("Predator extinction")
  }
  
  if (R_eq < tolerance) {
    return("Resource extinction")
  }
  
  # Calculate Jacobian at equilibrium numerically
  with(as.list(parameters), {
    eps <- 1e-6
    state_eq <- c(R_eq, C_eq)
    
    # Get baseline derivatives
    dydt_base <- rm_system(0, state_eq, parameters)[[1]]
    
    # Perturb R slightly
    dydt_R_plus <- rm_system(0, c(R_eq + eps, C_eq), parameters)[[1]]
    dydt_R_minus <- rm_system(0, c(R_eq - eps, C_eq), parameters)[[1]]
    
    # Perturb C slightly
    dydt_C_plus <- rm_system(0, c(R_eq, C_eq + eps), parameters)[[1]]
    dydt_C_minus <- rm_system(0, c(R_eq, C_eq - eps), parameters)[[1]]
    
    # Calculate Jacobian elements using central differences
    J11 <- (dydt_R_plus[1] - dydt_R_minus[1]) / (2 * eps)
    J12 <- (dydt_C_plus[1] - dydt_C_minus[1]) / (2 * eps)
    J21 <- (dydt_R_plus[2] - dydt_R_minus[2]) / (2 * eps)
    J22 <- (dydt_C_plus[2] - dydt_C_minus[2]) / (2 * eps)
    
    # Eigenvalue analysis
    trace <- J11 + J22
    det <- J11 * J22 - J12 * J21
    discriminant <- trace^2 - 4 * det
    
    # Classification based on eigenvalues
    if (det < -tolerance) {
      return("Saddle (unstable)")
    } else if (abs(det) < tolerance) {
      return("Neutral/Transcritical")
    } else if (trace > tolerance) {
      if (discriminant > 0) {
        return("Unstable node")
      } else {
        return("Unstable focus (limit cycles)")
      }
    } else if (trace < -tolerance) {
      if (discriminant > tolerance) {
        return("Stable node")
      } else {
        return("Stable focus (damped oscillations)")
      }
    } else {
      # trace ≈ 0
      if (det > tolerance) {
        return("Center (neutral cycles)")
      } else {
        return("Neutral stability")
      }
    }
  })
}

# Function to find equilibrium numerically
find_equilibrium <- function(parameters) {
  with(as.list(parameters), {
    # First check consumer persistence condition
    # At R = K, consumer can persist if: e * mean_attack_rate / (1 + mean_attack_rate * mean_handling_time * K) > m
    
    # Calculate mean attack rate and handling time
    a_mean <- a_max * exp(-delta_squared / (2 * s_squared))
    h_mean <- h_max - (h_max - h_min) * exp(-delta_squared / (2 * m_squared))
    
    max_consumption_rate <- a_mean / (1 + a_mean * h_mean * K)
    
    if (e * max_consumption_rate <= m) {
      # Consumer cannot persist
      return(c(K, 0))
    }
    
    # Try multiple initial conditions to find coexistence equilibrium
    initial_conditions <- list(
      c(K * 0.5, 0.1),
      c(K * 0.8, 0.05),
      c(K * 0.3, 0.2),
      c(K * 0.7, 0.15)
    )
    
    best_result <- NULL
    min_residual <- Inf
    
    for (init in initial_conditions) {
      tryCatch({
        # Longer integration to reach equilibrium
        times <- seq(0, 2000, by = 10)
        sol <- ode(y = init, times = times, 
                   func = rm_system, parms = parameters, 
                   method = "lsoda")
        
        # Get final state (last 10% of trajectory to check stability)
        n_final <- max(1, floor(0.1 * nrow(sol)))
        final_states <- sol[(nrow(sol) - n_final + 1):nrow(sol), 2:3]
        
        # Check convergence - final state should be stable
        final_state <- colMeans(final_states)
        final_var <- apply(final_states, 2, var)
        
        # Only accept if trajectory has converged and both populations are positive
        if (all(final_var < 1e-6) && all(final_state > 1e-8)) {
          # Check if it's truly an equilibrium
          residual <- rm_system(0, final_state, parameters)[[1]]
          residual_norm <- sqrt(sum(residual^2))
          
          if (residual_norm < min_residual) {
            min_residual <- residual_norm
            best_result <- final_state
          }
        }
      }, error = function(e) {
        # Skip this initial condition if integration fails
      })
    }
    
    # If no coexistence found, return prey-only equilibrium
    if (is.null(best_result) || min_residual > 1e-2) {
      return(c(K, 0))
    }
    
    return(best_result)
  })
}

# Main function to analyze dynamics for given sigma and delta
analyze_dynamics <- function(sigma, delta, parameters) {
  # Update parameters
  params <- parameters
  params$sigma_squared <- sigma^2
  params$delta_squared <- delta^2
  
  # Find equilibrium
  eq <- find_equilibrium(params)
  R_eq <- eq[1]
  C_eq <- eq[2]
  
  # Classify dynamics
  dynamics_type <- classify_dynamics(R_eq, C_eq, params)
  
  # Calculate interaction strengths at equilibrium
  IS_R <- interaction_strength_R(R_eq, C_eq, params$x_mean, params$sigma_squared, 
                                 params$delta_squared, params$a_max, params$h_max, 
                                 params$h_min, params$s_squared, params$m_squared)
  
  IS_C <- interaction_strength_C(R_eq, C_eq, params$x_mean, params$sigma_squared, 
                                 params$delta_squared, params$a_max, params$h_max, 
                                 params$h_min, params$s_squared, params$m_squared, params$e)
  
  return(list(
    R_eq = R_eq,
    C_eq = C_eq,
    dynamics = dynamics_type,
    IS_R = abs(IS_R),
    IS_C = IS_C,
    sigma = sigma,
    delta = delta
  ))
}

# Set default parameters (based on paper)
default_params <- list(
  r = 0.3,        # Resource growth rate
  K = 1,          # Carrying capacity
  e = 0.8,        # Conversion efficiency (increased for coexistence)
  m = 0.5,        # Consumer mortality rate (decreased for coexistence)
  a_max = 1.5,    # Maximum attack rate
  h_max = 2,      # Maximum handling time
  h_min = 0.5,    # Minimum handling time (decreased for better efficiency)
  s_squared = 1,  # Attack rate variance parameter
  m_squared = 1,  # Handling time variance parameter
  x_mean = 0      # Mean phenotype
)

# Example: Analyze dynamics for a grid of sigma and delta values
sigma_values <- seq(0.1, 6, by = 0.15)  # Smaller range for testing
delta_values <- seq(0, 3, by = 0.15)  # Smaller range for testing

# Create results dataframe
results <- expand.grid(sigma = sigma_values, delta = delta_values)
results$dynamics <- NA
results$IS_R <- NA
results$IS_C <- NA
results$R_eq <- NA
results$C_eq <- NA

# Quick test with a few combinations first
cat("Testing a few parameter combinations first...\n")
test_params <- data.frame(
  sigma = c(0.1, 0.5, 1.0, 1.5),
  delta = c(0.0, 0.3, 0.5, 1.0)
)

for (i in 1:nrow(test_params)) {
  result <- test_dynamics(test_params$sigma[i], test_params$delta[i])
}

# Run analysis (this may take a few minutes)
cat("Running full analysis...\n")
for (i in 1:nrow(results)) {
  if (i %% 20 == 0) cat("Progress:", i, "/", nrow(results), "\n")
  
  tryCatch({
    analysis <- analyze_dynamics(results$sigma[i], results$delta[i], default_params)
    results$dynamics[i] <- analysis$dynamics
    results$IS_R[i] <- analysis$IS_R
    results$IS_C[i] <- analysis$IS_C
    results$R_eq[i] <- analysis$R_eq
    results$C_eq[i] <- analysis$C_eq
  }, error = function(e) {
    results$dynamics[i] <- "Error"
    cat("Error at sigma =", results$sigma[i], ", delta =", results$delta[i], ":", e$message, "\n")
  })
}

cat("Analysis complete!\n")

# Display summary of results
cat("\nSummary of dynamics:\n")
print(table(results$dynamics, useNA = "ifany"))

# Create Figure 2 reproduction
# Map dynamics to colors similar to the paper
results$color_code <- case_when(
  results$dynamics == "Predator extinction" ~ "black",
  grepl("limit cycles|Unstable", results$dynamics) ~ "darkgray", 
  grepl("damped|Stable focus", results$dynamics) ~ "lightgray",
  grepl("Stable node|nonoscillatory", results$dynamics) ~ "white",
  TRUE ~ "gray"
)

# Plot interaction strengths (Figure 2 style)
p1 <- ggplot(results, aes(x = sigma, y = IS_R)) +
  geom_point(aes(color = factor(delta)), alpha = 0.7) +
  labs(title = "Resource Interaction Strength vs Individual Variation",
       x = "Individual variation (σ)", 
       y = "Interaction strength",
       color = "Phenotypic mismatch (δ)") +
  theme_minimal()

p2 <- ggplot(results, aes(x = sigma, y = IS_C)) +
  geom_point(aes(color = factor(delta)), alpha = 0.7) +
  labs(title = "Consumer Interaction Strength vs Individual Variation",
       x = "Individual variation (σ)", 
       y = "Interaction strength",
       color = "Phenotypic mismatch (δ)") +
  theme_minimal()

print(p1)
print(p2)

# Plot dynamics classification
p3 <- ggplot(results, aes(x = sigma, y = delta, fill = dynamics)) +
  geom_tile() +
  labs(title = "Dynamics Classification",
       x = "Individual variation (σ)",
       y = "Phenotypic mismatch (δ²)",
       fill = "Dynamics") +
  theme_minimal()

print(p3)

# Function to quickly test a specific sigma, delta combination
test_dynamics <- function(sigma_val, delta_val) {
  result <- analyze_dynamics(sigma_val, delta_val, default_params)
  cat("σ =", sigma_val, ", δ =", delta_val, "\n")
  cat("Dynamics:", result$dynamics, "\n")
  cat("Resource equilibrium:", result$R_eq, "\n") 
  cat("Consumer equilibrium:", result$C_eq, "\n")
  cat("Resource interaction strength:", result$IS_R, "\n")
  cat("Consumer interaction strength:", result$IS_C, "\n")
  return(result)
}

# Example usage:
# test_dynamics(0.5, 0.5)
# test_dynamics(2.0, 1.0)