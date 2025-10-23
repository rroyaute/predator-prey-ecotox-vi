library(deSolve)
library(ggplot2)

# Parameters from Gibert & Brassil 2014
params_default <- list(
  r = 0.3,        # Resource growth rate
  K = 1,          # Carrying capacity  
  e = 0.5,        # Conversion efficiency
  m = 0.1,        # Consumer mortality (they use 'd' in paper)
  amax = 2,       # Maximum attack rate
  hmin = 1,       # Minimum handling time (gmin in paper)
  hmax = 2,       # Maximum handling time (gmax in paper)
  tau = 1,        # Width for attack rate (s in paper)
  nu = 1          # Width for handling time (μ in paper)
)

# Calculate mean attack rate and handling time with Jensen's inequality
calculate_mean_params <- function(sigma, d, params) {
  if (sigma <= 0) {
    # No variation case
    a_mean <- params$amax * exp(-d^2 / (2 * params$tau^2))
    h_mean <- params$hmax - (params$hmax - params$hmin) * exp(-d^2 / (2 * params$nu^2))
    return(list(a = a_mean, h = h_mean))
  }
  
  # Integration for Jensen's inequality effect
  n_points <- 1000
  x_range <- max(5 * sigma, abs(d) + 5 * sigma)
  x <- seq(-x_range, x_range, length.out = n_points)
  dx <- diff(x)[1]
  
  # Trait distribution
  p_x <- dnorm(x, mean = 0, sd = sigma)
  
  # Attack rate and handling time as functions of trait
  a_x <- params$amax * exp(-(x - d)^2 / (2 * params$tau^2))
  h_x <- params$hmax - (params$hmax - params$hmin) * exp(-(x - d)^2 / (2 * params$nu^2))
  
  # Mean values
  a_mean <- sum(a_x * p_x) * dx
  h_mean <- sum(h_x * p_x) * dx
  
  return(list(a = a_mean, h = h_mean))
}

# Find equilibrium points
find_equilibrium <- function(sigma, d, params = params_default) {
  # Get mean parameters
  mp <- calculate_mean_params(sigma, d, params)
  a <- mp$a
  h <- mp$h
  
  # Consumer isocline: R* = m/(e*a - m*h)
  # This is the resource level where dC/dt = 0
  denominator <- params$e * a - params$m * h
  
  if (denominator <= 0) {
    # Consumer cannot persist
    return(list(exists = FALSE, R = params$K, C = 0))
  }
  
  R_star <- params$m / denominator
  
  # Check if R* is feasible
  if (R_star <= 0 || R_star >= params$K) {
    return(list(exists = FALSE, R = params$K, C = 0))
  }
  
  # Resource isocline at R*: C* = (r/a)(1 - R*/K)(1 + a*h*R*)
  C_star <- (params$r / a) * (1 - R_star/params$K) * (1 + a * h * R_star)
  
  if (C_star <= 0) {
    return(list(exists = FALSE, R = params$K, C = 0))
  }
  
  return(list(exists = TRUE, R = R_star, C = C_star, a = a, h = h))
}

# Stability analysis at equilibrium
analyze_stability <- function(sigma, d, params = params_default) {
  eq <- find_equilibrium(sigma, d, params)
  
  if (!eq$exists) {
    return("extinction")
  }
  
  R_star <- eq$R
  C_star <- eq$C
  a <- eq$a
  h <- eq$h
  
  # Jacobian matrix elements
  # For Rosenzweig-MacArthur model at equilibrium
  
  # Common term
  denom <- (1 + a * h * R_star)^2
  
  # Jacobian elements
  J11 <- params$r * (1 - 2*R_star/params$K) - (a * C_star) / denom
  J12 <- -a * R_star / (1 + a * h * R_star)
  J21 <- params$e * a * C_star / denom
  J22 <- 0  # Always 0 at equilibrium for this model
  
  # Eigenvalue analysis
  trace <- J11  # Since J22 = 0
  det <- -J12 * J21  # Since J22 = 0
  
  # For stability, we need det > 0
  if (det <= 0) {
    return("unstable")
  }
  
  # discriminant of characteristic equation
  disc <- trace^2 - 4*det
  
  if (disc > 0) {
    # Real eigenvalues
    lambda1 <- (trace + sqrt(disc))/2
    lambda2 <- (trace - sqrt(disc))/2
    
    if (lambda1 < 0 && lambda2 < 0) {
      return("stable")
    } else {
      return("unstable")
    }
  } else {
    # Complex eigenvalues
    real_part <- trace/2
    
    if (real_part < -0.001) {
      # Stable spiral - need simulation to distinguish damped vs limit cycle
      return(check_oscillation_type(R_star, C_star, a, h, params))
    } else if (real_part > 0.001) {
      return("unstable")
    } else {
      # Center - structurally unstable, treat as limit cycle
      return("limit_cycle")
    }
  }
}

# Simulate to distinguish oscillation types
check_oscillation_type <- function(R_star, C_star, a, h, params) {
  ode_system <- function(t, state, parms) {
    R <- state[1]
    C <- state[2]
    
    dR <- parms$r * R * (1 - R/parms$K) - parms$a * R * C / (1 + parms$a * parms$h * R)
    dC <- parms$e * parms$a * R * C / (1 + parms$a * parms$h * R) - parms$m * C
    
    list(c(dR, dC))
  }
  
  # Parameters for simulation
  sim_params <- params
  sim_params$a <- a
  sim_params$h <- h
  
  # Initial condition (small perturbation)
  y0 <- c(R = R_star * 1.01, C = C_star * 1.01)
  times <- seq(0, 5000, by = 0.5)
  
  sol <- tryCatch({
    ode(y0, times, ode_system, sim_params, method = "lsoda", atol = 1e-10, rtol = 1e-10)
  }, error = function(e) NULL)
  
  if (is.null(sol) || nrow(sol) < 1000) {
    return("unknown")
  }
  
  # Analyze last quarter of simulation
  n <- nrow(sol)
  last_quarter <- sol[(3*n/4):n, ]
  R_vals <- last_quarter[, 2]
  
  # Find peaks
  dR <- diff(R_vals)
  peaks <- which(diff(sign(dR)) == -2) + 1
  
  if (length(peaks) < 3) {
    # No clear oscillations
    cv <- sd(R_vals) / mean(R_vals)
    if (cv < 0.001) {
      return("stable")
    } else {
      return("damped_oscillations")
    }
  }
  
  # Check peak heights
  peak_heights <- R_vals[peaks]
  
  # Are peaks decreasing (damped) or constant (limit cycle)?
  if (length(peak_heights) >= 5) {
    # Use last 5 peaks
    last_peaks <- tail(peak_heights, 5)
    peak_trend <- cor(1:5, last_peaks)
    
    if (abs(peak_trend) < 0.1 && sd(last_peaks)/mean(last_peaks) < 0.01) {
      return("limit_cycle")
    } else if (peak_trend < -0.3) {
      return("damped_oscillations")
    } else {
      # Check coefficient of variation
      cv <- sd(last_peaks) / mean(last_peaks)
      if (cv < 0.02) {
        return("limit_cycle")
      } else {
        return("damped_oscillations")
      }
    }
  } else {
    # Not enough peaks, use amplitude
    amplitude <- max(R_vals) - min(R_vals)
    relative_amp <- amplitude / mean(R_vals)
    
    if (relative_amp > 0.1) {
      return("limit_cycle")
    } else {
      return("damped_oscillations")
    }
  }
}

# Create regime map
create_regime_map <- function(sigma_range = seq(0, 3, by = 0.1),
                              d_range = seq(0, 2, by = 0.1),
                              params = params_default) {
  
  grid <- expand.grid(sigma = sigma_range, d = d_range)
  grid$regime <- NA_character_
  
  cat("Computing regime map...\n")
  pb <- txtProgressBar(max = nrow(grid), style = 3)
  
  for (i in 1:nrow(grid)) {
    grid$regime[i] <- analyze_stability(grid$sigma[i], grid$d[i], params)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  return(grid)
}

# Plotting function
plot_regime_map <- function(grid) {
  # Make sure regime is a factor with all possible levels
  grid$regime <- factor(grid$regime,
                        levels = c("extinction", "stable", 
                                   "damped_oscillations", "limit_cycle",
                                   "unstable", "unknown"))
  
  # Colors similar to the paper
  colors <- c(
    "extinction" = "black",
    "stable" = "white",
    "damped_oscillations" = "#CCCCCC",  # Light gray
    "limit_cycle" = "#666666",           # Dark gray
    "unstable" = "red",
    "unknown" = "yellow"
  )
  
  # Create the plot with proper axes
  p <- ggplot(grid, aes(x = sigma, y = d, fill = regime)) +
    geom_tile() +
    scale_fill_manual(values = colors, drop = FALSE, na.value = "gray50") +
    scale_x_continuous(expand = c(0, 0), 
                       breaks = seq(0, 2.5, 0.5),
                       labels = function(x) x^2) +  # Show as sigma^2
    scale_y_continuous(expand = c(0, 0),
                       breaks = seq(0, 1.5, 0.3),
                       labels = function(x) x^2) +  # Show as d^2
    labs(x = expression("Individual variation ("*sigma^2*")"),
         y = expression("Phenotypic mismatch ("*d^2*")"),
         fill = "Regime",
         title = "Predator-Prey Dynamics (cf. Gibert & Brassil 2014, Fig. 3)") +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    coord_cartesian(xlim = c(0, 2.5), ylim = c(0, 1.5))
  
  return(p)
}

# Test function for debugging
test_point <- function(sigma, d, params = params_default) {
  cat(sprintf("\nTesting sigma = %.2f, d = %.2f\n", sigma, d))
  cat(paste(rep("=", 40), collapse=""), "\n")
  
  # Mean parameters
  mp <- calculate_mean_params(sigma, d, params)
  cat(sprintf("Mean attack rate: %.4f\n", mp$a))
  cat(sprintf("Mean handling time: %.4f\n", mp$h))
  
  # Equilibrium
  eq <- find_equilibrium(sigma, d, params)
  if (eq$exists) {
    cat(sprintf("Equilibrium: R* = %.4f, C* = %.4f\n", eq$R, eq$C))
  } else {
    cat("No coexistence equilibrium\n")
  }
  
  # Stability
  regime <- analyze_stability(sigma, d, params)
  cat(sprintf("Regime: %s\n", regime))
  
  return(invisible(list(params = mp, equilibrium = eq, regime = regime)))
}

# Test with debug output
debug_counter <<- 0
cat("\nRe-testing points with corrected Jacobian:\n")
cat(paste(rep("=", 50), collapse=""), "\n")

test_results <- list()
test_params <- list(
  list(0.5, 0.0, "Low variation, no mismatch"),
  list(1.0, 0.5, "Moderate both"),
  list(0.5, 1.5, "Low variation, high mismatch"),
  list(2.0, 0.2, "High variation, low mismatch")
)

for (tp in test_params) {
  cat(sprintf("\n%s (σ=%.1f, d=%.1f):\n", tp[[3]], tp[[1]], tp[[2]]))
  result <- test_point(tp[[1]], tp[[2]])
  test_results[[length(test_results) + 1]] <- result
}

# Remove debug counter
rm(debug_counter)

# Generate and plot the regime map
cat("\n\nGenerating full regime map...\n")
regime_map <- create_regime_map(
  sigma_range = seq(0, sqrt(6), length.out = 25),
  d_range = seq(0, sqrt(3), length.out = 20),
  params = params_default
)

# Plot
p <- plot_regime_map(regime_map)
print(p)

# Summary
cat("\nRegime distribution:\n")
print(table(regime_map$regime))

# Test what happens at high d values
cat("\n\nDiagnosing high d behavior:\n")
cat(paste(rep("=", 40), collapse=""), "\n")

# Test with increasing d values
for (d_test in c(0.5, 1.0, 1.5, 2.0)) {
  mp <- calculate_mean_params(sigma = 0.5, d = d_test, params_default)
  invasion_check <- params_default$e * mp$a - params_default$m * mp$h
  cat(sprintf("d = %.1f: mean_a = %.4f, mean_h = %.4f, invasion = %.4f\n", 
              d_test, mp$a, mp$h, invasion_check))
  
  if (invasion_check > 0) {
    R_star <- params_default$m / invasion_check
    cat(sprintf("  -> R* = %.4f (K = %.1f)\n", R_star, params_default$K))
  } else {
    cat("  -> Cannot invade (extinction)\n")
  }
}

# Generate and plot with alternative parameters
cat("\n\nTrying with alternative parameters...\n")
regime_map_alt <- create_regime_map(
  sigma_range = seq(0, 2.5, by = 0.15),
  d_range = seq(0, 2.0, by = 0.1),
  params = params_alternative
)

p_alt <- plot_regime_map(regime_map_alt)
print(p_alt)

cat("\nRegime distribution (alternative params):\n")
print(table(regime_map_alt$regime))