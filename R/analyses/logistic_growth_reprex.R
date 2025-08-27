n_sim <- 1e4
r_sim <- rnorm(n_sim, 10, .1)
hist(r_sim)

logistic_model_vi <- odin({
  # States:
  # R: population density of resources (prey)
  # Parameters:
  # r = per capita growth rate (Resource)
  # K = Resource carrying capacity
  
  deriv(R[]) <- r[i]*R[i]*(1 - R[i]/K)
  
  # Initial conditions
  initial(R[]) <- R0[i]
  
  # Define inital values
  R0[] <- user() # Initial population density
  
  # Reasonable default values (following Gibert & Brassil 2014)
  r[] <- user()
  K <- user(1)
  
  # Dimensions
  dim(r) <- user()
  n_val <- length(r)
  dim(R0) <- n_val
  dim(R) <- n_val
  
})

pars <- list(
  # initial values
  R0 = rep(100, n_sim),
  # range of a and h values
  r = r_sim)

mod <- logistic_model_vi$new(user = pars)

t <- seq(0, 2000, length.out = 10001)
y <- mod$run(t)
str(y)
plot(y)
