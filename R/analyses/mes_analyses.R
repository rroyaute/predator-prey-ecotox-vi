library(odin); library(tidyverse)

# 1. Rosenzweig-MacArthur Consummer-Resource model without variation ----
# Calculate parameters a and h according to Gibert & Brassil 2014
a_max <- 2
h_min <- 1
h_max <- 2
tau <- 1
nu <- 1
d_a <- .5
d_h <- .5

a <- a_max * exp(-d_a^2/(2 * tau^2))
h <- h_max - (h_max - h_min) * exp(-d_h^2/(2 * nu^2))

RM_model <- odin({
  # States:
  # R: population density of resources (prey)
  # C: population density of consumers (predators)
  # Parameters:
  # r = per capita growth rate (Resource)
  # K = Resource carrying capacity
  # a = attack rate
  # h = handling time
  # e = conversion efficiency
  # d = Consummer death rate
  
  deriv(R) <- r*R*(1 - R/K) - (a*R*C)/(1 + a*h*R)
  deriv(C) <- e*(a*R*C)/(1 + a*h*R) - d*C
  
  # Initial conditions
  initial(R) <- R0
  initial(C) <- C0
  
  # Define values
  R0 <- user(1) # Initial resource population density
  C0 <- user(0.05) # Initial consummer population density
  
  # Reasonable default values (following Gibert & Brassil 2014)
  r <- user(.3)
  K <- user(1)
  a <- user(1.764994)
  h <- user(1.117503)
  e <- user(.5)
  d <- user(.1)
  
  })

mod <- RM_model$new()

t <- seq(0, 1000, by = .5)
y <- mod$run(t)
plot(y)

y %>% 
  as.data.frame() %>% 
  pivot_longer(cols = c("R", "C"),
               names_to = "Density_type",
               values_to = "Density") %>% 
  ggplot(aes(x = t, y = Density, 
             group = Density_type, 
             color = Density_type)) +
  geom_line() +
  theme_bw()



# 2. Adding individual variation -----
# Defining phenotype x which follows a normal distribution and influence
# a and h
a_max <- 2
h_min <- 1
h_max <- 2
tau <- 1
nu <- 1
d_a <- .5
d_h <- .5
theta_a <- -sqrt(.5)
theta_h <- -sqrt(.5)

# Suppose x follows centered normal distribution
n_sim <- 1e4
x <- rnorm(n = n_sim, mean = 0, sd = .3)
hist(x)

# Store a and h vectors
a_sim <- a_max * exp(-(x - theta_a)^2/(2 * tau^2))
h_sim <- h_max - (h_max - h_min) * exp(-(x - theta_h)^2/(2 * nu^2))
# plot(a_sim,h_sim)

RM_model_vi <- odin({
  # States:
  # R: population density of resources (prey)
  # C: population density of consumers (predators)
  # Parameters:
  # r = per capita growth rate (Resource)
  # K = Resource carrying capacity
  # a = attack rate (vector)
  # h = handling time  (vector)
  # e = conversion efficiency
  # d = Consummer death rate
  
  deriv(R[]) <- r*R[i]*(1 - R[i]/K) - (a[i]*R[i]*C[i])/(1 + a[i]*h[i]*R[i])
  deriv(C[]) <- e*(a[i]*R[i]*C[i])/(1 + a[i]*h[i]*R[i]) - d*C[i]
  
  # Initial conditions
  initial(R[]) <- R0[i]
  initial(C[]) <- C0[i]
  
  # Define inital values
  R0[] <- user() # Initial resource population density
  C0[] <- user() # Initial consummer population density
  
  # Reasonable default values (following Gibert & Brassil 2014)
  r <- user(.3)
  K <- user(1)
  a[] <- user() # vector of a values
  h[] <- user() # vector of h values
  e <- user(.5)
  d <- user(.1)
  
  # Dimensions
  dim(a) <- user()
  n_val <- length(a)
  dim(R0) <- n_val
  dim(C0) <- n_val
  dim(R) <- n_val
  dim(C) <- n_val
  dim(h) <- n_val
  
})

pars <- list(
  # initial values
  R0 = rep(1, n_sim),
  C0 = rep(0.05, n_sim),
  # range of a and h values
  a = a_sim, h = h_sim)

RM_model_vi_mod <- RM_model_vi$new(user = pars)

t <- seq(0, 1000, by = .5)
y <- RM_model_vi_mod$run(t)
plot(y)

# 3. Modeling a dose-response on x ----
# Assume x at optimum theta_a & theta_h at dose = 0
n_sim <- 1e4
log_x <- rnorm(n = n_sim, mean = 0, sd = .3)
hist(log_x); hist(exp(log_x)) 



