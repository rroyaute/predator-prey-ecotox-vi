library(odin); library(tidyverse); library(ggthemes)
library(patchwork)

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

fig_pred_prey_rm <- y %>% 
  as.data.frame() %>% 
  pivot_longer(cols = c("R", "C"),
               names_to = "Density_type",
               values_to = "Density") %>% 
  ggplot(aes(x = t, y = Density, 
             group = Density_type, 
             color = Density_type)) +
  geom_line(linewidth = .7) +
  scale_color_wsj() +
  theme_bw()
fig_pred_prey_rm
ggsave("outputs/figs/fig_pred_prey_rm.jpeg", fig_pred_prey_rm)


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
theta_a <- 1
theta_h <- 1
# TODO use lognormal formulation from Appendix 5
curve(from = -10, to = 10, expr = a_max * exp(-(x - theta_a)^2) / (2 * tau^2))
curve(from = -10, to = 10, expr = h_max - (h_max - h_min) * exp(-(x - theta_h)^2) / (2 * nu^2))

n_sim <- 1e4
mu_log <- 0
sd_log <- .3
log_x <- rnorm(n = n_sim, mean = mu_log, sd = sd_log)
hist(log_x); hist(exp(log_x)) 

# dose gradient
c <- seq(0, 100, by = 1)
c[1] <- 1e-4 # replace 0 with small non-negative value
# dr parameters
max <- 0
min <- -10
ec50 <- 40
beta <- -1.5

log_x_hat <- max + (min - max) / (1 + exp((ec50 - c) * exp(beta)))
plot(c, log_x_hat, type = "l")
plot(c, exp(log_x_hat), type = "l")

df_sim <- data.frame(Concentration = c,
                     log_x_hat) %>% 
  mutate(log_x = rnorm(n(), log_x_hat, sd_log)) %>% 
  mutate(x_hat = exp(log_x_hat),
         x = exp(log_x)) %>% 
  mutate(a = a_max * exp(-(x - theta_a)^2) / (2 * tau^2),
         h = h_max - (h_max - h_min) * exp(-(x - theta_h)^2) / (2 * nu^2))

fig_drc_x <- df_sim %>% 
  ggplot(aes(x = Concentration, y = x)) +
  geom_point() +
  geom_line(aes(x = Concentration, y = x_hat), 
            color = "red", linewidth = 1.5) +
  ylab("Behavioral trait (x)") +
  ggtitle("Effect of contaminant concentration \n on behavioral expression") +
  theme_bw()
fig_drc_x  

fig_drc_a <- df_sim %>% 
  ggplot(aes(x = Concentration, y = a)) +
  geom_point() +
  ylab("Attack rate (a)") +
  ggtitle("Effect of contaminant concentration \n on attack rates") +
  theme_bw()
fig_drc_a

fig_drc_h <- df_sim %>% 
  ggplot(aes(x = Concentration, y = h)) +
  geom_point() +
  ylab("Handling time (h)") +
  ggtitle("Effect of contaminant concentration \n on handling time") +
  theme_bw()
fig_drc_h

ggsave("outputs/figs/fig_drc_x.jpeg", fig_drc_x)
ggsave("outputs/figs/fig_drc_a.jpeg", fig_drc_a)
ggsave("outputs/figs/fig_drc_h.jpeg", fig_drc_h)


# 4. Add a dose response on sigma^2 as well -----
sigma0 <- log(sd_log)
s1 <-  .06
s2 <- -.002
log_sigma <- sigma0 + s1 * c + s2 * c^2
plot(c, log_sigma, type = "l")
plot(c, exp(log_sigma), type = "l")

# Combine mean and variance dose-response
df_sim_vi <- data.frame(Concentration = c,
                     log_x_hat,
                     sd_log = exp(log_sigma)) %>% 
  mutate(log_x = rnorm(n(), log_x_hat, sd_log)) %>% 
  mutate(x_hat = exp(log_x_hat),
         x = exp(log_x)) %>% 
  mutate(a = a_max * exp(-(x - theta_a)^2) / (2 * tau^2),
         h = h_max - (h_max - h_min) * exp(-(x - theta_h)^2) / (2 * nu^2))

fig_drc_x_vi <- df_sim_vi %>% 
  ggplot(aes(x = Concentration, y = x)) +
  geom_point() +
  geom_line(aes(x = Concentration, y = x_hat), 
            color = "red", linewidth = 1.5) +
  ylab("Behavioral trait (x)") +
  ggtitle("Effect of contaminant concentration \n on behavioral expression",
          subtitle = "Case where both mean and behavioral variance is affected") +
  theme_bw()
fig_drc_x_vi  

fig_drc_a_vi <- df_sim_vi %>% 
  ggplot(aes(x = Concentration, y = a)) +
  geom_point() +
  ylab("Attack rate (a)") +
  ggtitle("Effect of contaminant concentration \n on attack rates",
          subtitle = "Case where both mean and behavioral variance is affected") +
  theme_bw()
fig_drc_a_vi

fig_drc_h_vi <- df_sim_vi %>% 
  ggplot(aes(x = Concentration, y = h)) +
  geom_point() +
  ylab("Handling time (h)") +
  ggtitle("Effect of contaminant concentration \n on handling time",
          subtitle = "Case where both mean and behavioral variance is affected") +
  theme_bw()
fig_drc_h_vi

fig_drc_sigma_vi <- df_sim_vi %>% 
  ggplot(aes(x = Concentration, y = sd_log^2)) +
  geom_line(color = "red", linewidth = 1.5) +
  ylab("Behavioral variance (sigma^2)") +
  ggtitle("Effect of contaminant concentration \n on behavioral variance") +
  theme_bw()
fig_drc_sigma_vi  

ggsave("outputs/figs/fig_drc_x_vi.jpeg", fig_drc_x_vi)
ggsave("outputs/figs/fig_drc_a_vi.jpeg", fig_drc_a_vi)
ggsave("outputs/figs/fig_drc_h_vi.jpeg", fig_drc_h_vi)
ggsave("outputs/figs/fig_drc_sigma_vi.jpeg", fig_drc_sigma_vi)

# Combine in one big plot
# fig_drc <- fig_drc_x / fig_drc_a / fig_drc_h
fig_drc_vi <- fig_drc_x_vi / fig_drc_a_vi / fig_drc_h_vi
fig_drc_vi_4pan <- fig_drc_vi | fig_drc_sigma_vi
ggsave("outputs/figs/fig_drc_vi_4pan.jpeg", fig_drc_vi_4pan, 
       height = 8, width = 8.5)
