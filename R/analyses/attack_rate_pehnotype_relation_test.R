n <- 10^6
mu <- 100
sd <- .1

plot(hist(rlnorm(n, log(mu), sd)), xlim = c(0,200))
plot(hist(rlnorm(n, log(mu/2), sd)), xlim = c(0,200))
plot(hist(rlnorm(n, log(mu/4), sd)), xlim = c(0,200))


a_max <- 100
theta_a <- 100
tau <- .25
x <- 0:200
a_x <- a_max * exp(-(log(x) - log(theta_a))^2 / (2 * tau^2))
plot(a_x ~ x, type = "l")

mu_x <- 50
sigma <- 5
x <- rnorm(n/10, mu_x, sigma*4)
a_x <- a_max * exp(-(log(x) - log(theta_a))^2 / (2 * tau^2))
mean(a_x)
plot(a_x ~ x, xlim = c(0,200))


