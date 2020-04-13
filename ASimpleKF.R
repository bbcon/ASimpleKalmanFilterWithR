# Setup
rm(list=ls())
require(tidyverse)

theme_set(theme_bw() + theme(legend.position = "bottom"))

# Parameters of the state space model
phi <- 0.8
sigma_sq_v <- 1
sigma_sq_w <- 1


# Generating data
n_periods <- 80

set.seed(1)

noise <- rnorm(mean = 0, sd = sigma_sq_w, n=n_periods)
state_noise <- rnorm(n=n_periods, mean = 0, sd = sigma_sq_v)
arima_yt <- arima.sim(list(order=c(1,0,0), ar=phi, rand.gen = state_noise), n=n_periods)
yt <- arima_yt + noise
yt <- as.vector(yt)


# Recursion algorithm

# notation: l. means the conditional expectation on t-1 while f. means the actual forecast

# Initialization
f.xi <- 0
f.P <- 100

l.xi <- 0
l.y <- 0
l.P <- 0
h_t <- 0
K_t <- 0
eta <- 0

# loop
for(i in c(2:length(yt))){
  l.xi[i] <- phi*f.xi[i-1]
  l.y[i] <- l.xi[i]
  l.P[i] <- phi*phi*f.P[i-1] + sigma_sq_v
  h_t[i] <- l.P[i] + sigma_sq_w
  K_t[i] <- l.P[i]*(h_t[i])^(-1) 
  eta[i] <- yt[i]-l.y[i]
  f.xi[i] <- l.xi[i] + K_t[i]*eta[i]
  f.P[i] <- l.P[i] - K_t[i]*l.P[i]
}

