###Use this file to generate results

#source the Balanced estimand Functions.R file 
source()

##Parameters ####

#parameters of scenario 1: 
delta1 <- -0.5 
delta2 <- 0.1
sd.l <- 0.3
omega1 <- -7
omega2 <- -0.01
omega3 <- -7
alpha1 <- 0
alpha2 <- 0.5
alpha3 <- 2
alpha4 <- 0.1
alpha5 <- -0.5
lambda1 <- -5 
lambda2 <- -0.02 
sd.y <- 0.3
rho <- 0.9

#Total sample size 
n <- 1000

##Generate data####
data <- simulate.data(n)

##Estimate the balanced estimand####
IPW(data,rho)
