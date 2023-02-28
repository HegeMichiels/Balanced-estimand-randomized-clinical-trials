###This file contains all functions required to estimate the balanced estimand mu = mu_1-mu_0, 
### "The treatment effect that would have been observed if all patients switch as under control treatment"

##Library ####
library(clusterPower)
library(nleqslv)
library(ggplot2)
library(reshape2)
library(plyr)

##Functions ####

##Function to simulate data
#Input: n = total sample size
#Output: dataset with variables R,C,L,S and Y, as defined in the paper
simulate.data = function(n){
  R = rbinom(n, 1, 0.5) #Treatment arm
  C = rnorm(n,0,1) #Baseline covariate 
  L.one = rnorm(n, mean = delta1 + delta2*C, sd = sd.l) #Severity of disease
  #Treatment group
  S.one = rbinom(n,size = 1,prob = expit(omega1+omega2*C+omega3*L.one)) #Switching under treatment
  Y.one = rnorm(n, mean = alpha1 + alpha2*S.one+alpha3*L.one+alpha4*C , sd = sd.y) #Outcome under treatment
  #Control group
  S.zero = rbinom(n,size = 1, prob = expit(lambda1 +lambda2*C + rho*omega3*L.one)) #Switching under treatment
  Y.zero = rnorm(n,mean = alpha5+alpha1+alpha2*S.zero+alpha3*L.one+alpha4*C,sd=sd.y) #Outcome under control
  #Observed data
  L = ifelse(R ==1, yes = L.one,no= NA )
  S = R*S.one + (1-R)*S.zero
  Y = R*Y.one + (1-R)*Y.zero
  data = data.frame(R,C, L,S,Y)
  return(data)
}

##Function to estimate the parameters mu, mu1 and mu0 using the proposed IPW estimator
#Inputs: data = dataset with variables R,C,L,S and Y, as defined in the paper
#        rho = sensitivity parameter
#Output as vector: mu = balanced effect estimate (=mu1-mu0)
#                  mu1 = expected outcome in treatment arm, with switching as under control
#                  mu0 = expected outcome in the control arm
IPW = function(data,rho){
  #pi
  prob.R = predict(glm(R~C,data=data,family=binomial),newdata = data.frame(C=data$C),type="response") #probability to be treated, conditional on C
  #l(L^1,C)
  fit.l = glm(S~L+C, data=subset(data,R==1), family = binomial(link="logit")) #Logistic regression model for P(S=1|L,R=1,C)
  prob.l = ifelse(data$R ==1,yes =  predict(fit.l,newdata = data.frame(C=data$C,L=data$L), type= "response"), no = 0) #Estimated probabilities P(S=1|L,R=1,C)
  omega3.est = coef(fit.l)["L"] #coefficient of L in logistic regression model for P(S=1|L,R=1,C)
  omegas = c(coef(fit.l)["(Intercept)"],coef(fit.l)["C"] )
  #q1 
  q1 = ifelse(data$R == 1,yes = IPW.calculate.q1(data,rho,omega3.est), no = 0) #q1 = (rho-1)*omega3*L
  #q0
  lambda = nleqslv(fn = IPW.calculate.q0, x = c(omegas[1],omegas[2]),data=data, prob.l = prob.l,q1=q1,prob.R=prob.R,omegas=omegas)$x #Estimation of the parameters lambda
  q0 = lambda1-omegas[1] + (lambda2-omegas[2])*data$C 
  #W 
  weights = ifelse(data$R == 1, yes = IPW.weights.W(data,prob.l,q0,q1), no =0) #Weights W, as defined in the main paper
  data$W = weights
  #parameters
  mu1 = IPW.mu1(data,prob.R) #Estimation of mu1
  mu0 = IPW.mu0(data,prob.R) #Estimation of mu0
  mu = mu1-mu0 #Estimation of mu
  #output
  output = c(mu,mu1,mu0)
  names(output) = c("mu","mu1","mu0")
  return(output)
}

##Function that returns q1 = (rho-1)*omega3*L
#This function is needed to estimate P(S=1|L^1,R=0,C)
#Inputs: data = dataset with variables R,C,L,S and Y, as defined in the paper
#        rho = sensitivity parameter
#        omega3 = coefficient of L in logistic regression model for P(S=1|L,R=1,C)
#Output as vector: estimates for q1, for every patient in the dataset
IPW.calculate.q1 = function(data,rho,omega3){
  L = data$L
  return((rho-1)*omega3*L)
}

##Function that return the obtained values for the estimating equations, for certain parameters lambda1 en lambda2
#Inputs: x = vector with values for lambda1 and lambda2
#        data = dataset with variables R,C,L,S and Y, as defined in the paper
#        prob.l = vector with estimated probabilities P(S=1|L,R=1,C) for every patient
#        q1 = vector with estimates for q1 for every patient 
#        prob.R = vector with the estimated probabilities to be treated, conditional on C, for every patient
#Output as vector: obtained values for the estimating equations, for the given lambda1 and lambda2 values
IPW.calculate.q0 = function(x, data, prob.l,q1,prob.R,omegas){
  R = data$R
  C = data$C
  S = data$S
  L = data$L
  lambda1 = as.numeric(x[1])
  lambda2 = as.numeric(x[2])
  omega1 = omegas[1]
  omega2 = omegas[2]
  q0 = lambda1-omega1 + (lambda2-omega2)*C
  main = (1-R)*(1-S)/(1-prob.R) - (R/prob.R)*(1-S)/(prob.l*(exp(q0+q1)-1)+1)
  return(c(sum(main),sum(C*main)))
}


##Function that returns the estimated weights W = P(S|L^1,R=0,C)/P(S|L^1,R=1,C)
#Inputs: data = dataset with variables R,C,L,S and Y, as defined in the paper
#        prob.l = vector with estimated probabilities P(S=1|L,R=1,C) for every patient
#        q0 = vector with estimates for q0 for every patient
#        q1 = vector with estimates for q1 for every patient
#Output as vector: weights W for every patient
IPW.weights.W = function(data,prob.l,q0,q1){
  R = data$R
  C = data$C
  S = data$S
  L = data$L
  Y = data$Y
  weights = exp(S*(q0+q1))/(prob.l*(exp(q0+q1)-1)+1)
  return(weights)
}

##Function that estimates mu1
#Inputs: data = dataset with variables R,C,L,S and Y, as defined in the paper and weights W
#        prob.R = vector with the estimated probabilities to be treated, conditional on C, for every patient
#Output: estimate for mu1
IPW.mu1 = function(data,prob.R){
  R = data$R
  C = data$C
  S = data$S
  L = data$L
  Y = data$Y
  W = data$W
  mu1 = mean(Y*R*W/prob.R)/mean(R*W/prob.R)
  return(mu1)
}

##Function that estimates mu0
#Inputs: data = dataset with variables R,C,L,S and Y, as defined in the paper
#        prob.R = vector with the estimated probabilities to be treated, conditional on C, for every patient
#Output: estimate for mu0
IPW.mu0 = function(data,prob.R){
  R = data$R
  Y = data$Y
  mu0 = mean((1-R)*Y/(1-prob.R))/mean((1-R)/(1-prob.R))
  return(mu0)
}
