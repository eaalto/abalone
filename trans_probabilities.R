################################################################
###  codice per calcolare la matrice di transizione (Bardos et al. 2005)
################################################################
# Original code by Marisa Rossetto
# Some comments, reformatting by Emil Aalto

estimate5_3  # a parameter-estimation class loaded in from "Fitting_Bardos"

g     = as.list(coef(estimate5_3))$gi
sigma = as.list(coef(estimate5_3))$sigma
Ln    = as.list(coef(estimate5_3))$Ln
beta  = as.list(coef(estimate5_3))$beta
gamma = as.list(coef(estimate5_3))$gamma
dt    = 1
m     = 3
n     = 3

# local functions
G <- function(l1) { Ln/(1+(beta*l1/Ln)^m) }
H <- function(l1) { sigma/(1+(gamma*l1/Ln)^n) }
lambda <- function(l1) { G(l1)/(H(l1))^2 }
ro <- function(l1) { (G(l1))^2/(H(l1))^2 }
k <- function(dl,l1) { ((l1+dl)*l1^(-exp(-g*dt)))^(1/(1-exp(-g*dt))) }

pLinf <- function(l1,simul) { rgamma(simul, rate=(G(l1)/H(l1)^2), shape=G(l1)^2/H(l1)^2)+l1 }

#calcola deltaL dalla 3 di Bardos et al. 2005
deltaL <- function(a,l1,dt) { (a*(l1/a)^exp(-g*dt))-l1 }

Trans_matrix <- function(limits) {
  p = matrix(0,nrow=length(limits)-1,ncol=length(limits)-1);

  ##OK!!
  for(i in 1:(length(limits)-1)) {
	for (j in 1:(length(limits)-1))	{
	  integrand <- function(l1) { 
	    pgamma(k(limits[i+1]-l1,l1)-l1, rate=lambda(l1), shape=ro(l1)) - 
		  pgamma(k(limits[i]-l1,l1)-l1, rate=lambda(l1), shape=ro(l1));
	  }
      p[i,j] = round((integrate(integrand,lower=limits[j],upper=limits[j+1]))$value/(limits[j+1]-limits[j]), 4);
	}
  }
  p[length(limits)-1,length(limits)-1] = 1;
  return(p);
}
