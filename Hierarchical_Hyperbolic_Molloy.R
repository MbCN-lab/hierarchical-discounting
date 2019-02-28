###################################################################################
## This is code providing the hierarchical Bayesian hyperbolic discounting model from 
## "Hierarchies improve individual assessment of temporal discounting behavior" 
## by Molloy et al. (in preparation)
###################################################################################

# Load in required packages 
library("rjags")

# #######################
# # Objects in the order they appear ##
# subjects is the number of subjects
# trials is a vector of length subjects containing the number of trials each subject completed
# Pd is a matrix of (maximum) trials by subjects containing the observed choice
# muPd is a matrix of (maximum) trials by subjects containing the probability of choosing a larger-later option
# m is a vector of length subjects containing each subject's m estimate
# Vd is a matrix of (maximum) trials by subjects containing the subjective (discounted) value of the larger-later option
# V1 is a matrix of (maximum) trials by subjects containing the subjective (discounted) value of the smaller-sooner option, if delay = 0, this is equivalent to the objective value
# r is a matrix of (maximum) trials by subjects containing the objective value of larger-later option
# k is a vector of length subjects containing the subject's estimated discounting rate (k)
# time is a matrix of (maximum) trials by subjects containing the delay in days of the larger-later option
# mean_pd is a vector of length subjects calculating the overall probability of choosing larger-later
# k_hyper and m_hyper are the estimated hyperparameters governing the means of k and m respectively
# k_sigma and m_sigma are the estimated hyperparameters governing the precisions of k and m respectively
# ####################

# JAGS Model 
model.hyperbolic = "
model{

# Likelihood
for (s in 1:subjects){
for(t in 1:trials[s]){
Pd[t,s] ~ dbern(muPd[t,s]*.99999999+0.00000001) 
muPd[t,s] = 1/(1+exp(-m[s]*(Vd[t,s]-V1[t,s]))) 
}
}

# Hyperbolic discounting function
for(s in 1:subjects){
for(t in 1:trials[s]){
Vd[t,s] = r[t,s]/(1+k[s]*time[t,s])
}
}

# Overall probability of choosing larger-later 
mean_pd[s]=sum(muPd[1:trials[s],s])/trials[s]

# Priors
for(s in 1:subjects){
k[s] ~ dnorm(k_hyper,sigma_k) T(0,)
m[s] ~ dnorm(m_hyper,sigma_m) T(0,)
}

k_hyper ~ dunif(0,10)
m_hyper ~ dunif(0,10)

sigma_k ~ dexp(1)
sigma_m ~ dexp(1)

}
"

# Save observed data as a list
dat.hyperbolic = list(Pd=Pd,time=time,V1=V1,r=r,trials=trials,subjects=subjects)

# Initialization
model.hyperbolic = jags.model(textConnection(model.hyperbolic), data = dat.hyperbolic,
                              n.chains = 3, n.adapt = 1500)

# Burn-in
update(model.hyperbolic, n.iter = 4000, progress.bar = "text") 

# Posterior sampling
hyperbolic.km = coda.samples(model = model.hyperbolic,
                             variable.names = c("k_hyper","m_hyper","k","m","muPd"),  # Choose which variables to save out
                             n.iter = 5000) 

# Converts JAGS object into array for convenience 
post.samples = function(out){
  n.chains = length(out)
  nrow = nrow(out[[1]])
  ncol = ncol(out[[1]])
  post = array(NA, c(nrow, ncol, n.chains))
  for (i in 1:n.chains){
    post[,,i] = out[[i]]
  }
  colnames(post) = colnames(out[[1]])
  return(post)
}

post_km=post.samples(hyperbolic.km)

