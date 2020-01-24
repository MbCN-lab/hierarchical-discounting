###################################################################################
## This is code providing the simulation code for simulating and fitting Bayesian 
## hyperbolic discounting model from 
## "Hierarchies improve individual assessment of temporal discounting behavior" 
## by Molloy et al. (in preparation)
###################################################################################

# Load in required packages 


library("rjags")

########################
## Objects in the order they appear ##
# trials.seq is a vector containing the lengths of trials to simulate
# sim.k and sim.m is a matrix of the simulated k and m values by the length of trials.seq
# est.k and est.m are arrays where x= length of sim.k, y =  length of sim.m, and z= length of trials.seq containing means of model fit posteriors for k and m respectively
# sd.k and sd.m are arrays where x= length of sim.k, y =  length of sim.m, and z= length of trials.seq containing standard deviations of model fit posteriors for k and m respectively
# sim.trials gives the integer of trials to simulate
# sim.r is a vector of length sim.trials containing simulated delayed values
# v1s is a vector of possible initial values
# sim.V1 is a vector of length sim.trials containing simulated initial values
# delays is a vector of possible delay durations (in days) 
# sim.time is a vector of length sim.trials containing simulated delays
# sim.hyper.out is sim.trials by 2 matrix containing the calculated value of delay and simulated choice 
# mean_pd is a vector of length subjects calculating the overall probability of choosing larger-later
# k_hyper and m_hyper are the estimated hyperparameters governing the means of k and m respectively
# k_sigma and m_sigma are the estimated hyperparameters governing the precisions of k and m respectively
#####################


trials.seq=c(30,50,70,150)

sim.k=sim.m=matrix(NA,nrow=100,ncol=length(trials.seq))
est.k=est.m=array(NA,c(100,100,length(trials.seq)))
sd.k=sd.m=array(NA,c(100,100,length(trials.seq)))

for(x in 1:length(trials.seq)){
  sim.trials=trials.seq[x]
  
  #larger-later reward always 50 in Finn et al.
  sim.r=rep(50,sim.trials)
  
  #initial rewards
  v1s=seq(2.5,47.5,2.5)
  sim.V1=sample(v1s,sim.trials,replace=TRUE)
  
  #delay
  delays=c(7,14,30,90,365)
  sim.time=sample(delays,sim.trials,replace=TRUE)
  
  #sequence of true k and m
  sim.k[,x]=seq(0.01,1,0.01)
  sim.m[,x]=seq(0.05,5,0.05)
  
  for(k in 1:length(sim.k[,x])){
    for(m in 1:length(sim.m[,x])){
      
      print(paste("trials:",trials.seq[x],"k,m",k,m))
      sim.hyper.out=matrix(NA,nrow=sim.trials,ncol=2)
      colnames(sim.hyper.out)<-c("sim.Vd","sim.ch")
      
      sim.hyper <- function(sim.m,sim.trials,sim.V1,sim.r,sim.k,sim.time){
        for(t in 1:sim.trials){
          sim.hyper.out[t,"sim.Vd"] = sim.r[t]/(1+sim.k*sim.time[t])
          muPd = 1/(1+exp(-sim.m*(sim.hyper.out[t,"sim.Vd"]-sim.V1[t])))
          sim.hyper.out[t,"sim.ch"] = rbinom(1,1,muPd)
        }
        return(sim.hyper.out)
      }
      
      #Simulate choice data to pass into model
      sim.hyper.out=sim.hyper(sim.m[m],sim.trials,sim.V1,sim.r,sim.k[k],sim.time)
      
      #JAGS model code
      model.hyperbolic = "
    model{
    
    # likelihood
    for(t in 1:trials){
    Pd[t] ~ dbern(muPd[t]*.99999999+0.00000001) 
    muPd[t] = 1/(1+exp(-m*(Vd[t]-V1[t]))) 
    
    }
    
    # hyperbolic discounting function
    for(t in 1:trials){
    Vd[t] = r[t]/(1+k*time[t])
    }
    
    # Prior
    k ~ dunif(0,10)
    m ~ dunif(0,10)
    
    }
    "
      
      dat.hyperbolic = list(Pd=sim.hyper.out[,"sim.ch"],time=sim.time,V1=sim.V1,r=sim.r,trials=sim.trials)
      
      # Initialization
      model.hyperbolic = jags.model(textConnection(model.hyperbolic), data = dat.hyperbolic,
                                    n.chains = 3, n.adapt = 3000)
      
      # Burn-in
      update(model.hyperbolic, n.iter = 4000, progress.bar = "text") 
      
      # Posterior sampling
      hyperbolic.out = coda.samples(model = model.hyperbolic,
                                    variable.names = c("k","m"), 
                                    n.iter = 6000) 
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
      
      post=post.samples(hyperbolic.out)
      
      # Visually check chains
      # plot.ts(post[,"k",],plot.type=c("single"), main = paste("k",i), col=c("red","blue","green"))
      # plot.ts(post[,"m",],plot.type=c("single"), main = paste("m",i), col=c("red","blue","green"))
      
      # Save out means and standard deviations of posteriors
      est.k[k,m,x]=mean(post[,"k",])
      est.m[k,m,x]=mean(post[,"m",])
      sd.k[k,m,x]=sd(post[,"k",])
      sd.m[k,m,x]=sd(post[,"m",])
    }
    
  }
  
  # Plot simulated vs estimated k and m
  plot(rep(sim.k[,x],100),est.k[,,x],xlab="Simulated k",ylab="Estimated k",main=paste("Trials =",trials.seq[x]))
  abline(0,1)
  plot(rep(sim.m[,x],100),est.m[,,x],xlab="Simulated m",ylab="Estimated m",main=paste("Trials =",trials.seq[x]))
  abline(0,1)
  
}
