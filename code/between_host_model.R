#### BETWEEN-HOST MODEL ####  -----------------------------------------------------------------------------

source("model_functions.R")
source("within_host_model.R")


# IBM for one time step.
betweenHost <- function(host.statuses.input, 
                        host.contacts, 
                        coef.logit,
                        host.trans.probs,
                        N,
                        times,
                        omega1,
                        omega2, 
                        bottle, 
                        eta,
                        theta,
                        Rmax,
                        delta1,
                        delta2,
                        alpha1,
                        alpha2,
                        beta1,
                        beta2,
                        epsilon,
                        sigma1,
                        sigma2,
                        gamma,
                        rho,
                        tau,
                        mu)  {
  temp <- host.statuses.input
    
  ### 1. DETERMINE PATHOGEN TRANSMISSION ###
  
  for (i in 1:nrow(host.contacts)) {
    contacting <- host.contacts[i,1]    # Identify the contacting host.
    contacted <- host.contacts[i,2]     # Identify the contacted host.
    
    path.load.p1 <- temp[contacting, "P1"]
    path.load.p2 <- temp[contacting, "P2"]
    
    
    # Determine if transmission occurs for P1.
    if (host.trans.probs[i,"P1"] <= logitToProb(coef.logit$P1, path.load.p1)) {
      
      # Calculate the amount of pathogen transmitted to contacted host.
      infecting.amount.p1 <- sum(rbinom(n = temp[contacting, "P1"], size = 1, prob = bottle))
      # Remove transmitted pathogen from contacting host. 
      temp[contacting, "P1"] <- temp[contacting, "P1"] - infecting.amount.p1
      # Add transmitted pathogen to contacted host.
      temp[contacted, "P1"] <- temp[contacted, "P1"] + infecting.amount.p1

    }
    
    
    # Determine if transmission occurs for P2.
    if (host.trans.probs[i,"P2"] <= logitToProb(coef.logit$P2, path.load.p2)) {
      
      # Calculate the amount of pathogen transmitted to contacted host.
      infecting.amount.p2 <- sum(rbinom(n = temp[contacting, "P2"], size = 1, prob = bottle))
      # Remove transmitted pathogen from contacting host. 
      temp[contacting, "P2"] <- temp[contacting, "P2"] - infecting.amount.p2
      # Add transmitted pathogen to contacted host.
      temp[contacted, "P2"] <- temp[contacted, "P2"] + infecting.amount.p2
      
    }
    
  
  }
  
  
  ### 2. RESOLVE WITHIN-HOST DYNAMICS ###
  
  # Resolve within-host dynamics.
  temp <- withinHost(temp = temp, 
                     N = N, 
                     times = times,
                     theta = theta,
                     Rmax = Rmax,
                     delta1 = delta1,
                     delta2 = delta2,
                     alpha1 = alpha1,
                     alpha2 = alpha2,
                     beta1 = beta1,
                     beta2 = beta2,
                     epsilon = epsilon,
                     sigma1 = sigma1,
                     sigma2 = sigma2,
                     gamma = gamma,
                     rho = rho,
                     tau = tau,
                     mu = mu)
  

  ### 3. RESOLVE DEMOGRAPHICS ###
  
  # Randomly select hosts in the population that will "die" and then replace their within-host statuses.
  for (i in 1:N) {
    if(rbinom(1, 1, prob = eta) == 1) {temp[i,] = c(0,0,1,1,0,0,Rmax)}
  }
  
  ####
  
  host.statuses.output <- temp
  
  return(host.statuses.output)
}

