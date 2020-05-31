#### WITHIN-HOST MODEL ####  -----------------------------------------------------------------------------
library("here")

source(here("code/model_functions.R"))

# Include event definition that represents pathogens being cleared from the host. Otherwise, an infected host will always remain infected,
# and host immunity will reach equilibrium and never decline. This condition is assessed at every time step.
eventfun <- function(t, y, parms) {
  with (as.list(y), {
    P1 <- ifelse(P1 < 0.01, 0, P1)
    P2 <- ifelse(P2 < 0.01, 0, P2)

    return(c(P1, P2, I1, I2, M1, M2, R))
  })
}

# Series of equations that describe within-host immune dynamics that is used by the ODE solver.
deriv <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dR <- theta * (1 - (R/Rmax)) - (delta1 * P1 * R) - (delta2 * P2 * R)

    dP1 <- (alpha1 * delta1 * P1 * R) - (beta1 * I1 * P1)
    dP2 <- (alpha2 * delta2 * P2 * R) - (beta2 * I2 * P2)
    
    # dP1 <- (alpha1 * delta1 * P1 * R) - (beta1 * I1 * P1) - (beta2 * I2 * P1)  # General immunity scenario.
    # dP2 <- (alpha2 * delta2 * P2 * R) - (beta2 * I2 * P2) - (beta1 * I1 * P2)  # General immunity scenario.

    dI1 <- (epsilon * beta1 * I1 * P1) + (sigma1 * M1 * P1) + (tau * M2 * P1) - (gamma * I1)
    dI2 <- (epsilon * beta2 * I2 * P2) + (sigma2 * M2 * P2) + (tau * M1 * P2) - (gamma * I2)
    
    # dI1 <- (epsilon * beta1 * I1 * P1) + (epsilon * beta2 * I2 * P1) + (sigma1 * M1 * P1) + (sigma2 * M2 * P1) + (tau * M2 * P1) + (tau * M1 * P1) - (gamma * I1) # General immunity scenario.
    # dI2 <- (epsilon * beta2 * I2 * P2) + (epsilon * beta1 * I1 * P2) + (sigma2 * M2 * P2) + (sigma1 * M1 * P2) + (tau * M1 * P2) + (tau * M2 * P2) - (gamma * I2) # General immunity scenario.
    
    dM1 <- (rho * I1) - (mu * M1)
    dM2 <- (rho * I2) - (mu * M2)
    
    # dM1 <- (rho * I1) + (rho * I2) - (mu * M1)   # General immunity scenario.
    # dM2 <- (rho * I2) + (rho * I1) - (mu * M2)   # General immunity scenario.
    
    return(list(c(P1 = dP1, P2 = dP2, I1 = dI1, I2 = dI2, M1 = dM1, M2 = dM2, R = dR)))
  })
}

# Function that resolves within-host dynamics.
withinHost <- function(temp, 
                       N, 
                       times,
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
                       tau,
                       gamma,
                       rho,
                       mu) {
  for (i in 1:N) {
    parms <- c("theta" = theta,
               "Rmax" = Rmax,
               "delta1" = delta1,
               "delta2" = delta2,
               "alpha1" = alpha1, 
               "alpha2" = alpha2, 
               "beta1" = beta1, 
               "beta2" = beta2, 
               "epsilon" = epsilon, 
               "sigma1" = sigma1,
               "sigma2" = sigma2,
               "tau" = tau, 
               "gamma" = gamma, 
               "rho" = rho, 
               "mu" = mu)
    init <- c(P1 = temp[i, "P1"], 
              P2 = temp[i, "P2"], 
              I1 = temp[i, "I1"], 
              I2 = temp[i, "I2"], 
              M1 = temp[i, "M1"], 
              M2 = temp[i, "M2"],
              R = temp[i, "R"])
    
    wnhost_results <- lsoda(y = init, times = times, func = deriv, parms = parms, events = list(func = eventfun,
                                                                                                time = times))
    temp[i,] <- wnhost_results[2, c("P1", "P2", "I1", "I2", "M1", "M2", "R")]
    # The output of lsoda is a matrix with the first column = time variable and rows are times. Only take results from the simulated step forward (i.e. second row).
  }
  
  return(temp)
  
}

