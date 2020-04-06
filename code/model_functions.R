#### FUNCTIONS ####

# Function to draw parameters for between-host and within-host models.
generateParams <- function() {  
  ##### SENSITIVITY ANALYSIS #####
  N <- 1000          
  t <- 2000           
  
  sd.var <- 0.25  # Standard deviation of the parameters when conducting sensitivity analysis.
  
  omega1 <- rnorm(1, mean = 0.01, sd = 0.01 * sd.var)
  omega2 <- rnorm(1, mean = 0.01, sd = 0.01 * sd.var) 
  kappa <- rnorm(1, mean = 2, sd = 2 * sd.var)
  bottle <- rnorm(1, mean = 0.05, sd = 0.05 * sd.var)
  eta <- rnorm(1, mean = 0.001, sd = 0.001 * sd.var)
  
  lo.n <- 0
  med.n <- 5000
  hi.n <- 10000
  
  lo.k.p1 <- 0
  lo.k.p2 <- 0
  med.k.p1 <- 0.025
  med.k.p2 <- 0.025
  hi.k.p1 <- 0.05
  hi.k.p2 <- 0.05
  
  ## Within-host model parameters ##
  times <- seq(0, 1, by = 1)  # The number of time steps conducted in the within-host model for every time step in the between-host model.
  
  theta <- rnorm(1, mean = 10, sd = 10 * sd.var)
  Rmax <- rnorm(1, mean = 100, sd = 100 * sd.var)
  delta1 <- rnorm(1, mean = 0.001, sd = 0.001 * sd.var)
  delta2 <- rnorm(1, mean = 0.001, sd = 0.001 * sd.var)
  alpha1 <- rnorm(1, mean = 100, sd = 100 * sd.var)
  alpha2 <- rnorm(1, mean = 100, sd = 100 * sd.var)
  beta1 <- rnorm(1, mean = 0.01, sd = 0.01 * sd.var)
  beta2 <- rnorm(1, mean = 0.01, sd = 0.01 * sd.var)
  epsilon <- rnorm(1, mean = 0.00001, sd = 0.00001 * sd.var)
  sigma1 <- rnorm(1, mean = 0.005, sd = 0.005 * sd.var)
  sigma2 <- rnorm(1, mean = 0.005, sd = 0.005 * sd.var)
  tau <- rnorm(1, mean = 0.0025, sd = 0.0025 * sd.var)
  gamma <- rnorm(1, mean = 0.05, sd = 0.05 * sd.var)
  rho <- rnorm(1, mean = 0.0005, sd = 0.0005 * sd.var)
  mu <- rnorm(1, mean = 0.0005, sd = 0.0005 * sd.var)
  
  # Coexistence decomposition parameters.
  equil.times <- c(1900:2000) # Timeframe of equilibrium.
  spat.equil.times <- 10    # Amount of time given to reach spatial equilibrium.
  
  

  parameter.list <- list("N" = N,
                         "t" = t,
                         "sd.var" = sd.var,
                         "omega1" = omega1,
                         "omega2" = omega2,
                         "kappa" = kappa,
                         "bottle" = bottle,
                         "eta" = eta,
                         "lo.n" = lo.n,
                         "med.n" = med.n,
                         "hi.n" = hi.n,
                         "lo.k.p1" = lo.k.p1,
                         "lo.k.p2" = lo.k.p2,
                         "med.k.p1" = med.k.p1,
                         "med.k.p2" = med.k.p2,
                         "hi.k.p1" = hi.k.p1,
                         "hi.k.p2" = hi.k.p2,
                         "times" = times,
                         "theta" = theta,
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
                         "gamma" = gamma,
                         "rho" = rho,
                         "tau" = tau,
                         "mu" = mu,
                         "equil.times" = equil.times,
                         "spat.equil.times" = spat.equil.times)
  
  return(parameter.list)
                         
}


# Function that utilizes initial parameters and enters them into array
# initP1 = initial number of hosts infected with P1
# initP1.load = initial load of P1 in P1 infected hosts
# initP2 = initial number of hosts infected with P2
# initP2.load = initial load of P2 in P2 infected hosts
initialCondit <- function(results, initP1, initP1.load, initP2, initP2.load) { 
  temp <- results[,,1]
  temp[1:initP1,"P1"] <- ifelse(initP1 != 0, initP1.load, 0)
  temp[(initP1 + 1):(initP1 + initP2),"P2"] <- ifelse(initP2 != 0, initP2.load, 0)
  
  return(temp)
  
}


# Function to generate logit function to determine probability of between-host transmission.
# XX.k represents probability of transmission at XX pathogen load, while XX.n represents that pathogen load 
calcLogit <- function(lo.k, lo.n, med.k, med.n, hi.k, hi.n) {
  
  a.k <- rbinom(n = 10000, size = 1, prob = lo.k)
  a.n <- rep(x = lo.n, times = 10000)
  
  b.k <- rbinom(n = 10000, size = 1, prob = med.k)
  b.n <- rep(x = med.n, times = 10000)
  
  c.k <- rbinom(10000, size = 1, prob = hi.k)
  c.n <- rep(x = hi.n, times = 10000)
  
  temp1 <- append(a.k, values = c(b.k, c.k))
  temp2 <- append(a.n, values = c(b.n, c.n))
  
  temp3 <- data.frame("trans" = temp1, "load" = temp2)
  
  logit <- glm(trans ~ load, family = binomial(link = "logit"), data = temp3)
  coef.logit <- coef(logit)
  
  return(coef.logit)
  
}


# Function to convert log odds to probability for between-host transmission.
logitToProb <- function(coef.logit, path.load) {
  odds <- exp(coef.logit[1] + coef.logit[2] * path.load)
  prob <- odds / (1 + odds)
  return(prob)
}
    

#############################################################################################

#### DECOMPOSITION FUNCTIONS ####

# Function to determine the spatial equilibrium of invading pathogen and then to calculate Nt+1 of the pathogen.
calculateLDGR <- function(resident.equilibrium,
                         host.status,
                         status.categ,
                         invader, 
                         equil.times, 
                         spat.equil.times, 
                         iabun,
                         invade.abundance,
                         t,
                         host.contacts.eq,
                         coef.logit,
                         host.trans.probs.eq,
                         N,
                         times,
                         bottle, 
                         eta,
                         omega1,
                         omega2,
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
                         mu) {
  
  space.equil <- array(NA, dim = c(N, 1, spat.equil.times), dimnames = list(c(seq(N)), invader, seq(spat.equil.times)))
  space.equil[,,1] <- invade.abundance
  
  spat.equil.check <- array(NA, dim = c(N, spat.equil.times, length(equil.times)))
  
  resident.invader.t0 <- array(host.status, dim = c(N, length(status.categ), length(equil.times)), 
                               dimnames = list(c(seq(N)), status.categ, c(1:(length(equil.times)))))
  resident.invader.t1 <- array(host.status, dim = c(N, length(status.categ), length(equil.times)), 
                               dimnames = list(c(seq(N)), status.categ, c(1:(length(equil.times)))))
  
  
  for (y in 1:length(equil.times)) {
    invade.abundance.t0 <- invade.abundance
    
    # Calculate spatial equilibrium for current time.
    for (r in 1:(spat.equil.times-1)) {
      temp.space.equib <- resident.equilibrium[,,y]
      temp.space.equib[,invader] <- invade.abundance.t0
      
      model.outcome <- betweenHost(host.statuses.input = temp.space.equib, 
                                   host.contacts = host.contacts.eq[[y]], 
                                   coef.logit = coef.logit,
                                   host.trans.probs = host.trans.probs.eq[[y]],
                                   N = N,
                                   times = times,
                                   omega1 = omega1,
                                   omega2 = omega2, 
                                   bottle = bottle, 
                                   eta = eta,
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
      
      
      invade.abundance.t1 <- model.outcome[,invader] / sum(model.outcome[,invader]) * iabun
      
      space.equil[,,r+1] <- invade.abundance.t1
      spat.equil.check[,r,y] <- round(((invade.abundance.t1 - invade.abundance.t0) * 100 / invade.abundance.t0), digits = 3)
      
      invade.abundance.t0 <- invade.abundance.t1
      
    }
    
    # Using the invading abundance at spatial equilibrium, calculate one more time step.
    resident.invader.t0[,,y] <- resident.equilibrium[,,y]
    resident.invader.t0[,invader,y] <- space.equil[,,spat.equil.times]
    resident.invader.t1[,,y] <- betweenHost(host.statuses.input = resident.invader.t0[,,y], 
                                            host.contacts = host.contacts.eq[[y]], 
                                            coef.logit = coef.logit,
                                            host.trans.probs = host.trans.probs.eq[[y]],
                                            N = N,
                                            times = times,
                                            omega1 = omega1,
                                            omega2 = omega2, 
                                            bottle = bottle, 
                                            eta = eta,
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
    
    
    # if(y %% 1 == 0) {
    #   print(paste(y, "of", length(equil.times), "finished"))
    # }
    
  }
  
  return.list <- list("t0" = resident.invader.t0, "t1" = resident.invader.t1, "check" = spat.equil.check)
  return(return.list)
  
}

#########################

calculateCoexist <- function(resident.equilibrium,
                             host.status,
                             status.categ,
                             equil.times, 
                             spat.equil.times, 
                             iabun,
                             invade.abundance,
                             t,
                             host.contacts.eq,
                             coef.logit,
                             host.trans.probs.eq,
                             N,
                             times,
                             bottle, 
                             eta,
                             omega1,
                             omega2,
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
                             mu) {
  
  resident.invader.t0 <- array(host.status, dim = c(N, length(status.categ), length(equil.times)), 
                               dimnames = list(c(seq(N)), status.categ, c(1:(length(equil.times)))))
  resident.invader.t1 <- array(host.status, dim = c(N, length(status.categ), length(equil.times)), 
                               dimnames = list(c(seq(N)), status.categ, c(1:(length(equil.times)))))
  
  
  for (y in 1:length(equil.times)) {
    # Using the invading abundance at spatial equilibrium, calculate one more time step.
    resident.invader.t0[,,y] <- resident.equilibrium[,,y]
    resident.invader.t1[,,y] <- betweenHost(host.statuses.input = resident.invader.t0[,,y], 
                                            host.contacts = host.contacts.eq[[y]], 
                                            coef.logit = coef.logit,
                                            host.trans.probs = host.trans.probs.eq[[y]],
                                            N = N,
                                            times = times,
                                            omega1 = omega1,
                                            omega2 = omega2, 
                                            bottle = bottle, 
                                            eta = eta,
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
    
    
    if(y %% 1 == 0) {
      print(paste(y, "of", length(equil.times), "finished"))
    }
    
  }
  
  return.list <- list("t0" = resident.invader.t0, "t1" = resident.invader.t1)
  return(return.list)
  
}


# Function to calculate growth rate.
calculateR <- function(pop.t0, pop.t1) {
  tempt0 <- apply(pop.t0, c(2,3), sum)
  tempt1 <- apply(pop.t1, c(2,3), sum)
  
  tempR <- log(tempt1/tempt0)
  R.bar.result <- apply(tempR, 1, mean)

  return(R.bar.result)
}



#############################################################################################

#### DIAGNOSTIC FUNCTIONS #####


# Function that graphs  model results.
graphBetween <- function(t, results, i) {
  numSus <- numP1 <- numP2 <- numCoI <- c()
  times <- c(1:t)
  
  for (a in 1:t) {
    numP1[a] <- sum(results[,"P1",a] > 100)
    numP2[a] <- sum(results[,"P2",a] > 100)
    numSus[a] <- sum(results[,"P1",a] == 0 & results[,"P2",a] == 0)
    numCoI[a] <- sum(results[,"P1",a] > 0 & results[,"P2",a] > 100)
}
  
  numPop <- as.data.frame(cbind(numSus, numP1, numP2, times))
  
  ggplot(numPop, aes(x = times)) +
    geom_line(aes(y = numSus), color = "blue", size = 1) +
    geom_line(aes(y = numP1), color = "red", size = 1) +
    geom_line(aes(y = numP2), color = "green", size = 1) +
    geom_line(aes(y = numCoI), color = "black", size = 1) +
    # ggtitle(paste("SIMULATION #", i, sep = "")) +
    theme_classic()
  
}

graphWithin <- function(t, results) {
  totP1 <- totP2 <- totI1 <- totI2 <- totM1 <- totM2 <- totR <- c()
  times <- c(1:t)
  
  for (a in 1:t) {
    totP1[a] <- sum(results[,"P1",a])
    totP2[a] <- sum(results[,"P2",a])
    totI1[a] <- sum(results[,"I1",a])
    totI2[a] <- sum(results[,"I2",a])
    totM1[a] <- sum(results[,"M1",a])
    totM2[a] <- sum(results[,"M2",a])
    totR[a] <- sum(results[,"R",a])
  }
  
  totWithin <- as.data.frame(cbind(totP1, totP2, totI1, totI2, totM1, totM2, totR, times))
  
  ggplot(totWithin, aes(x = times)) +
    geom_line(aes(y = totP1), color = "red", size = 1) +
    geom_line(aes(y = totI1), color = "green", size = 1) +
    geom_line(aes(y = totM1), color = "blue", size = 1) +
    geom_line(aes(y = totP2), color = "red", linetype = 2, size = 1) +
    geom_line(aes(y = totI2), color = "green", linetype = 2, size = 1) +
    geom_line(aes(y = totM2), color = "blue", linetype = 2, size = 1) +
    # geom_line(aes(y = totR), color = "black", linetype = 1, size = 1) +
    # ggtitle(paste("SIMULATION #", i, sep = "")) +
    theme_classic()
  
}

wnHostDiag <- function(N, host, t, results) {
  times <- c(1:t)
  
  diagHost <- t(results[host,,times])
  diagHost <- as.data.frame(diagHost)
  
  ggplot(diagHost, aes(x = times)) +
    geom_line(aes(y = P1), color = "red", size = 1) +
    geom_line(aes(y = I1), color = "green", size = 1) +
    geom_line(aes(y = M1), color = "blue", size = 1) +
    geom_line(aes(y = P2), color = "red", linetype = 2, size = 1) +
    geom_line(aes(y = I2), color = "green", linetype = 2, size = 1) +
    geom_line(aes(y = M2), color = "blue", linetype = 2, size = 1) +
    geom_line(aes(y = R), color = "black", linetype = 1, size = 1) +
    ggtitle(paste("HOST #", host, sep = "")) +
    theme_classic()
}
  
    
