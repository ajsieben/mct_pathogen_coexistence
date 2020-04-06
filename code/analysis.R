library("ggplot2")
library("gridExtra")
library("tidyverse")
library("reshape2")

# setwd("./../..")
source("model_functions.R")


##### FIGURE 2: SENSITIVITY ANALYSIS ######

# tot.p1.p2.resident.equilibrium <- 
#   tot.p1.resident <- 
#   tot.p2.resident <- 
#   list()

tot.parameter.list <- 
  tot.p1.invader.decomp.indiv <- 
  tot.p1.invader.decomp.compare <- 
  tot.p2.invader.decomp.indiv <-
  tot.p2.invader.decomp.compare <- 
  list()

numSims <- 511
setwd("./results/sens_analysis")


for (b in 1:numSims) {
  if (file.exists(paste("decomp_results_", b, ".RData", sep = "")) == TRUE) {
    # load(paste("raw_results_", b, ".RData", sep = ""))
    load(paste("decomp_results_", b, ".RData", sep = ""))

    # tot.p1.p2.resident.equilibrium[[b]] <- raw.results$results.p1.p2.equilibrium
    # tot.p1.resident[[b]] <- raw.results$results.p1.resident
    # tot.p1.resident.no.var[[b]] <- raw.results$results.p1.resident.no.var
    # tot.p1.resident.kappa.var[[b]] <- raw.results$results.p1.resident.kappa.var
    # tot.p1.resident.sigma.var[[b]] <- raw.results$results.p1.resident.sigma.var
    # tot.p2.resident[[b]] <- raw.results$results.p2.resident
    # tot.p2.resident.no.var[[b]] <- raw.results$results.p2.resident.no.var
    # tot.p2.resident.kappa.var[[b]] <- raw.results$results.p2.resident.kappa.var
    # tot.p2.resident.sigma.var[[b]] <- raw.results$results.p2.resident.sigma.var
    
    decomp.results$parameter.list$equil.times <- length(decomp.results$parameter.list$equil.times)
    tot.parameter.list[[b]] <- decomp.results$parameter.list
    tot.p1.invader.decomp.indiv[[b]] <- decomp.results$p1.invader.decomp.indiv
    tot.p1.invader.decomp.compare[[b]] <- decomp.results$p1.invader.decomp.compare
    tot.p2.invader.decomp.indiv[[b]] <- decomp.results$p2.invader.decomp.indiv
    tot.p2.invader.decomp.compare[[b]] <- decomp.results$p2.invader.decomp.compare
    
  }
  
}

tot.parameter.list <- lapply(tot.parameter.list, unlist)
param.names <- append(names(tot.parameter.list[[1]]), names(tot.p1.invader.decomp.compare[[1]]))

sens.results <- matrix(nrow = numSims, ncol = length(param.names), dimnames = list(seq(numSims), param.names))

for (c in 1:numSims) {
  temp <- append(tot.parameter.list[[c]], tot.p1.invader.decomp.compare[[c]])
  ifelse(is.null(temp) == FALSE, sens.results[c,] <- temp, sens.results[c,] <- NA)
}

sens.results <- as.data.frame(sens.results)
sens.results <- na.omit(sens.results)


#####

fit <- glm(LDGR ~ theta +
                  Rmax +
                  avg.delta1 +
                  avg.delta2 +
                  alpha1 +
                  alpha2 +
                  beta1 +
                  beta2 +
                  epsilon +
                  avg.sigma1 +
                  avg.sigma2 +
                  cross +
                  gamma +
                  rho +
                  mu +
                  avg.omega1 +
                  avg.omega2 +
                  avg.contacts +
                  eta,
                  # eta.p1 +
                  # eta.p2,
                  
                  data = sens.results)
                  
                  

fit.coef <- summary(fit)$coef[-1,1]
sd.err <- summary(fit)$coef[-1,2]

sx <- as.numeric(lapply(fit$model[-1], sd))
sy <- as.numeric(lapply(fit$model[1], sd))
beta <- fit.coef * sx/sy
stand.sd.err <- sd.err * sx/sy
stand.sd.dev <- stand.sd.err * sqrt(nrow(sens.results))

params <- names(beta)
coefficients <- as.numeric(unname(beta))
fit.stand <- data.frame("params" = params, "coefficients" = coefficients, "se" = stand.sd.err)
fit.stand$params <- factor(fit.stand$params, levels = fit.stand$params)


ggplot(data = fit.stand, aes(x = params, y = coefficients)) + 
  geom_bar(stat = "identity", fill = "gray80", colour = "black") +
  geom_errorbar(aes(ymin = coefficients - se, ymax = coefficients + se), width = 0.2) +
  labs(title = "Sensitivity Analysis w/ SE") +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


setwd("./../..")
# write.csv(fit.stand, file = "sens_analysis.csv")


############################################################################################
############################################################################################

###### FIGURE 3: NEUTRAL SCENARIO #######

#### NEUTRALITY WITH GENERAL IMMUNITY #####

tot.parameter.list <- 
  tot.p1.invader.decomp.indiv <- 
  tot.p1.invader.decomp.compare <- 
  tot.p2.invader.decomp.indiv <-
  tot.p2.invader.decomp.compare <- 
  list()

numSims <- 1
setwd("./results/neutral_general")


for (b in 1:numSims) {
  if (file.exists(paste("decomp_results_", b, ".RData", sep = "")) == TRUE) {
    # load(paste("raw_results_", b, ".RData", sep = ""))
    load(paste("decomp_results_", b, ".RData", sep = ""))

    # tot.p1.p2.resident.equilibrium[[b]] <- raw.results$results.p1.p2.equilibrium
    # tot.p1.resident[[b]] <- raw.results$results.p1.resident
    # tot.p1.resident.no.var[[b]] <- raw.results$results.p1.resident.no.var
    # tot.p1.resident.kappa.var[[b]] <- raw.results$results.p1.resident.kappa.var
    # tot.p1.resident.sigma.var[[b]] <- raw.results$results.p1.resident.sigma.var
    # tot.p2.resident[[b]] <- raw.results$results.p2.resident
    # tot.p2.resident.no.var[[b]] <- raw.results$results.p2.resident.no.var
    # tot.p2.resident.kappa.var[[b]] <- raw.results$results.p2.resident.kappa.var
    # tot.p2.resident.sigma.var[[b]] <- raw.results$results.p2.resident.sigma.var
    
    decomp.results$parameter.list$equil.times <- length(decomp.results$parameter.list$equil.times)
    tot.parameter.list[[b]] <- decomp.results$parameter.list
    tot.p1.invader.decomp.indiv[[b]] <- decomp.results$p1.invader.decomp.indiv
    tot.p1.invader.decomp.compare[[b]] <- decomp.results$p1.invader.decomp.compare
    tot.p2.invader.decomp.indiv[[b]] <- decomp.results$p2.invader.decomp.indiv
    tot.p2.invader.decomp.compare[[b]] <- decomp.results$p2.invader.decomp.compare
    
  }
  
}

temp1.p1 <- as.data.frame(tot.p1.invader.decomp.compare)
temp1.p2 <- as.data.frame(tot.p2.invader.decomp.compare)

temp2.p1 <- data.frame(t(temp1.p1), row.names = NULL)
temp2.p2 <- data.frame(t(temp1.p2), row.names = NULL)

temp.mean.p1 <- unlist(lapply(temp2.p1, mean))
temp.mean.p2 <- unlist(lapply(temp2.p2, mean))


std.err <- function(x) sd(x)/sqrt(length(x))
temp.se.p1 <- unlist(lapply(temp2.p1, sd))
temp.se.p2 <- unlist(lapply(temp2.p2, sd))

ldgr.iden <- as.factor(c(1,0,0,0,0))
decomp <- c("LDGR", "e0", "eL", "eD", "eL*D")


neutral.general.decomp.p1 <- data.frame("decomp" = decomp, "mean" = temp.mean.p1, "se" = temp.se.p1, ldgr.iden)
neutral.general.decomp.p2 <- data.frame("decomp" = decomp, "mean" = temp.mean.p2, "se" = temp.se.p2, ldgr.iden)
neutral.general.decomp.p1$decomp <- factor(neutral.general.decomp.p1$decomp, levels = neutral.general.decomp.p1$decomp)
neutral.general.decomp.p2$decomp <- factor(neutral.general.decomp.p2$decomp, levels = neutral.general.decomp.p2$decomp)



setwd("./../..")

######################################################################

#### NEUTRALITY WITH SPECIFIC IMMUNITY #####

tot.parameter.list <- 
  tot.p1.invader.decomp.indiv <- 
  tot.p1.invader.decomp.compare <- 
  tot.p2.invader.decomp.indiv <-
  tot.p2.invader.decomp.compare <- 
  list()

numSims <- 1
setwd("./results/neutral_specific")


for (b in 1:numSims) {
  if (file.exists(paste("decomp_results_", b, ".RData", sep = "")) == TRUE) {
    # load(paste("raw_results_", b, ".RData", sep = ""))
    load(paste("decomp_results_", b, ".RData", sep = ""))
    
    # tot.p1.p2.resident.equilibrium[[b]] <- raw.results$results.p1.p2.equilibrium
    # tot.p1.resident[[b]] <- raw.results$results.p1.resident
    # tot.p1.resident.no.var[[b]] <- raw.results$results.p1.resident.no.var
    # tot.p1.resident.kappa.var[[b]] <- raw.results$results.p1.resident.kappa.var
    # tot.p1.resident.sigma.var[[b]] <- raw.results$results.p1.resident.sigma.var
    # tot.p2.resident[[b]] <- raw.results$results.p2.resident
    # tot.p2.resident.no.var[[b]] <- raw.results$results.p2.resident.no.var
    # tot.p2.resident.kappa.var[[b]] <- raw.results$results.p2.resident.kappa.var
    # tot.p2.resident.sigma.var[[b]] <- raw.results$results.p2.resident.sigma.var
    
    decomp.results$parameter.list$equil.times <- length(decomp.results$parameter.list$equil.times)
    tot.parameter.list[[b]] <- decomp.results$parameter.list
    tot.p1.invader.decomp.indiv[[b]] <- decomp.results$p1.invader.decomp.indiv
    tot.p1.invader.decomp.compare[[b]] <- decomp.results$p1.invader.decomp.compare
    tot.p2.invader.decomp.indiv[[b]] <- decomp.results$p2.invader.decomp.indiv
    tot.p2.invader.decomp.compare[[b]] <- decomp.results$p2.invader.decomp.compare
    
  }
  
}

temp1.p1 <- as.data.frame(tot.p1.invader.decomp.compare)
temp1.p2 <- as.data.frame(tot.p2.invader.decomp.compare)

temp2.p1 <- data.frame(t(temp1.p1), row.names = NULL)
temp2.p2 <- data.frame(t(temp1.p2), row.names = NULL)

temp.mean.p1 <- unlist(lapply(temp2.p1, mean))
temp.mean.p2 <- unlist(lapply(temp2.p2, mean))

temp.se.p1 <- unlist(lapply(temp2.p1, sd))
temp.se.p2 <- unlist(lapply(temp2.p2, sd))

neutral.specific.decomp.p1 <- data.frame("decomp" = decomp, "mean" = temp.mean.p1, "se" = temp.se.p1, ldgr.iden)
neutral.specific.decomp.p2 <- data.frame("decomp" = decomp, "mean" = temp.mean.p2, "se" = temp.se.p2, ldgr.iden)
neutral.specific.decomp.p1$decomp <- factor(neutral.specific.decomp.p1$decomp, levels = neutral.specific.decomp.p1$decomp)
neutral.specific.decomp.p2$decomp <- factor(neutral.specific.decomp.p2$decomp, levels = neutral.specific.decomp.p2$decomp)


setwd("./../..")

#################

neutral.general <- ggplot(data = neutral.general.decomp.p1, aes(x = decomp, y = mean, fill = ldgr.iden)) +
                          geom_bar(stat = "identity", position = "dodge", colour = "black") +
                          geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
                          scale_fill_manual(values = c("gray80", "gray40")) +
                          scale_y_continuous(limits = c(-3,9), breaks = c(-3, 0, 3, 6, 9)) +
                          theme_classic()


neutral.specific <- ggplot(data = neutral.specific.decomp.p1, aes(x = decomp, y = mean, fill = ldgr.iden)) +
                           geom_bar(stat = "identity", position = "dodge", colour = "black") +
                           geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
                           scale_fill_manual(values = c("gray80", "gray40")) +
                           scale_y_continuous(limits = c(-3,9), breaks = c(-3, 0, 3, 6, 9)) +
                           theme_classic()
  

grid.arrange(neutral.specific, neutral.general, ncol = 2)



############################################################################################
############################################################################################

###### FIGURE 4: COMPETITION/COLONIZATION SCENARIO #######

#### COMPETITION/COLONIZATION WITH GENERAL IMMUNITY #####

tot.parameter.list <- 
  tot.p1.invader.decomp.indiv <- 
  tot.p1.invader.decomp.compare <- 
  tot.p2.invader.decomp.indiv <-
  tot.p2.invader.decomp.compare <- 
  list()

numSims <- 1
setwd("./results/comp_col_general")


for (b in 1:numSims) {
  if (file.exists(paste("decomp_results_", b, ".RData", sep = "")) == TRUE) {
    # load(paste("raw_results_", b, ".RData", sep = ""))
    load(paste("decomp_results_", b, ".RData", sep = ""))
    
    # tot.p1.p2.resident.equilibrium[[b]] <- raw.results$results.p1.p2.equilibrium
    # tot.p1.resident[[b]] <- raw.results$results.p1.resident
    # tot.p1.resident.no.var[[b]] <- raw.results$results.p1.resident.no.var
    # tot.p1.resident.kappa.var[[b]] <- raw.results$results.p1.resident.kappa.var
    # tot.p1.resident.sigma.var[[b]] <- raw.results$results.p1.resident.sigma.var
    # tot.p2.resident[[b]] <- raw.results$results.p2.resident
    # tot.p2.resident.no.var[[b]] <- raw.results$results.p2.resident.no.var
    # tot.p2.resident.kappa.var[[b]] <- raw.results$results.p2.resident.kappa.var
    # tot.p2.resident.sigma.var[[b]] <- raw.results$results.p2.resident.sigma.var
    
    decomp.results$parameter.list$equil.times <- length(decomp.results$parameter.list$equil.times)
    tot.parameter.list[[b]] <- decomp.results$parameter.list
    tot.p1.invader.decomp.indiv[[b]] <- decomp.results$p1.invader.decomp.indiv
    tot.p1.invader.decomp.compare[[b]] <- decomp.results$p1.invader.decomp.compare
    tot.p2.invader.decomp.indiv[[b]] <- decomp.results$p2.invader.decomp.indiv
    tot.p2.invader.decomp.compare[[b]] <- decomp.results$p2.invader.decomp.compare
    
  }
  
}

temp1.p1 <- as.data.frame(tot.p1.invader.decomp.compare)
temp1.p2 <- as.data.frame(tot.p2.invader.decomp.compare)

temp2.p1 <- data.frame(t(temp1.p1), row.names = NULL)
temp2.p2 <- data.frame(t(temp1.p2), row.names = NULL)

temp.mean.p1 <- unlist(lapply(temp2.p1, mean))
temp.mean.p2 <- unlist(lapply(temp2.p2, mean))

temp.se.p1 <- unlist(lapply(temp2.p1, sd))
temp.se.p2 <- unlist(lapply(temp2.p2, sd))

ldgr.iden <- as.factor(c(1,0,0,0,0))
decomp <- c("LDGR", "e0", "eL", "eD", "eL*D")


comp.col.general.decomp.p1 <- data.frame("decomp" = decomp, "mean" = temp.mean.p1, "se" = temp.se.p1, ldgr.iden)
comp.col.general.decomp.p2 <- data.frame("decomp" = decomp, "mean" = temp.mean.p2, "se" = temp.se.p2, ldgr.iden)
comp.col.general.decomp.p1$decomp <- factor(comp.col.general.decomp.p1$decomp, levels = comp.col.general.decomp.p1$decomp)
comp.col.general.decomp.p2$decomp <- factor(comp.col.general.decomp.p2$decomp, levels = comp.col.general.decomp.p2$decomp)



setwd("./../..")

######################################################################

#### COMPETITION/COLONIZATION WITH SPECIFIC IMMUNITY #####

tot.parameter.list <- 
  tot.p1.invader.decomp.indiv <- 
  tot.p1.invader.decomp.compare <- 
  tot.p2.invader.decomp.indiv <-
  tot.p2.invader.decomp.compare <- 
  list()

numSims <- 1
setwd("./results/comp_col_specific")


for (b in 1:numSims) {
  if (file.exists(paste("decomp_results_", b, ".RData", sep = "")) == TRUE) {
    # load(paste("raw_results_", b, ".RData", sep = ""))
    load(paste("decomp_results_", b, ".RData", sep = ""))
    
    # tot.p1.p2.resident.equilibrium[[b]] <- raw.results$results.p1.p2.equilibrium
    # tot.p1.resident[[b]] <- raw.results$results.p1.resident
    # tot.p1.resident.no.var[[b]] <- raw.results$results.p1.resident.no.var
    # tot.p1.resident.kappa.var[[b]] <- raw.results$results.p1.resident.kappa.var
    # tot.p1.resident.sigma.var[[b]] <- raw.results$results.p1.resident.sigma.var
    # tot.p2.resident[[b]] <- raw.results$results.p2.resident
    # tot.p2.resident.no.var[[b]] <- raw.results$results.p2.resident.no.var
    # tot.p2.resident.kappa.var[[b]] <- raw.results$results.p2.resident.kappa.var
    # tot.p2.resident.sigma.var[[b]] <- raw.results$results.p2.resident.sigma.var
    
    decomp.results$parameter.list$equil.times <- length(decomp.results$parameter.list$equil.times)
    tot.parameter.list[[b]] <- decomp.results$parameter.list
    tot.p1.invader.decomp.indiv[[b]] <- decomp.results$p1.invader.decomp.indiv
    tot.p1.invader.decomp.compare[[b]] <- decomp.results$p1.invader.decomp.compare
    tot.p2.invader.decomp.indiv[[b]] <- decomp.results$p2.invader.decomp.indiv
    tot.p2.invader.decomp.compare[[b]] <- decomp.results$p2.invader.decomp.compare
    
  }
  
}

temp1.p1 <- as.data.frame(tot.p1.invader.decomp.compare)
temp1.p2 <- as.data.frame(tot.p2.invader.decomp.compare)

temp2.p1 <- data.frame(t(temp1.p1), row.names = NULL)
temp2.p2 <- data.frame(t(temp1.p2), row.names = NULL)

temp.mean.p1 <- unlist(lapply(temp2.p1, mean))
temp.mean.p2 <- unlist(lapply(temp2.p2, mean))

temp.se.p1 <- unlist(lapply(temp2.p1, sd))
temp.se.p2 <- unlist(lapply(temp2.p2, sd))

comp.col.specific.decomp.p1 <- data.frame("decomp" = decomp, "mean" = temp.mean.p1, "se" = temp.se.p1, ldgr.iden)
comp.col.specific.decomp.p2 <- data.frame("decomp" = decomp, "mean" = temp.mean.p2, "se" = temp.se.p2, ldgr.iden)
comp.col.specific.decomp.p1$decomp <- factor(comp.col.specific.decomp.p1$decomp, levels = comp.col.specific.decomp.p1$decomp)
comp.col.specific.decomp.p2$decomp <- factor(comp.col.specific.decomp.p2$decomp, levels = comp.col.specific.decomp.p2$decomp)


setwd("./../..")

#################

comp.col.general.p1 <- ggplot(data = comp.col.general.decomp.p1, aes(x = decomp, y = mean, fill = ldgr.iden)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  scale_fill_manual(values = c("gray80", "gray40")) +
  scale_y_continuous(limits = c(-8,11), breaks = c(-5, 0, 5, 10)) +
  theme_classic()

comp.col.general.p2 <- ggplot(data = comp.col.general.decomp.p2, aes(x = decomp, y = mean, fill = ldgr.iden)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  scale_fill_manual(values = c("gray80", "gray40")) +
  scale_y_continuous(limits = c(-8,11), breaks = c(-5, 0, 5, 10)) +
  theme_classic()



comp.col.specific.p1 <- ggplot(data = comp.col.specific.decomp.p1, aes(x = decomp, y = mean, fill = ldgr.iden)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  scale_fill_manual(values = c("gray80", "gray40")) +
  scale_y_continuous(limits = c(-8,11), breaks = c(-5, 0, 5, 10)) +
  theme_classic()

comp.col.specific.p2 <- ggplot(data = comp.col.specific.decomp.p2, aes(x = decomp, y = mean, fill = ldgr.iden)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  scale_fill_manual(values = c("gray80", "gray40")) +
  scale_y_continuous(limits = c(-8,11), breaks = c(-5, 0, 5, 10)) +
  theme_classic()


grid.arrange(comp.col.specific.p1, comp.col.specific.p2,
             comp.col.general.p1, comp.col.general.p2, ncol = 2)




#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################


##### BETWEEN HOST DIAGNOSTICS #####
# results.p1.p2.equilibrium <- raw.results$results.p1.p2.equilibrium
# temp <- dim(results.p1.p2.equilibrium)
# t <- temp[3]
# 
# graphBetween(t = t, results = results.p1.p2.equilibrium)


##### WITHIN INDIVIDUAL HOST DIAGNOSTICS ######

# diag.data <- raw.results$results.p1.p2.equilibrium
# 
# wnHostDiag(N = N,
#            host = sample(1:N, 1, replace = TRUE),
#            t = t,
#            results = diag.data)
# 
# 
# pops <- c("P1","P2","I1","I2","M1","M2")
# for (i in pops) {
#   hist(diag.data[,i,], main = paste(i))
# }


########################### SUPPLEMENTARY FIGURE: TIME STEP ANALYSIS ################################

#### TIME STEP 3 #####

tot.parameter.list <- 
  tot.p1.invader.decomp.indiv <- 
  tot.p1.invader.decomp.compare <- 
  tot.p2.invader.decomp.indiv <-
  tot.p2.invader.decomp.compare <- 
  list()

numSims <- 31
setwd("./results/time_step_3")


for (b in 1:numSims) {
  if (file.exists(paste("decomp_results_", b, ".RData", sep = "")) == TRUE) {
    # load(paste("raw_results_", b, ".RData", sep = ""))
    load(paste("decomp_results_", b, ".RData", sep = ""))
    
    # tot.p1.p2.resident.equilibrium[[b]] <- raw.results$results.p1.p2.equilibrium
    # tot.p1.resident[[b]] <- raw.results$results.p1.resident
    # tot.p1.resident.no.var[[b]] <- raw.results$results.p1.resident.no.var
    # tot.p1.resident.kappa.var[[b]] <- raw.results$results.p1.resident.kappa.var
    # tot.p1.resident.sigma.var[[b]] <- raw.results$results.p1.resident.sigma.var
    # tot.p2.resident[[b]] <- raw.results$results.p2.resident
    # tot.p2.resident.no.var[[b]] <- raw.results$results.p2.resident.no.var
    # tot.p2.resident.kappa.var[[b]] <- raw.results$results.p2.resident.kappa.var
    # tot.p2.resident.sigma.var[[b]] <- raw.results$results.p2.resident.sigma.var
    
    decomp.results$parameter.list$equil.times <- length(decomp.results$parameter.list$equil.times)
    tot.parameter.list[[b]] <- decomp.results$parameter.list
    tot.p1.invader.decomp.indiv[[b]] <- decomp.results$p1.invader.decomp.indiv
    tot.p1.invader.decomp.compare[[b]] <- decomp.results$p1.invader.decomp.compare
    tot.p2.invader.decomp.indiv[[b]] <- decomp.results$p2.invader.decomp.indiv
    tot.p2.invader.decomp.compare[[b]] <- decomp.results$p2.invader.decomp.compare
    
  }
  
}

temp1.p1 <- as.data.frame(tot.p1.invader.decomp.compare)
temp2.p1 <- data.frame(t(temp1.p1), row.names = NULL)
temp.mean.p1 <- unlist(lapply(temp2.p1, mean))
temp.se.p1 <- unlist(lapply(temp2.p1, sd))

ldgr.iden <- as.factor(c(1,0,0,0,0))
decomp <- c("LDGR", "e0", "eL", "eD", "eL*D")

time.step.3.decomp.p1 <- data.frame("decomp" = decomp, "mean" = temp.mean.p1, "se" = temp.se.p1, ldgr.iden)
time.step.3.decomp.p1$decomp <- factor(time.step.3.decomp.p1$decomp, levels = time.step.3.decomp.p1$decomp)


setwd("./../..")



#### TIME STEP 5 #####

tot.parameter.list <- 
  tot.p1.invader.decomp.indiv <- 
  tot.p1.invader.decomp.compare <- 
  tot.p2.invader.decomp.indiv <-
  tot.p2.invader.decomp.compare <- 
  list()

numSims <- 31
setwd("./results/time_step_5")


for (b in 1:numSims) {
  if (file.exists(paste("decomp_results_", b, ".RData", sep = "")) == TRUE) {
    # load(paste("raw_results_", b, ".RData", sep = ""))
    load(paste("decomp_results_", b, ".RData", sep = ""))
    
    # tot.p1.p2.resident.equilibrium[[b]] <- raw.results$results.p1.p2.equilibrium
    # tot.p1.resident[[b]] <- raw.results$results.p1.resident
    # tot.p1.resident.no.var[[b]] <- raw.results$results.p1.resident.no.var
    # tot.p1.resident.kappa.var[[b]] <- raw.results$results.p1.resident.kappa.var
    # tot.p1.resident.sigma.var[[b]] <- raw.results$results.p1.resident.sigma.var
    # tot.p2.resident[[b]] <- raw.results$results.p2.resident
    # tot.p2.resident.no.var[[b]] <- raw.results$results.p2.resident.no.var
    # tot.p2.resident.kappa.var[[b]] <- raw.results$results.p2.resident.kappa.var
    # tot.p2.resident.sigma.var[[b]] <- raw.results$results.p2.resident.sigma.var
    
    decomp.results$parameter.list$equil.times <- length(decomp.results$parameter.list$equil.times)
    tot.parameter.list[[b]] <- decomp.results$parameter.list
    tot.p1.invader.decomp.indiv[[b]] <- decomp.results$p1.invader.decomp.indiv
    tot.p1.invader.decomp.compare[[b]] <- decomp.results$p1.invader.decomp.compare
    tot.p2.invader.decomp.indiv[[b]] <- decomp.results$p2.invader.decomp.indiv
    tot.p2.invader.decomp.compare[[b]] <- decomp.results$p2.invader.decomp.compare
    
  }
  
}

temp1.p1 <- as.data.frame(tot.p1.invader.decomp.compare)
temp2.p1 <- na.omit(data.frame(t(temp1.p1), row.names = NULL))
temp.mean.p1 <- unlist(lapply(temp2.p1, mean))
temp.se.p1 <- unlist(lapply(temp2.p1, sd))

ldgr.iden <- as.factor(c(1,0,0,0,0))
decomp <- c("LDGR", "e0", "eL", "eD", "eL*D")

time.step.5.decomp.p1 <- data.frame("decomp" = decomp, "mean" = temp.mean.p1, "se" = temp.se.p1, ldgr.iden)
time.step.5.decomp.p1$decomp <- factor(time.step.5.decomp.p1$decomp, levels = time.step.5.decomp.p1$decomp)


setwd("./../..")

#### TIME STEP 7 #####

tot.parameter.list <- 
  tot.p1.invader.decomp.indiv <- 
  tot.p1.invader.decomp.compare <- 
  tot.p2.invader.decomp.indiv <-
  tot.p2.invader.decomp.compare <- 
  list()

numSims <- 31
setwd("./results/time_step_7")


for (b in 1:numSims) {
  if (file.exists(paste("decomp_results_", b, ".RData", sep = "")) == TRUE) {
    # load(paste("raw_results_", b, ".RData", sep = ""))
    load(paste("decomp_results_", b, ".RData", sep = ""))
    
    # tot.p1.p2.resident.equilibrium[[b]] <- raw.results$results.p1.p2.equilibrium
    # tot.p1.resident[[b]] <- raw.results$results.p1.resident
    # tot.p1.resident.no.var[[b]] <- raw.results$results.p1.resident.no.var
    # tot.p1.resident.kappa.var[[b]] <- raw.results$results.p1.resident.kappa.var
    # tot.p1.resident.sigma.var[[b]] <- raw.results$results.p1.resident.sigma.var
    # tot.p2.resident[[b]] <- raw.results$results.p2.resident
    # tot.p2.resident.no.var[[b]] <- raw.results$results.p2.resident.no.var
    # tot.p2.resident.kappa.var[[b]] <- raw.results$results.p2.resident.kappa.var
    # tot.p2.resident.sigma.var[[b]] <- raw.results$results.p2.resident.sigma.var
    
    decomp.results$parameter.list$equil.times <- length(decomp.results$parameter.list$equil.times)
    tot.parameter.list[[b]] <- decomp.results$parameter.list
    tot.p1.invader.decomp.indiv[[b]] <- decomp.results$p1.invader.decomp.indiv
    tot.p1.invader.decomp.compare[[b]] <- decomp.results$p1.invader.decomp.compare
    tot.p2.invader.decomp.indiv[[b]] <- decomp.results$p2.invader.decomp.indiv
    tot.p2.invader.decomp.compare[[b]] <- decomp.results$p2.invader.decomp.compare
    
  }
  
}

temp1.p1 <- as.data.frame(compact(tot.p1.invader.decomp.compare))
temp2.p1 <- data.frame(t(temp1.p1), row.names = NULL)
temp.mean.p1 <- unlist(lapply(temp2.p1, mean, na.rm = TRUE))
temp.se.p1 <- unlist(lapply(temp2.p1, sd))

ldgr.iden <- as.factor(c(1,0,0,0,0))
decomp <- c("LDGR", "e0", "eL", "eD", "eL*D")

time.step.7.decomp.p1 <- data.frame("decomp" = decomp, "mean" = temp.mean.p1, "se" = temp.se.p1, ldgr.iden)
time.step.7.decomp.p1$decomp <- factor(time.step.7.decomp.p1$decomp, levels = time.step.7.decomp.p1$decomp)


setwd("./../..")



#### TIME STEP 10 #####

tot.parameter.list <- 
  tot.p1.invader.decomp.indiv <- 
  tot.p1.invader.decomp.compare <- 
  tot.p2.invader.decomp.indiv <-
  tot.p2.invader.decomp.compare <- 
  list()

numSims <- 31
setwd("./results/time_step_10")


for (b in 1:numSims) {
  if (file.exists(paste("decomp_results_", b, ".RData", sep = "")) == TRUE) {
    # load(paste("raw_results_", b, ".RData", sep = ""))
    load(paste("decomp_results_", b, ".RData", sep = ""))
    
    # tot.p1.p2.resident.equilibrium[[b]] <- raw.results$results.p1.p2.equilibrium
    # tot.p1.resident[[b]] <- raw.results$results.p1.resident
    # tot.p1.resident.no.var[[b]] <- raw.results$results.p1.resident.no.var
    # tot.p1.resident.kappa.var[[b]] <- raw.results$results.p1.resident.kappa.var
    # tot.p1.resident.sigma.var[[b]] <- raw.results$results.p1.resident.sigma.var
    # tot.p2.resident[[b]] <- raw.results$results.p2.resident
    # tot.p2.resident.no.var[[b]] <- raw.results$results.p2.resident.no.var
    # tot.p2.resident.kappa.var[[b]] <- raw.results$results.p2.resident.kappa.var
    # tot.p2.resident.sigma.var[[b]] <- raw.results$results.p2.resident.sigma.var
    
    decomp.results$parameter.list$equil.times <- length(decomp.results$parameter.list$equil.times)
    tot.parameter.list[[b]] <- decomp.results$parameter.list
    tot.p1.invader.decomp.indiv[[b]] <- decomp.results$p1.invader.decomp.indiv
    tot.p1.invader.decomp.compare[[b]] <- decomp.results$p1.invader.decomp.compare
    tot.p2.invader.decomp.indiv[[b]] <- decomp.results$p2.invader.decomp.indiv
    tot.p2.invader.decomp.compare[[b]] <- decomp.results$p2.invader.decomp.compare
    
  }
  
}

temp1.p1 <- as.data.frame(compact(tot.p1.invader.decomp.compare))
temp2.p1 <- data.frame(t(temp1.p1), row.names = NULL)
temp.mean.p1 <- unlist(lapply(temp2.p1, mean, na.rm = TRUE))
temp.se.p1 <- unlist(lapply(temp2.p1, sd))

ldgr.iden <- as.factor(c(1,0,0,0,0))
decomp <- c("LDGR", "e0", "eL", "eD", "eL*D")


time.step.10.decomp.p1 <- data.frame("decomp" = decomp, "mean" = temp.mean.p1, "se" = temp.se.p1, ldgr.iden)
time.step.10.decomp.p1$decomp <- factor(time.step.10.decomp.p1$decomp, levels = time.step.10.decomp.p1$decomp)


setwd("./../..")

#####################

time.step.3 <- ggplot(data = time.step.3.decomp.p1, aes(x = decomp, y = mean, fill = ldgr.iden)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  scale_fill_manual(values = c("gray80", "gray40")) +
  # scale_y_continuous(limits = c(-8,15), breaks = c(-5, 0, 5, 10, 15)) +
  theme_classic()

time.step.5 <- ggplot(data = time.step.5.decomp.p1, aes(x = decomp, y = mean, fill = ldgr.iden)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  scale_fill_manual(values = c("gray80", "gray40")) +
  # scale_y_continuous(limits = c(-8,15), breaks = c(-5, 0, 5, 10, 15)) +
  theme_classic()



time.step.7 <- ggplot(data = time.step.7.decomp.p1, aes(x = decomp, y = mean, fill = ldgr.iden)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  scale_fill_manual(values = c("gray80", "gray40")) +
  # scale_y_continuous(limits = c(-8,15), breaks = c(-5, 0, 5, 10, 15)) +
  theme_classic()

time.step.10 <- ggplot(data = time.step.10.decomp.p1, aes(x = decomp, y = mean, fill = ldgr.iden)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  scale_fill_manual(values = c("gray80", "gray40")) +
  # scale_y_continuous(limits = c(-8,15), breaks = c(-5, 0, 5, 10, 15)) +
  theme_classic()


grid.arrange(time.step.3, time.step.5, 
             time.step.7, time.step.10, ncol = 2)

