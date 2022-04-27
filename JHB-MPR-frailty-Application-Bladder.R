##############################################################################
#                                                                            #
# JHB-MPR-frailty-Code-Share                                                 #
#                                                                            #
##############################################################################

# This R script file implements the h-likelihood procedure of 
# "Multi-Parameter Regression Survival Modelling with Random Effects" 
# by Fatima-Zahra Jaouimaa, Il Do Ha, and Kevin Burke
#
# The necessary functions are sourced from the github repository 
# "https://github.com/fatimZJ/JHB-MPR-Frailty", wherein a description of the
# purpose of each function can be found.
#
# The code below fits the BVNF model (bivariate normal frailties), the
# ScF model (scale frailty) and ShF model (shape frailty) to the bladder
# cancer dataset as per Section 4.3 of the paper.
#
# The code runs sequentially and requires the frailtyHL and mpr pacakges.
# An internet connection is also needed for the purpose of loading in the
# underlying R script from github.


##############################################################################
#                                                                            #
# Application to the Bladder Cancer dataset                                  #
#                                                                            #
##############################################################################


################ packages and github script

## clear current environment
rm(list = ls(all = TRUE))

## loading required packages
library(frailtyHL)
library(mpr)

## sources the functions from github
url <- "https://raw.githubusercontent.com/fatimZJ/JHB-MPR-Frailty/main/JHB-MPR-frailty-Code-Share-BVN-ScF-ShF.R"
script <- paste(readLines(url), collapse = "\n")
eval(parse(text = script), envir= .GlobalEnv)


################ loading and manipulating data

## load data
data(bladder0, package = "frailtyHL")

## omit observations with zero survival time
bladder <- bladder0[(bladder0$Surtime != 0), ]

## making factor variables
bladder[, c("Center", "Chemo", "Tustat")] <-
  lapply(bladder[, c("Center", "Chemo", "Tustat")], factor)

## changing the survival time from days to year
bladder$Surtime <- bladder$Surtime/365

## display data sctructure
str(bladder)

################ initialization

## model matrix
X <- model.matrix(~ Chemo + Tustat, data = bladder) 

## data is a data.frame used in the fitting function and whose columns
## must have the following order: the first is the survival time, the second
## is the censoring indicator, the next group are the covariates, and the
## last column is the centre ID.
data <- data.frame(bladder$Surtime, bladder$Status,
                      X, as.factor(bladder$Center))

## some constants that are needed
k <- dim(X)[2]
q <- nlevels(bladder$Center)

################ setting initial values

## this portion is used to get some "good" initial values that are used later
## when fitting the frailty models

## exponential dist parameter estimate 
exp_rate_mle <- sum(bladder$Status) / sum(bladder$Surtime)

init0 <- c(log(exp_rate_mle), 0.01)

# first fit a no frailty, no covariate Weibull model
weibull <- mpr(Surv(Surtime, Status) ~ 1, data = bladder, init = init0)

init0 <- c(weibull$coefficients$beta, rep(0.01, 2), 
           weibull$coefficients$alpha, rep(0.01, 2))

# then fit a Weibull MPR model (but still there is no frailty)
weibull_mpr <- mpr(Surv(Surtime, Status)~ Chemo + Tustat,
                   data = bladder, init = init0)


################ model fitting

## BVNF model: Weibull MPR with bivariate normal frailties

# initial values for fixed effects and frailties
thetav_init <- c(weibull_mpr$coefficients$beta,
                 weibull_mpr$coefficients$alpha, 
                 rep(0.01, (q*2)))

# initial values for frailty dispersion parameters
disper_init <- c(0.1, 0.1, -0.1)

BVNF <- NA

BVNF <- mpr_frailty(data = data, frailty = "BVNF", k = k, q = q,
                    thetav_init = thetav_init, disper_init = disper_init)

## note that ".b" indicates the scale (beta) coefficients and ".a" indicates
## the shape (alpha) coefficients

round(BVNF[[1]],2)
round(BVNF[[2]],2)
round(BVNF[[3]],2)


## ScF model: Weibull MPR with scale frailties

thetav_init <- c(weibull_mpr$coefficients$beta,
                 weibull_mpr$coefficients$alpha, 
                 rep(0.01, q))

disper_init <- 0.1

ScF <- NA

ScF <- mpr_frailty(data = data, frailty = "ScF", k = k, q = q, 
                  thetav_init = thetav_init, disper_init = disper_init)

round(ScF[[1]],2)
round(ScF[[2]],2)
round(ScF[[3]],2)


## ShF model: Weibull MPR with shape frailties

ShF <- NA

ShF <- mpr_frailty(data = data, frailty = "ShF", k = k, q = q, 
                   thetav_init = thetav_init, disper_init = disper_init)

round(ShF[[1]],2)
round(ShF[[2]],2)
round(ShF[[3]],2)


## END


