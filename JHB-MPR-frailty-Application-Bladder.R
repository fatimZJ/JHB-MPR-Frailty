################################################################################
#                                                                              #
# Application to the Bladder Cancer dataset                                    #
#                                                                              #
################################################################################
## loading required packages
library(devtools)
library(frailtyHL)
library(mpr)

## sources the functions from github
SourceURL <- "https://raw.githubusercontent.com/fatimZJ/JHB-MPR-Frailty/main/JHB-MPR-frailty-Code-Share-BVN-ScF-ShF.R"
source_url(SourceURL)

## loading the dataset
data(bladder0, package = "frailtyHL")

str(bladder0)
## some data manipulation 
bladder <- bladder0[(bladder0$Surtime != 0), ]

bladder[,c("Center", "Chemo", "Tustat")] <- lapply(bladder[,c("Center", "Chemo",
                                                              "Tustat")], factor)
str(bladder)
# checking for NAs
apply(bladder, 2, function(x) sum(is.na(x)))

# changing the survival time from days to year
bladder$Surtime <- bladder$Surtime/365

# model matrix
X <- model.matrix(~Chemo + Tustat, data = bladder) 

# defining the data.frame that gets passed into the model fitting function- first
# column is the survival time, the second is the censoring indicator, followed by 
# the covariates, and the last column corresponds to the centre,

surdata <- data.frame(bladder$Surtime, bladder$Status, X, as.factor(bladder$Center))

# defining other parameters 
k <- dim(X)[2]
q <- nlevels(bladder$Center)
maxiter <- 5000
halfmax <- 500
tol <- 1e-4

## generating starting values 
ti <- surdata[,1]
deltai <- surdata[,2]
lam.hat <- sum(deltai)/sum(ti)

bet0 <- log(lam.hat)
alp0 <- 0.01
init0 <- c(bet0, alp0)

# fitting standard no frailty, no covariate Weibull MPR model (to get starting values later)
nocov.weibfit <- mpr(Surv(Surtime,Status)~list(~1,~1), data = bladder, init = init0)

init <- c(nocov.weibfit$coefficients$beta, rep(0.01, (k - 1)), 
          nocov.weibfit$coefficients$alpha, rep(0.01, (k - 1)))

# fitting a no frailty weibull MPR model with covariates

weib.fit <- mpr(Surv(Surtime, Status)~list(~Chemo + Tustat,~Chemo + Tustat), data = bladder, 
                init = init)

thetav.init <- c(weib.fit$coefficients$beta, weib.fit$coefficients$alpha, 
                 rep(0.01, (q*2)))
disper.init <- c(0.1, 0.1, -0.1)

## fitting the Weibull MPR model assuming BVN random effects- BVNF 

estimates.BVNF <- NA

estimates.BVNF <- HLfit.algo(surdata = surdata, k = k, q = q, thetav.init = thetav.init, 
                        tol = tol, disper.init = disper.init, maxiter = maxiter, 
                        halfmax = halfmax, frailtystruc = "BVNF")

round(estimates.BVNF[[1]],2)
round(estimates.BVNF[[2]],2)
round(estimates.BVNF[[3]],2)

## fitting the 
thetav.init <- c(weib.fit$coefficients$beta, weib.fit$coefficients$alpha, 
                 rep(0.01, q))
disper.init <- 0.1

##Fitting the Weibull MPR model assuming scale random effects model- ScF
estimates.ScF <- NA

estimates.ScF <- HLfit.algo(surdata = surdata, k = k, q = q, thetav.init = thetav.init, tol = tol,
                   disper.init = disper.init, maxiter = maxiter, halfmax = halfmax,
                   frailtystruc = "ScF")

round(estimates.ScF[[1]],2)
round(estimates.ScF[[2]],2)
round(estimates.ScF[[3]],2)

##Fitting the Weibull MPR model assuming scale random effects model- ScF
estimates.ShF <- NA

estimates.ShF <- HLfit.algo(surdata = surdata, k = k, q = q, thetav.init = thetav.init, tol = tol,
                   disper.init = disper.init, maxiter = maxiter, halfmax = halfmax,
                   frailtystruc = "ShF")

round(estimates.ShF[[1]],2)
round(estimates.ShF[[2]],2)
round(estimates.ShF[[3]],2)


