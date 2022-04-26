################################################################################
#                                                                              #
# JHB-MPR-frailty-Code-Share                                                   #
#                                                                              #
################################################################################

# This R script file implements the h-likelihood procedure of 
# "Multi-Parameter Regression Survival Modelling with Random Effects" 

# The models are fitted to the Bladder Cancer dataset analysed in Section 4.3 
# of the paper.

# The necessary functions are defined first along with detailed explanations of 
# the purpose of each one of them. Then, these functions are applied to the said 
# dataset.

# The code runs cleanly sequentially (and requires installing the frailtyHL and 
# mpr pacakges). 

library(frailtyHL)
library(mpr)

################################################################################
#                                                                              #
# Functions                                                                    #
#                                                                              #
################################################################################
## * Function: weibhlike
## * Description: Defines the h-likelihood function.
## * Arguments: thetav = a vector of fixed effects followed by the random effects, 
##              surdata = a data.frame containing survival data where the first 
##              column is the survival time, the second is the censoring indicator, 
##              and the remaining columns correspond to covariates, with the last 
##              column corresponding to the centre,
##              k = number of scale (or shape) covariates plus one (for intercept), 
##              Disper = the frailty dispersion parameters
## * Note: This function is used as part of other functions. 

weibhlike <- function(thetav, surdata, k, q, disper, frailtystruc) {
  
  # if (!frailtystruc %in% c("BVNF", "ScF", "ShF" )){#(all(frailtystruc != c("BVNF", "ScF", "ShF" ))) {
  #   stop("Must specify the fariltystruc to be one of BVNF, ScF or ShF.")
  # } #  set it in the overall algorithm
  # 
  p <- length(thetav)
  
  tij <- surdata[, 1]
  deltai <- surdata[, 2]
  X <- as.matrix(surdata[, 3:(k + 2)])
  Z <- model.matrix(~(surdata[,ncol(surdata)])+0)
  
  beta <- thetav[1:k]
  alpha <- thetav[(k + 1):(k*2)]
  
  xbi <- X %*% beta
  xai <- X %*% alpha
  tau <- exp(xbi)
  gam <- exp(xai)
  
  if (frailtystruc == "BVNF") {
    
    sigbet <- disper[1]
    sigalp <- disper[2]
    rho <- disper[3]
    
    vbi <- thetav[((k*2) + 1):((k*2) + q)]
    vai <- thetav[((k*2) + q + 1):p]
    
    zvbet <- Z %*% vbi
    zvalp <- Z %*% vai
    
    eta.b <- tau*exp(zvbet)
    eta.a <- gam*exp(zvalp)
    
    hlike2 <- sum((-1/2)*((log(4*(pi^2)*(sigbet^2)*(sigalp^2)*(1 - (rho^2)))) + 
                            ((1/(1 - (rho^2)))*(((vbi^2)/(sigbet^2)) + ((vai^2)/(sigalp^2)) - 
                                                  ((2*rho*vbi*vai)/(sigbet*sigalp))))))
  } else if (frailtystruc == "ScF" | frailtystruc == "ShF") {
    sigma <- disper
    v <- thetav[-(1:(k*2))]
    zv <- Z %*% v
    
    if (frailtystruc == "ScF") {
      eta.b <- tau*exp(zv) #frailty in the scale parameter 
      eta.a <- gam
    } else if (frailtystruc == "ShF") {
      eta.b <- tau
      eta.a <- gam*exp(zv) #frailty in the shape parameter 
    }
    
    hlike2 <- sum((-1/2)*(log(2*pi*(sigma^2)) + ((v^2)/(sigma^2))))
  }
  
  hlike1 <- sum((deltai*(log(eta.b) + log(eta.a) + ((eta.a - 1)*log(tij)))) - 
                  (eta.b*(tij^eta.a)))
  
  hlike <- hlike1 + hlike2
  
  hlike
}
################################################################################
## * Function: Hmat
## * Description: Defines the H matrix (i.e. matrix of (-) second derivatives or
##                observed information matrix)
## * Arguments: thetav = a vector of fixed and random effects, 
##              surdata = a data.frame containing survival data where the first 
##              column is the survival time, the second is the censoring indicator, 
##              and the remaining columns correspond to covariates, with the last 
##              column corresponding to the centre,
##              k = number of scale (or shape) covariates plus one (for intercept), 
##              Disper = the frailty dispersion parameters
## * Note: This function is used as part of other functions. 

Hmat <- function(thetav, surdata, k, q, disper, frailtystruc){
  
  tij <- surdata[, 1]
  deltai <- surdata[, 2]
  X <- as.matrix(surdata[, 3:(k + 2)])
  Z <- model.matrix(~(surdata[,ncol(surdata)])+0)
  
  beta <- thetav[1:k]
  alpha <- thetav[(k + 1):(k*2)]
  
  xbi <- X %*% beta
  xai <- X %*% alpha
  tau <- exp(xbi)
  gam <- exp(xai)
  
  if (frailtystruc == "BVNF") {
    
    p <- length(thetav)
    sigbet <- disper[1]
    sigalp <- disper[2]
    rho <- disper[3]
    
    vbi <- thetav[((k*2) + 1):((k*2) + q)]
    vai <- thetav[((k*2) + q + 1):p]
    
    zvbet <- Z %*% vbi
    zvalp <- Z %*% vai
    
    eta.b <- tau*exp(zvbet)
    eta.a <- gam*exp(zvalp)
    
    # some terms that appear several times
    wb <- eta.b*(tij^eta.a)
    wba <- wb*eta.a*log(tij)
    wa <- (eta.a*log(tij))*((wb * ((eta.a*log(tij)) + 1)) - deltai)
    
    # (-) Second partial derivatives 
    Ib <- t(X * c(wb)) %*% X
    Iba <- t(X * c(wba)) %*% X
    Ibvb <- t(X * c(wb)) %*% Z
    Ibva <- t(X * c(wba)) %*% Z
    
    Iab <- t(Iba)
    Ia <- t(X * c(wa)) %*% X
    Iavb <- t(X * c(wba)) %*% Z
    Iava <- t(X * c(wa)) %*% Z
    
    Ivbb <- t(Ibvb)
    Ivba <- t(Iavb)
    Ivb <- (t(Z*c(wb)) %*% Z) + (diag(length(c(vbi)))/((1 - (rho^2))*(sigbet^2)))
    
    Ivbva <- (t(Z * c(wba)) %*% Z) - ((diag(length(c(vbi)))) * (rho/((1 - (rho^2))*sigbet*sigalp)))
    
    Ivab <- t(Ibva) 
    Ivaa <- t(Iava)
    Ivavb <- Ivbva
    Iva <- (t(Z * c(wa)) %*% Z) + (diag(length(c(vai)))/((1 - (rho^2))*(sigalp^2)))
    
    # the negative information matrix
    H <- rbind(cbind(Ib, Iba, Ibvb, Ibva),
               cbind(Iab, Ia, Iavb, Iava),
               cbind(Ivbb, Ivba, Ivb, Ivbva),
               cbind(Ivab, Ivaa, Ivavb, Iva))
  }
  
  if (frailtystruc == "ScF" | frailtystruc == "ShF") {
    sigma <- disper
    v <- thetav[-(1:(k*2))]
    zv <- Z %*% v
    if (frailtystruc == "ScF") {
      eta.b <- tau*exp(zv) #frailty in the scale parameter 
      eta.a <- gam
    } else if (frailtystruc == "ShF") {
      eta.b <- tau
      eta.a <- gam*exp(zv) #frailty in the shape parameter 
    }
    # some terms that appear several times
    wb <- eta.b*(tij^eta.a)
    wba <- wb*eta.a*log(tij)
    wa <- (eta.a*log(tij))*((wb * ((eta.a*log(tij)) + 1)) - deltai)
    
    # (-) Second partial derivatives 
    Ib <- t(X * c(wb)) %*% X
    Iba <- t(X * c(wba)) %*% X
    Iab <- t(Iba)
    Ia <- t(X * c(wa)) %*% X
    
    if (frailtystruc == "ScF") {
      Ibv <- t(X * c(wb)) %*% Z
      Iav <- t(X * c(wba)) %*% Z
      Iv <- (t(Z*c(wb)) %*% Z) + (diag(length(c(v)))/(sigma^2))
      
    } else if (frailtystruc == "ShF") {
      Ibv <- t(X * c(wba)) %*% Z 
      Iav <- t(X * c(wa)) %*% Z
      Iv <- (t(Z * c(wa)) %*% Z) + (diag(length(c(v)))/(sigma^2))
    }
    
    Ivb <- t(Ibv)
    Iva <- t(Iav)
    H <- rbind(cbind(Ib, Iba, Ibv),
               cbind(Iab, Ia, Iav),
               cbind(Ivb, Iva, Iv))
    
  }
  
  H
}

## The vector of first derivatives
Uvec <-  function(thetav, surdata, k, q, disper, frailtystruc){
  
  tij <- surdata[, 1]
  deltai <- surdata[, 2]
  X <- as.matrix(surdata[, 3:(k + 2)])
  Z <- model.matrix(~(surdata[,ncol(surdata)])+0)
  
  beta <- thetav[1:k]
  alpha <- thetav[(k + 1):(k*2)]
  
  xbi <- X %*% beta
  xai <- X %*% alpha
  tau <- exp(xbi)
  gam <- exp(xai)
  
  if (frailtystruc == "BVNF") {
    
    p <- length(thetav)
    sigbet <- disper[1]
    sigalp <- disper[2]
    rho <- disper[3]
    
    vbi <- thetav[((k*2) + 1):((k*2) + q)]
    vai <- thetav[((k*2) + q + 1):p]
    
    zvbet <- Z %*% vbi
    zvalp <- Z %*% vai
    
    eta.b <- tau*exp(zvbet)
    eta.a <- gam*exp(zvalp)
    
    ub <- deltai - (eta.b*(tij^eta.a))
    Ub <- t(X) %*% ub
    
    Uvb <- (t(Z) %*% ub) - ((1/(1 - (rho^2)))*((vbi/(sigbet^2)) - 
                                                 ((rho*vai)/(sigalp*sigbet))))
    
    ua <- (deltai*(1 + (eta.a*log(tij)))) - (eta.b*eta.a*(tij^eta.a)*log(tij))
    Ua <- t(X) %*% ua
    
    Uva <- (t(Z) %*% ua) - ((1/(1 - (rho^2)))*((vai/(sigalp^2)) - 
                                                 ((rho*vbi)/(sigalp*sigbet))))
    U <- c(Ub, Ua, Uvb, Uva)
    
  } else if  (frailtystruc == "ScF" | frailtystruc == "ShF") {
    
    sigma <- disper
    v <- thetav[-(1:(k*2))]
    zv <- Z %*% v
    
    if (frailtystruc == "ScF") {
      eta.b <- tau*exp(zv) #frailty in the scale parameter 
      eta.a <- gam
    } else if (frailtystruc == "ShF") {
      eta.b <- tau
      eta.a <- gam*exp(zv) #frailty in the shape parameter 
    }
    
    ub <- deltai - (eta.b*(tij^eta.a))
    Ub <- t(X) %*% ub
    ua <- (deltai*(1 + (eta.a*log(tij)))) - (eta.b*(tij^eta.a)*eta.a*log(tij)) 
    Ua <- t(X) %*% ua
    
    if (frailtystruc == "ScF") {
      Uv <- (t(Z) %*% ub) - (v/(sigma^2))
    } else if (frailtystruc == "ShF") {
      Uv <- (t(Z) %*% ua) - (v/(sigma^2))
    }
    U <- c(Ub, Ua, Uv)
  }
  U
}

################################################################################
## * Function: pbvh
## * Description: Defines the adjusted profile likelihood for
## * Arguments: tran.disper = transformed frailty dispersion parameters
##              thetav = a vector of fixed and random effects, 
##              surdata = a data.frame containing survival data where the first 
##              column is the survival time, the second is the censoring indicator, 
##              and the remaining columns correspond to covariates, with the last 
##              column corresponding to the centre,
##              k = number of scale (or shape) covariates plus one (for intercept),
##              q = the number of clusters or centres,
## * Note: This function is used in the estimation of the dispersion parameters  
##         (using nlm for the maximisation).

pbvh <-  function(tran.disper, surdata, thetav, k, q, frailtystruc){
  
  if (frailtystruc == "BVNF") {
    sigma <- exp(tran.disper[c(1,2)])
    rho1 <- tran.disper[3]
    rho <- (2/(1 + exp(-rho1))) - 1
    disper <- c(sigma, rho)
  } else if (frailtystruc == "ScF" | frailtystruc == "ShF") {
    disper <- exp(tran.disper)
  }
  
  h <- weibhlike(thetav = thetav, surdata = surdata, k = k, q = q, 
                 frailtystruc = frailtystruc, disper = disper)
  H <- Hmat(thetav = thetav, surdata = surdata, k = k, q = q, 
            frailtystruc = frailtystruc, disper = disper)
  
  logdetH <- determinant(H/(2*pi))$modulus[1]
  
  plik <- h - ((1/2)*logdetH)
  
  plik <- -plik
  
  plik
}
################################################################################
## * Function: hlikeNReq
## * Description: Fits the h-likelihood, i.e., for a given set of dispersion 
##                parameters this function estimates the fixed and random effects 
##                using the iterative estimation equations (Newton-Raphson) 
##                from the paper.
## * Arguments: surdata = a data.frame containing survival data where the first 
##              column is the survival time, the second is the censoring indicator, 
##              and the remaining columns correspond to covariates, with the last 
##              column corresponding to the centre,
##              init = a vector of initial values for the fixed and random effects, 
##              k = number of scale (or shape) covariates plus one (for intercept), 
##              q = the number of clusters or centres,
##              disper = frailty dispersion parameters,
##              tol = the tolerance defining when the optimization algorithm stops
##              (i.e., the largest change in the parameter vector is less than tol), iterlim = the maximum number of iterations of the
##              iterlim = the maximum number of iterations the optimization 
##              algorithm is allowed to take,
##              halfmax = the maximum number of step-halving steps carried out in 
##              the algorithm. 
## * Note: This function implements the first step of the fitting algorithm 
##         given in Section 2.4 of the paper.
## * Output: A list containing the estimates, the variance covariance matrix.

# hlikeNReq <- function(surdata, thetav.init, k, q, disper, frailtystruc, tol, 
#                       maxiter, halfmax){
#   
#   thetav <- thetav.init
#   sigbet <- disper[1]
#   sigalp <- disper[2]
#   rho <- disper[3]
#   
#   p <- length(thetav)
#   
#   tij <- surdata[, 1]
#   deltai <- surdata[, 2]
#   X <- as.matrix(surdata[, 3:(k + 2)])
#   Z <- model.matrix(~(surdata[,ncol(surdata)])+0)
#   
#   iter <- 0
#   thetavold <- Inf
#   
#   while (max(abs(thetav - thetavold)) > tol & iter < maxiter) {
#     
#     thetavold <- thetav
#     beta <- thetav[1:k]
#     alpha <- thetav[(k + 1):(k*2)]
#     vbi <- thetav[((k*2) + 1):((k*2) + q)]
#     vai <- thetav[((k*2) + q + 1):p]
#     
#     zvbet <- Z %*% vbi
#     zvalp <- Z %*% vai
#     
#     xbi <- X %*% beta
#     xai <- X %*% alpha
#     tau <- exp(xbi)
#     gam <- exp(xai)
#     
#     eta.b <- tau*exp(zvbet)
#     eta.a <- gam*exp(zvalp)
#     
#     # first derivatives
#     ub <- deltai - (eta.b*(tij^eta.a))
#     Ub <- t(X) %*% ub
#     
#     Uvb <- (t(Z) %*% ub) - ((1/(1 - (rho^2)))*((vbi/(sigbet^2)) - ((rho*vai)/(sigalp*sigbet))))
#     
#     ua <- (deltai*(1 + (eta.a*log(tij)))) - (eta.b*eta.a*(tij^eta.a)*log(tij))
#     Ua <- t(X) %*% ua
#     
#     Uva <- (t(Z) %*% ua) - ((1/(1 - (rho^2)))*((vai/(sigalp^2)) - ((rho*vbi)/(sigalp*sigbet))))
#     
#     dhdthetav <- c(Ub, Ua, Uvb, Uva)
#     
#     # (-) second derivatives
#     H <- Hmat(surdata = surdata, k = k, q = q, thetav = thetav, disper = disper,
#               frailtystruc = frailtystruc)
#     
#     # solving for the difference between parameters
#     pardiff <- solve(H,dhdthetav)
#     
#     hlikeold <- weibhlike(thetav = thetav, surdata = surdata, k = k, q = q, 
#                           disper = disper, frailtystruc = frailtystruc)
#     
#     hlike <- -Inf
#     
#     j <- 0
#     del <- 1
#     
#     # step halving when Netwon-Raphson step is too large (keep halving step 
#     # until the h-likelihood is increased).
#     while (hlike < hlikeold & j < halfmax) {
#       
#       del <- del/(2^j)
#       
#       thetav <- as.vector(thetavold + del*pardiff)
#       
#       hlike <- weibhlike(thetav = thetav, surdata = surdata, k = k, q = q, 
#                          disper = disper, frailtystruc = frailtystruc)
#       
#       hlike <- ifelse(is.na(hlike), -Inf, hlike)
#       j <- j + 1
#     }
#     iter <- iter + 1
#   }
#   # store quantities of interest
#   vcovmat <- solve(Hmat(surdata = surdata, k = k, thetav = thetav, q = q,
#                         disper = disper, frailtystruc = frailtystruc))
#   #variance covariance matrix can be estimated later too
#   
#   list(thetav = thetav, vcovmat = vcovmat)
# }

## Newton Raphson to fit the h-likelihood
hlikeNReq <- function(surdata, thetav.init, k, disper, maxiter, tol, halfmax, q, 
                      frailtystruc){
  thetav <- thetav.init
  thetaold <- Inf
  iter <- 0
  
  while (max(abs(thetav - thetaold)) > tol & iter < maxiter) {
    
    thetaold <- thetav
    # first derivatives
    U <- Uvec(surdata = surdata, k = k, q = q, thetav = thetav, disper = disper, 
              frailtystruc = frailtystruc)
    # (-) second derivatives
    H <- Hmat(surdata = surdata, k = k, q = q, thetav = thetav, disper = disper,
              frailtystruc = frailtystruc)
    # solving for the difference between parameters
    pardiff <- solve(H,U)
    
    hlikeold <- weibhlike(thetav = thetav, surdata = surdata, k = k, q = q, 
                          disper = disper, frailtystruc = frailtystruc)
    hlike <- -Inf
    
    j <- 0
    del <- 1
    
    # step halving when Netwon-Raphson step is too large
    while (hlike < hlikeold & j < halfmax) {
      
      del <- del/(2^j)
      
      thetav <- as.vector(thetaold + del*pardiff)
      
      hlike <- weibhlike(thetav = thetav, surdata = surdata, k = k, q = q, 
                         disper = disper, frailtystruc = frailtystruc)
      
      hlike <- ifelse(is.na(hlike), -Inf, hlike)
      j <- j + 1
    }
    iter <- iter + 1
  }
  # store quantities of interest
  vcovmat <- solve(Hmat(surdata = surdata, k = k, thetav = thetav, q = q,
                        disper = disper, frailtystruc = frailtystruc))
  #variance covariance matrix can be estimated later too
  list(thetav = thetav, vcovmat = vcovmat)
}

################################################################################
## * Function: HLfit.algo
## * Description: Gives the overall model fitting algorithm, summarised in Section
##                2.4.
## * Arguments: surdata = a data.frame containing survival data where the first 
##              column is the survival time, the second is the censoring indicator, 
##              and the remaining columns correspond to covariates, with the last 
##              column corresponding to the centre,
##              init = a vector of initial values for the fixed and random effects, 
##              k = number of scale (or shape) covariates plus one (for intercept), 
##              disper0 = the initial values for frailty dispersion parameter(s),
##              tol = the tolerance defining when the optimization algorithm stops
##              (i.e., the largest change in the parameter vector is less than tol), iterlim = the maximum number of iterations of the
##              iterlim = the maximum number of iterations the optimization 
##              algorithm is allowed to take,
##              halfmax = the maximum number of step-halving steps carried out in 
##              the algorithm. 
## * Note: This function implements the overall h-likelihood fitting algorithm 
##         summarised in Section 2.4 of the paper.
## * Output: A list containing estimates of the different sets of parameters, 
##           their standard errors and the number of iterations the algorithm took.

HLfit.algo <- function(surdata, thetav.init, k, q, disper.init, frailtystruc,
                       tol, maxiter, halfmax){
  # extract initial values
  lam <- c(thetav.init, disper.init)
  lamold <- Inf
  iter = 0
  
  while (max(abs(lamold - lam)) > tol & iter <= maxiter) {
    
    lamold <- lam
    
    if (frailtystruc == "BVNF") {
      thetav <- lam[1:((k*2) + (q*2))]
      disper <- lam[((k*2) + (q*2) + 1):(length(lam))]
    } else if (frailtystruc == "ScF" | frailtystruc == "ShF") {
      thetav <- lam[1:((k*2) + q)]
      disper <- lam[-(1:((k*2) + q))]
    }
    
    # STEP 1: maximise h to obtain estimates of the fixed and random effects 
    hlike.fit <- hlikeNReq(surdata = surdata, thetav.init = thetav, 
                           k = k, q = q, disper = disper, frailtystruc = frailtystruc,
                           maxiter = maxiter, tol = tol, halfmax = halfmax)
    
    # update thetav
    thetav <- hlike.fit$thetav
    # transform the dispersion parameters to ensure the algorithm remains in the 
    # desired parameter space.
    if (frailtystruc == "BVNF") {
      tran.rho <- log((1 + disper[3])/(1 - disper[3]))
      
      tran.disper <- c(log(disper[c(1,2)]), tran.rho)
    } else if (frailtystruc == "ScF" | frailtystruc == "ShF") {
      tran.disper <- log(disper)
    }
    # STEP 2: maximise the adjusted profile 
    pbvh.nlmfit <- nlm(pbvh, tran.disper, thetav = thetav, frailtystruc = frailtystruc,
                       surdata = surdata, k = k, q = q,iterlim = maxiter, 
                       hessian = TRUE)
    # bring dispersion parameters back to original scale
    if (frailtystruc == "BVNF") {
      sigma <- exp(pbvh.nlmfit$est[c(1,2)])
      
      tran.rho <- pbvh.nlmfit$est[3]
      
      rho <- (2/(1 + exp(-tran.rho))) - 1
      
      lam <- c(thetav, sigma, rho)
    } else if (frailtystruc == "ScF" | frailtystruc == "ShF") {
      sigma <- exp(pbvh.nlmfit$est)
      lam <- c(thetav, sigma)
    }
    iter <- iter + 1
    # print(c(sigma, rho))
  }
  # calculate the standard errors 
  SE.thetav <- sqrt((diag(hlike.fit$vcovmat)))
  SE.disper <- sqrt(diag(solve(pbvh.nlmfit$hess)))
  
  # storing quantities of interest in a list 
  fixed.effects <- cbind(thetav[1:k], SE.thetav[1:k], thetav[(k + 1):(k*2)], 
                         SE.thetav[(k + 1):(k*2)])
  colnames(fixed.effects) <- c("coef.b", "se.b","coef.a", "se.a")
  
  ######## different for differ models
  if (frailtystruc == "BVNF") {
    # using the delta method to get the standard errors of the dispersion parameters
    # in their original scale
    
    SE.sig.beta <-  SE.disper[1] * sigma[1]
    SE.sig.alpha <- SE.disper[2] * sigma[2]
    SE.rho <- SE.disper[3] * ((2*exp(-tran.rho))/((1 + exp(-tran.rho))^2))
    
    random.effects <- cbind(thetav[((k*2) + 1):((k*2) + q)], 
                            SE.thetav[((k*2) + 1):((k*2) + q)],
                            thetav[((k*2) + q + 1):length(thetav)],
                            SE.thetav[((k*2) + q + 1):length(thetav)])
    colnames(random.effects) <- c("vbi", "se.vbi", "vai", "se.vai")
    
    dispersion.params <- cbind(c(sigma, rho), c(SE.sig.beta, SE.sig.alpha, SE.rho))
    rownames(dispersion.params) <- c("sigma.b", "sigma.a", "rho")
    
  } else if (frailtystruc == "ScF" | frailtystruc == "ShF") {
    SE.sig <-  SE.disper * sigma
    
    random.effects <- cbind(thetav[-(1:(k*2))], SE.thetav[-(1:(k*2))])
    dispersion.params <- cbind(sigma, SE.sig)
    
    if (frailtystruc == "ScF") {
      colnames(random.effects) <- c("vbi", "se.vbi")
      rownames(dispersion.params) <- c("sigma.b")
    } else if (frailtystruc == "ShF") {
      colnames(random.effects) <- c("vai", "se.vai")
      rownames(dispersion.params) <- c("sigma.a")
      
    }
    
  }
  
  colnames(dispersion.params) <- c("estimate", "se")
  rownames(random.effects) <- paste0("centre ", levels(surdata[, ncol(surdata)]))
  list(fixed.effects = fixed.effects, dispersion.params = dispersion.params,
       random.effects = random.effects, iter = iter)
}