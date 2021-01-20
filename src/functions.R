###################
## Model fitting ##
###################

library(expm)
library(mnormt)
library(ramcmc)

## Make state and observation matrices
make.matrices = function(ltheta) {
  
  ## Transform parameters
  theta = exp(ltheta)
  gamma   = theta[1]
  C       = theta[2:4]
  kappa   = theta[5:7]
  epsilon = theta[8]
  sigmaf  = theta[9]
  sigmat  = theta[10]
  fco2    = theta[11]
  
  if (length(ltheta) > 11) {
    nn = 5
  } else {
    nn = 4
  }
  
  ## State vector dimension
  if (length(ltheta) > 11) {
    nn = 5
    pp = 3
  } else {
    nn = 4
    pp = 1
  }
  
  ## Forecast matrix
  F = matrix(0, 2, nn)
  F[1,2] =   1
  F[2,1] =   1
  F[2,2] = - kappa[1]
  F[2,3] =   (1 - epsilon)*kappa[3]
  F[2,4] = - (1 - epsilon)*kappa[3]
  
  ## Transition matrix
  G = matrix(0, nn, nn)
  G[1,1] = - gamma
  G[2,1] =   1/C[1]
  G[2,2] = - (kappa[1] + kappa[2])/C[1]
  G[2,3] =   kappa[2]/C[1]
  G[3,2] =   kappa[2]/C[2]
  G[3,3] = - (kappa[2] + epsilon*kappa[3])/C[2]
  G[3,4] =   epsilon*kappa[3]/C[2]
  G[4,3] =   kappa[3]/C[3]
  G[4,4] = - kappa[3]/C[3]
  
  ## Transition uncertainty
  W = matrix(0, nn, nn)
  W[1,1] = sigmaf^2
  W[2,2] = sigmat^2/C[1]^2
  
  ## Transition intercept
  B = matrix(0, nn, pp)
  B[1,1] = fco2 * gamma
  
  if (nn > 4) {
    fvol   = theta[12]
    sigmad = theta[13]
    F[2,5] = 1
    G[2,5] = 1/C[1]
    G[5,5] = 1e-12
    W[5,5] = sigmad^2
    B[1,2] = fvol * gamma
    B[5,3] = 1
  }
  
  ## Discretization
  Gd = expm(G)
  Bd = solve(G, (Gd - diag(nn)) %*% B)
  H  = rbind(cbind(-G,W),cbind(matrix(0,nn,nn),t(G)))
  Hd = expm(H)
  Wd = t(Hd[(nn + 1):(2*nn), (nn + 1):(2*nn)]) %*% Hd[1:nn, (nn + 1):(2*nn)]
  
  ## Initial conditions
  
  return(list(F = F, G = Gd, W = Wd, B = Bd))
  
}

## Fit model using Kalman filter
fit.model = function(theta, y, X, V = NULL, nu = NULL, m0 = NULL, C0 = NULL) {
  
  if (is.null(nu)) {
    XX = X
  } else{
    XX = rbind(X, nu)
  }
  
  ## Make matrices
  matrices = make.matrices(theta)
  F = matrices$F
  G = matrices$G
  W = matrices$W
  B = matrices$B
  
  ## Extract dimensions
  rr = nrow(y)
  tt = ncol(y)
  nn = ncol(F)
  pp = nrow(XX)
  
  ## Intercept terms
  c = matrix(0, rr, tt)
  d = B %*% XX
  
  if (rr < nrow(F))
    F = F[1,,drop=FALSE]
  
  ## Forecast uncertainty
  if (is.null(V))
    V = diag(rep(1e-12, rr))
  
  ## Initial conditions
  if (is.null(m0))
    m0 = rep(0, nrow(G))
  if (is.null(C0)) {
    C = solve(diag(4^2) - 
                kronecker(G[1:4,1:4],G[1:4,1:4]), as.vector(W[1:4,1:4]))
    C = matrix(C, 4, 4)
    if (nn > 4) {
      C0 = matrix(0, 5, 5)
      C0[1:4,1:4] = C
      C0[5,5] = 1e-12
    } else {
      C0 = C
    }
  }
  
  ## Transform matrices for use in fkf
  F  = array(F, c(rr,nn,tt))
  G  = array(G, c(nn,nn,tt))
  V  = array(V, c(rr,rr,tt))
  W  = array(W, c(nn,nn,tt))
  
  ## Fit model
  filter.dlm(y, F, G, V, W, m0, C0, c, d)
  
}


###############################
## Latent parameters - theta ##
###############################

## Evaluate log prior
theta.prior = function (theta, phi, kappa = 1) {
  
  dmnorm(theta, phi$mu, phi$Sigma*kappa^2, TRUE)

}

## Evaluate log jacobian
theta.jacobian = function (theta, phi, kappa = 1) {
  
  0
  
}

## Evaluate log-likelihood
theta.likelihood = function(theta, y0, X0, y1, X1, nu) {
  
  ll = fit.model(theta[1:11], y0, X0, m0 = c(2*exp(theta[11]),0,0,0))$logLik +
    fit.model(theta, y1, X1, nu = nu)$logLik
  
  return(ifelse(is.finite(ll), ll, -Inf))
  
}

## Sample model parameters
sample.theta = function(theta, y0, X0, y1, X1, nu, phi, 
                        Sigma, adaptive = TRUE, jj = 1) {
  
  ll.theta  = theta.likelihood(theta, y0, X0, y1, X1, nu)
  lp.theta  = theta.prior     (theta, phi)
  lj.theta  = theta.jacobian  (theta, phi)
  
  u         = rnorm(length(theta))
  theta.    = as.numeric(theta + Sigma %*% u)
  names(theta.) = names(theta)
  
  ll.theta. = theta.likelihood(theta., y0, X0, y1, X1, nu)
  lp.theta. = theta.prior     (theta., phi)
  lj.theta. = theta.jacobian  (theta., phi)
  
  a = min(1, exp(ll.theta. + lp.theta. + lj.theta. - 
                   ll.theta - lp.theta - lj.theta))
  if (!is.finite(a))
    a = 0
  
  if (runif(1) < a) 
    theta = theta.
  
  if (adaptive)
    Sigma = adapt_S(Sigma, u, a, jj)
  
  list(theta = theta, acceptance = a, Sigma = Sigma)
  
}

## Sample all models
sample.thetas = function(thetas, y0s, X0, y1s, X1, nu, phi, 
                         Sigmas, adaptive = TRUE, jj = 1) {
  
  ## Initialise storage
  n.models    = length(thetas)
  acceptances = numeric(n.models)
  
  for (m in 1:n.models) {
    
    buffer = sample.theta(thetas[[m]], y0s[[m]], X0, y1s[[m]], X1, nu, phi, 
                          Sigmas[[m]], adaptive, jj)
    thetas[[m]]    = buffer$theta
    Sigmas[[m]]    = buffer$Sigma
    acceptances[m] = buffer$acceptance
    
  } ## m
  
  list(thetas = thetas, Sigmas = Sigmas, acceptances = acceptances)
  
}

unlist.theta = function(theta) {
  
  unlist(theta)
  
}


######################
## Parameters - phi ##
######################

sample.phi = function (phi, thetas, hyper.parameters) {
  
  n.models = length(thetas)
  
  theta  = simplify2array(thetas)
  Tau    = solve(phi$Sigma)
  mean   = hyper.parameters$mup %*% hyper.parameters$mum + 
    Tau %*% rowSums(theta)
  Cov    = solve(hyper.parameters$mup + n.models*Tau)
  mu     = rmnorm(1, Cov %*% mean, Cov)
  df     = hyper.parameters$Sigman + n.models
  Scale  = solve(hyper.parameters$Sigmas + tcrossprod(theta - mu))
  Sigma  = solve(rWishart(1, df, Scale)[,,1])
  Sigma  = 0.5*(Sigma + t(Sigma))
  
  list(mu = mu, Sigma = Sigma)
  
}

unlist.phi = function(phi) {
  
  c(phi$mu,sqrt(diag(phi$Sigma)),cov2cor(phi$Sigma)[lower.tri(phi$Sigma)])
  
}


#########################
## Observations - zeta ##
#########################

## Evaluate log prior
zeta.prior = theta.prior

## Evaluate log jacobian
zeta.jacobian = theta.jacobian

## Evaluate log-likelihood
zeta.likelihood = function(zeta, z, V, X, nu) {
  
  mask = 1:ncol(z)
  ll = fit.model(zeta, z, X[,mask,drop=FALSE], V = V, nu = nu[mask])$logLik
  
  return(ifelse(is.finite(ll), ll, -Inf))
  
}

## Sample observation parameters
sample.zeta = function(zeta, z, V, X, nu, phi, kappa,
                       Sigma, adaptive = TRUE, jj = 1) {
  
  ll.zeta  = zeta.likelihood(zeta, z, V, X, nu)
  lp.zeta  = zeta.prior     (zeta, phi, kappa)
  lj.zeta  = zeta.jacobian  (zeta, phi, kappa)
  
  u         = rnorm(length(zeta))
  zeta.    = as.numeric(zeta + Sigma %*% u)
  names(zeta.) = names(zeta)
  
  ll.zeta. = zeta.likelihood(zeta., z, V, X, nu)
  lp.zeta. = zeta.prior     (zeta., phi, kappa)
  lj.zeta. = zeta.jacobian  (zeta., phi, kappa)
  
  a = min(1, exp(ll.zeta. + lp.zeta. + lj.zeta. - ll.zeta - lp.zeta - lj.zeta))
  if (!is.finite(a))
    a = 0
  
  if (runif(1) < a) 
    zeta = zeta.
  
  if (adaptive)
    Sigma = adapt_S(Sigma, u, a, jj)
  
  list(zeta = zeta, acceptance = a, Sigma = Sigma)
  
}

unlist.zeta = function(zeta) {
  
  unlist(zeta)
  
}


###########################
## delta - Discrepancies ##
###########################

sample.delta = function(delta, theta, nu, y, X) {
  
  model = fit.model(theta, y, X, nu = nu)
  
  delta. = try(sample0.dlm(model)[5,], silent = TRUE)
  if (is.numeric(delta.)) {
    delta.
  } else {
    delta
  }
  
}

sample.deltas = function(deltas, thetas, nu, ys, X) {
  
  ## Initialise storage
  n.models = length(thetas)
  
  ## Loop over models
  for (m in 1:n.models)
    deltas[m,] = sample.delta(deltas[m,], thetas[[m]], nu, ys[[m]], X)

  ## Return samples
  return(deltas)
  
}


#############################
## nu - Shared discrepancy ##
#############################

sample.nu = function(nu, deltas, lthetas, sigman) {
  
  sigmad = exp(sapply(lthetas, function(x) x[13]))
  sigman = exp(sigman)
  
  ## Extract dimensions
  rr = nrow(deltas)
  tt = ncol(deltas)
  nn = rr + 1
  
  ## Make matrices
  F = matrix(0, rr, nn)
  diag(F[1:rr,1:rr]) = 1
  G = diag(nn)
  G[nn,nn] = 0
  V = diag(rep(1e-12, rr))
  W = matrix(sigman^2, nn, nn)
  diag(W[1:rr,1:rr]) = sigman^2 + sigmad^2
  
  ## Intercept terms
  c = matrix(0, rr, tt)
  d = matrix(0, nn, tt)
  
  ## Initial conditions
  m0 = numeric(nn)
  C0 = matrix(0, nn, nn)
  
  ## Transform matrices for use in fkf
  F  = array(F, c(rr,nn,tt))
  G  = array(G, c(nn,nn,tt))
  V  = array(V, c(rr,rr,tt))
  W  = array(W, c(nn,nn,tt))
  
  ## Fit model
  model = filter.dlm(deltas, F, G, V, W, m0, C0, c, d)
  
  ## Sample multi-model mean
  nu. = try(sample0.dlm(model)[nn,], silent = TRUE)
  
  if (is.numeric(nu.)) {
    nu.
  } else {
    nu
  }
  
}

#####################################
## sigma_n - Shared discrepancy sd ##
#####################################

## Evaluate log prior
sigman.prior = function (sigman, hyper.parameters) {
  
  dnorm(sigman, hyper.parameters$sigmanm, 1/sqrt(hyper.parameters$sigmanp), TRUE)
  
}

## Evaluate log jacobian
sigman.jacobian = function (sigman, hyper.parameters) {
  
  0
  
}

## Evaluate log-likelihood
sigman.likelihood = function(lsigman, nu) {
  
  ll = sum(dnorm(nu, 0, exp(lsigman), TRUE))
  
  return(ifelse(is.finite(ll), ll, -Inf))
  
}

## Sample multi-model mean variance
sample.sigman = function(sigman, nu, hyper.parameters, 
                         Sigma, adaptive = TRUE, jj = 1) {
  
  ll.sigman  = sigman.likelihood(sigman, nu)
  lp.sigman  = sigman.prior     (sigman, hyper.parameters)
  lj.sigman  = sigman.jacobian  (sigman, hyper.parameters)
  
  u          = rnorm(length(sigman))
  sigman.    = as.numeric(sigman + Sigma %*% u)
  names(sigman.) = names(sigman)
  
  ll.sigman. = sigman.likelihood(sigman., nu)
  lp.sigman. = sigman.prior     (sigman., hyper.parameters)
  lj.sigman. = sigman.jacobian  (sigman., hyper.parameters)
  
  a = min(1, exp(ll.sigman. + lp.sigman. + lj.sigman. 
                 - ll.sigman - lp.sigman - lj.sigman))
  if (!is.finite(a))
    a = 0
  
  if (runif(1) < a) 
    sigman = sigman.
  
  if (adaptive)
    Sigma = adapt_S(Sigma, u, a, jj, 0.44)
  
  list(sigman = sigman, acceptance = a, Sigma = Sigma)
  
}


#############
## Sampler ##
#############

sampler = function(n.samples, n.thin, y0, X0, y1, X1, z, Vz, kappa, 
                   n.theta, n.phi, n.zeta, initial.values, hyper.parameters, 
                   random.seed, Sigma, adaptive = FALSE, jj = 0, 
                   theta.names = NULL, phi.names = NULL) {
  
  ## Initialise progress bar
  pb = txtProgressBar(min = 0, max = n.samples, style = 3)
  
  ## Set seed
  assign(".Random.seed", random.seed, .GlobalEnv)
  
  ## Constants
  n.models = length(y0)
  
  ## Initialise storage
  phis    = array(NA, c(n.samples,n.phi))
  thetas  = array(NA, c(n.samples,n.theta*n.models))
  zetas   = array(NA, c(n.samples,n.zeta))
  nus     = array(NA, c(n.samples,ncol(X1)))
  sigmans = array(NA, n.samples)
  colnames(phis)   = phi.names
  colnames(thetas) = rep(theta.names, n.models)
  colnames(zetas)  = theta.names
  
  acceptance.theta  = array(NA, c(n.samples,n.models))
  acceptance.zeta   = array(NA, n.samples)
  acceptance.sigman = array(NA, n.samples)
  
  ## Initial values
  theta  = initial.values$theta
  phi    = initial.values$phi
  zeta   = initial.values$zeta
  nu     = initial.values$nu
  delta  = initial.values$delta
  sigman = initial.values$sigman
  
  ## Initialise sampling couter
  i = 0

  ## Sampling loop
  while (i < n.samples) {
    
    ## Increment sampling counter
    i = i + 1

    ## Thinning loop    
    for (j in 1:n.thin) {
    
      ## Increment adaptive counter
      jj = jj + 1
      
      ## theta
      theta = sample.thetas(theta, y0, X0, y1, X1, nu, phi, 
                            Sigma$theta, adaptive, jj)
      Sigma$theta = theta$Sigmas
      acc.theta   = theta$acceptances
      theta       = theta$thetas
      
      ## phi
      phi = sample.phi(phi, theta, hyper.parameters)
      
      ## delta
      delta = sample.deltas(delta, theta, nu, y1, X1)
      
      ## nu
      nu = sample.nu(nu, delta, theta, sigman)
      
      ## sigman
      sigman = sample.sigman(sigman, nu, hyper.parameters, 
                             Sigma$sigman, adaptive, jj)
      Sigma$sigman = sigman$Sigma
      acc.sigman   = sigman$acceptance
      sigman       = sigman$sigman
      
      ## zeta
      zeta = sample.zeta(zeta, z, Vz, X1, nu, phi, kappa, 
                         Sigma$zeta, adaptive, jj)
      Sigma$zeta = zeta$Sigma
      acc.zeta   = zeta$acceptance
      zeta       = zeta$zeta
      
    } ## j
    
    ## Store samples
    phis   [i,] = unlist.phi  (phi)
    thetas [i,] = unlist.theta(theta)
    zetas  [i,] = unlist.zeta (zeta)
    nus    [i,] = nu
    sigmans[i ] = sigman
    
    ## Store acceptance probabilities
    acceptance.theta [i,] = acc.theta
    acceptance.zeta  [i ] = acc.zeta
    acceptance.sigman[i ] = acc.sigman
    
    ## Update progress bar
    setTxtProgressBar(pb, i)
    
  } ## i
  
  ## Store final values
  final.values = list()
  final.values$theta  = theta
  final.values$phi    = phi
  final.values$zeta   = zeta
  final.values$nu     = nu
  final.values$delta  = delta
  final.values$sigman = sigman
  
  close(pb)
  
  ## Return results
  return(list(theta = thetas, phi = phis, zeta = zetas, nus = nus,
              sigman = sigmans, acceptance.theta = acceptance.theta, 
              acceptance.zeta = acceptance.zeta,
              acceptance.sigman = acceptance.sigman,
              random.seed = .Random.seed, Sigma = Sigma,
              initial.values = final.values))
  
} ## sampler

## Wrapper for running multiple chains simultaneously
wrapper = function(i, n.samples, n.thin, y0, X0, y1, X1, z, Vz, kappa,
                   n.theta, n.phi, n.zeta, initial.values, hyper.parameters, 
                   random.seeds, Sigma, adaptive, jj, 
                   theta.names, phi.names)
  sampler(n.samples, n.thin, y0, X0, y1, X1, z, Vz, kappa,
          n.theta, n.phi, n.zeta,
          initial.values[[i]], hyper.parameters, 
          random.seeds[[i]], Sigma[[i]], adaptive, jj, 
          theta.names, phi.names)
