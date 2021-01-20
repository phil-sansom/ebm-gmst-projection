library(MASS)
library(nloptr)
library(numDeriv)

## MCMC setup
set.seed(0)

initial.values = list()
Sigma = list()
for (i in 1:n.chains) {
  
  initial.values[[i]] = list()
  initial.values[[i]]$theta  = list()

  Sigma[[i]] = list()
  Sigma[[i]]$theta  = list()

} ## i

linits = log(c(2,5,20,100,1,2,1,1,0.5,0.5,3,20,0.05))
fn = function(theta, y0, X0, y1, X1, nu, phi) {
  
  - theta.likelihood(theta, y0, X0, y1, X1, nu) - 
    theta.prior     (theta, phi) -
    theta.jacobian  (theta, phi)
  
}

delta0 = matrix(0, n.models, ncol(X1))
nu0 = rep(0, ncol(X1))
phi0 = list(mu = log(c(2,5,20,100,1,2,1,1,0.5,0.5,3,20,0.05)),
            Sigma = diag(rep(1/(log(10)/3)^2, n.theta)))

for (m in 1:n.models) {
  
  ## theta
  map = bobyqa(linits, fn, y0 = y0[[m]], X0 = X0, y1 = y1[[m]], 
               X1 = X1, nu = nu0, phi = phi0, 
               control = list(maxeval = 2e4))$par
  vcov = solve(hessian(func = fn, x = map, y = y0[[m]], X = X0, 
                       y1 = y1[[m]], X1 = X1, nu = nu0, phi = phi0))
  vcov = 0.5*(vcov + t(vcov))
  Sigma0 = t(chol(2.38^2 * vcov / length(map)))
  
  for (i in 1:n.chains) {
    
    initial.values[[i]]$theta[[m]] = mvrnorm(1, map, vcov)
    Sigma[[i]]$theta[[m]] = Sigma0
    
  } ## i
  
  rm(map,Sigma0,vcov)
  
} ## m

fn = function(zeta, z, V, X, nu, phi, kappa) {
  
  - zeta.likelihood(zeta, z, V, X, nu) - 
    zeta.prior     (zeta, phi, kappa) -
    zeta.jacobian  (zeta, phi, kappa)
  
}

linits1 = log(0.1)
gn = function(sigman, nu, hyper.parameters) {
  
  - sigman.likelihood(sigman, nu) - 
    sigman.prior(sigman, hyper.parameters) -
    sigman.jacobian(sigman, hyper.parameters)
  
}

for (i in 1:n.chains) {
  
  ## phi
  initial.values[[i]]$phi = 
    sample.phi(phi0, initial.values[[i]]$theta, hyper.parameters)
  
  ## delta
  initial.values[[i]]$delta = 
    sample.deltas(delta0, initial.values[[i]]$theta, nu0, y1, X1)

  ## nu
  nu = apply(initial.values[[i]]$delta, 2, mean)
  nu = diff(c(0,nu))
  sigman = log(sd(nu))
  initial.values[[i]]$nu = 
    sample.nu(nu, initial.values[[i]]$delta, initial.values[[i]]$theta, sigman)
  
  ## sigman
  map = bobyqa(linits1, gn, nu = initial.values[[i]]$nu, 
               hyper.parameters = hyper.parameters,
               control = list(maxeval = 2e4))$par
  vcov = solve(hessian(func = gn, x = map, 
                       nu = initial.values[[i]]$nu, 
                       hyper.parameters = hyper.parameters,))
  vcov = 0.5*(vcov + t(vcov))
  Sigma[[i]]$sigman = t(chol(2.38^2 * vcov / length(map)))
  initial.values[[i]]$sigman = rnorm(1, map, sqrt(vcov))

  rm(map,vcov)
  
  ## zeta
  map = bobyqa(linits, fn, z = z, V = Vz, X = X1, 
               nu  = initial.values[[i]]$nu, phi = initial.values[[i]]$phi,
               kappa = kappa, control = list(maxeval = 2e4))$par
  vcov = solve(hessian(func = fn, x = map, z = z, V = Vz, X = X1, 
                       nu  = initial.values[[i]]$nu, 
                       phi = initial.values[[i]]$phi, kappa = kappa))
  vcov = 0.5*(vcov + t(vcov))
  Sigma[[i]]$zeta = t(chol(2.38^2 * vcov / length(map)))
  initial.values[[i]]$zeta = mvrnorm(1, map, vcov)
  
  rm(map,vcov)
  
} ## i
rm(fn,gn,linits,linits1,delta0,nu0,nu,sigman)

## Random seeds
random.seeds   = list()
for (i in 1:n.chains) {
  set.seed(i)
  random.seeds[[i]] = .Random.seed
} ## i
