library(coda)

load.samples = function(base.name, load.dir, n.blocks, n.thin = 1, 
                        t.blocks = n.blocks) {
  
  load(paste0(load.dir,"/",base.name,"-",
              formatC(1, width = nchar(t.blocks), flag = 0), ".RData"))
  
  n.chains  = length(samples)
  n.samples = nrow(samples[[1]]$theta)
  n.theta   = ncol(samples[[1]]$theta)
  n.phi     = ncol(samples[[1]]$phi)
  n.zeta    = ncol(samples[[1]]$zeta)
  n.nus     = ncol(samples[[1]]$nus)
  n.sigman  = ncol(samples[[1]]$sigman)
  n.models  = ncol(samples[[1]]$acceptance.theta)
  
  mask = seq(from = n.thin, to = n.samples, by = n.thin)
  m.samples = length(mask)
  t.samples = n.blocks*m.samples
  a.samples = n.blocks*n.samples
  
  ## Initialize storage
  theta  = list()
  phi    = list()
  zeta   = list()
  nus    = list()
  sigman = list()
  acceptance.theta  = list()
  acceptance.zeta   = list()
  acceptance.sigman = list()
  for (i in 1:n.chains) {
    
    theta [[i]] = matrix(NA, t.samples, n.theta)
    phi   [[i]] = matrix(NA, t.samples, n.phi  )
    zeta  [[i]] = matrix(NA, t.samples, n.zeta )
    nus   [[i]] = matrix(NA, t.samples, n.nus  )
    sigman[[i]] = matrix(NA, t.samples, 1      )
    
    acceptance.theta [[i]] = matrix(NA, t.samples, n.models)
    acceptance.zeta  [[i]] = matrix(NA, t.samples, 1       )
    acceptance.sigman[[i]] = matrix(NA, t.samples, 1       )
    
  } ## i
  
  slice = 1:m.samples
  for (i in 1:n.chains) {
    
    theta [[i]][slice,] = samples[[i]]$theta [mask,]
    phi   [[i]][slice,] = samples[[i]]$phi   [mask,]
    zeta  [[i]][slice,] = samples[[i]]$zeta  [mask,]
    nus   [[i]][slice,] = samples[[i]]$nus   [mask,]
    sigman[[i]][slice ] = samples[[i]]$sigman[mask ]
    
    acceptance.theta [[i]][slice,] = samples[[i]]$acceptance.theta [mask,]
    acceptance.zeta  [[i]][slice ] = samples[[i]]$acceptance.zeta  [mask ]
    acceptance.sigman[[i]][slice ] = samples[[i]]$acceptance.sigman[mask ]
    
  } ## i
  
  ## Load samples
  if (n.blocks > 1) {
    for (b in 2:n.blocks) {
      
      load(paste0(load.dir,"/",base.name,"-",
                  formatC(b, width = nchar(t.blocks), flag = 0), ".RData"))
      
      slice = ((b - 1)*m.samples + 1):(b*m.samples)
      for (i in 1:n.chains) {
        
        theta [[i]][slice,] = samples[[i]]$theta [mask,]
        phi   [[i]][slice,] = samples[[i]]$phi   [mask,]
        zeta  [[i]][slice,] = samples[[i]]$zeta  [mask,]
        nus   [[i]][slice,] = samples[[i]]$nus   [mask,]
        sigman[[i]][slice ] = samples[[i]]$sigman[mask ]
        
        acceptance.theta [[i]][slice,] = samples[[i]]$acceptance.theta [mask,]
        acceptance.zeta  [[i]][slice ] = samples[[i]]$acceptance.zeta  [mask ]
        acceptance.sigman[[i]][slice ] = samples[[i]]$acceptance.sigman[mask ]
        
      } ## i
    } ## b
  }
  
  for (i in 1:n.chains) {
    
    theta [[i]] = mcmc(theta [[i]], n.thin, a.samples, n.thin)
    phi   [[i]] = mcmc(phi   [[i]], n.thin, a.samples, n.thin)
    zeta  [[i]] = mcmc(zeta  [[i]], n.thin, a.samples, n.thin)
    nus   [[i]] = mcmc(nus   [[i]], n.thin, a.samples, n.thin)
    sigman[[i]] = mcmc(sigman[[i]], n.thin, a.samples, n.thin)
    
  } ## i
  
  theta  = mcmc.list(theta )
  phi    = mcmc.list(phi   )
  zeta   = mcmc.list(zeta  )
  nus    = mcmc.list(nus   )
  sigman = mcmc.list(sigman)
  
  return(list(theta = theta, phi = phi, zeta = zeta, nus = nus, sigman = sigman,
              acceptance.theta  = acceptance.theta, 
              acceptance.zeta   = acceptance.zeta,
              acceptance.sigman = acceptance.sigman))
  
}

