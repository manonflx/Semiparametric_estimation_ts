Whittle <- function(theta, phi, H, p, q, y){
  
  eta = numeric(length = (1 + p + q))
  eta = cbind(H, theta, phi)
  H <- eta[1]
  m = length(y)
  mhalfm <- (m-1) %/% 2L
  ### Fundamental frequencies 
  x <- 2*pi/m * (1:mhalfm)
  
  ### Periodogram 
  n <- length(y)
  per = (Mod(fft(y))^2/(2*pi*n)) [1:(n %/% 2 + 1)]
  periodogr.x = per[2:((n+1) %/% 2)]
  
  
  ### Spectral density
  if(p > 0) {
    px <- outer(x, 1:p)
    Rar <- cos(px) %*% theta
    Iar <- sin(px) %*% theta
    
    far <- (1-Rar)^2 + Iar^2
  } else {
    theta <- numeric(0)
    far <- 1
  }
  
  if(q > 0) {
    psi <- phi
    px <- outer(x, 1:q)
    Rma <- cos(px) %*% psi
    Ima <- sin(px) %*% psi
    
    fma <- (1+Rma)^2 + Ima^2
  } else {
    psi <- numeric(0)
    fma <- 1
  }
  
  spec <- fma/far * sqrt(2 - 2*cos(x))^(1-2*H)
  
  yf <- periodogr.x/spec
  B <- 2*(2*pi/n) * sum(yf) # to be optimized
}


Wasserstein <- function(theta, phi, H, p, q, y, theoritical, cdf = FALSE, mean = TRUE, k = 1 , weighted = FALSE, sinkhorn = FALSE, lambda = 1){
  
  eta = numeric(length = (1 + p + q))
  eta = cbind(H, theta, phi)
  H <- eta[1]
  m = length(y)
  mhalfm <- (m-1) %/% 2L
  ### Fundamental frequencies 
  x <- 2*pi/m * (1:mhalfm)
  
  ### Periodogram 
  n <- length(y)
  per = (Mod(fft(y))^2/(2*pi*n)) [1:(n %/% 2 + 1)]
  periodogr.x = per[2:((n+1) %/% 2)]
  
  
  ### Spectral density
  if(p > 0) {
    px <- outer(x, 1:p)
    Rar <- cos(px) %*% theta
    Iar <- sin(px) %*% theta
    
    far <- (1-Rar)^2 + Iar^2
  } else {
    theta <- numeric(0)
    far <- 1
  }
  
  if(q > 0) {
    psi <- phi
    px <- outer(x, 1:q)
    Rma <- cos(px) %*% psi
    Ima <- sin(px) %*% psi
    
    fma <- (1+Rma)^2 + Ima^2
  } else {
    psi <- numeric(0)
    fma <- 1
  }
  
  spec <- fma/far * sqrt(2 - 2*cos(x))^(1-2*H)
  
  ### Standardized periodogram ordinates  
  yper = periodogr.x
  yf <- 2*pi * yper/spec # follow asympt. an Exponential (1) distribution
  
  ### Compute the Wasserstein distance 1D 
  if (weighted == TRUE){
    weight = sort(yf) / sum(yf)
    WD = transport::wasserstein1d(theoritical, yf,  wa = NULL, wb = weight, p = 1)
  } else if (sinkhorn == TRUE) {
    x = as.matrix(yf/ sum(yf))
    z = as.matrix(theoritical/sum(theoritical))
    m = length(yf)
    size <- seq(0,1,length.out = sqrt(m))
    costm <- as.matrix(dist(expand.grid(size,rev(size)), diag=TRUE, upper=TRUE))
    
    SH = Barycenter::Sinkhorn(x, z, costm, lambda = lambda)$Distance
    
    #SH = T4transport::sinkhorn(theoritical, yf, p = 1, lambda = lambda)$distance
    
  } else if (cdf == TRUE) {
    ## Empirical CDF values 
    fun.ecdf <- ecdf(yf)
    my.ecdf <- fun.ecdf(sort(yf))
    
    WD =  pracma::trapz(sort(yf) ,  abs(1- exp(-sort(yf)) - my.ecdf))
    
  
  } else if (mean == TRUE) {
    ## Empirical CDF values 
    WD = numeric(k)
    for (i in 1:k){
      t = rexp(((n-1)/2), rate = 1)
      WD[i] = transport::wasserstein1d(t, yf,  wa = NULL, wb = NULL, p = 1)
    }
    
    WD = mean(WD)

  } else {
    WD = transport::wasserstein1d(theoritical, yf,  wa = NULL, wb = NULL, p = 1)
  } 
  
}

