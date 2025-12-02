# Function to simulate data under single-state Cormack-Jolly-Seber (CJS) model

sim_cjs <- function(
    T,    # number of sampling occasions
    n,    # number of new releases per occasion: vector of length T-1 or scalar 
    phi,  # survival probabilities: vector of length T-1, or scalar
    p     # detection probabilities: vector of length T, or scalar
) {
  # checks
  if (length(T) != 1) {
    stop("Number of sampling occasions must be an integer.")
  }
  if (as.integer(T) != T) {
    stop("Number of sampling occasions must be an integer.")
  }
  if (T < 2) {
    stop("At least two sampling occasions required.")
  }
  if ((length(n) != 1) & (length(n) != T-1)) {
    stop("Number of new releases per sampling occasion must be an integer, or a vector of T-1 integers.")
  }
  if (!all(as.integer(n) == n)) {
    stop("Number of new releases per sampling occasion must be an integer, or a vector of T-1 integers.")
  }
  if (!all(n >= 0)) {
    stop("Number of new releases per sampling occasion must be non-negative integer(s).")
  }
  if ((length(phi) != 1) & (length(phi) != T-1)) {
    stop("Survival probability must be a scalar, or a vector of length T-1.")
  }
  if (!all((0 <= phi) & (phi <= 1))) {
    stop("Survival probabilities must be between 0 and 1.")
  }
  if ((length(p) != 1) & (length(phi) != T)) {
    stop("Detection probability must be a scalar, or a vector of length T.")
  }
  if (!all((0 <= p) & (p <= 1))) {
    stop("Detection probabilities must be between 0 and 1.")
  }
  
  # promote scalar arguments to constant vectors
  if (length(n) == 1) {
    n <- rep(n, T-1)
  }
  if (length(phi) == 1) {
    phi <- rep(phi, T-1)
  }
  if (length(p) == 1) {
    p <- rep(p, T)
  }
  
  # initialise matrix of true states
  n_tot <- sum(n)
  fc <- rep(1:(T-1), n)  # first capture
  z <- matrix(NA, nrow = n_tot, ncol = T)
  for (i in 1:n_tot) {
    z[i, fc[i]] <- 1
  }
  
  # create transition matrices 
  tr <- list()
  for (t in 1:(T-1)) {
    tr[[t]] <- matrix(
      c(phi[t], 1 - phi[t],
             0,          1), 
      nrow = 2, ncol = 2, byrow = TRUE)
  }
  
  # simulate true states
  for (i in 1:n_tot) {
    for (t in fc[i]:(T-1)) {
      z[i,t+1] <- sample(1:2, 1, prob = tr[[t]][z[i,t], ])
    }
  }
  
  # create observation matrices 
  obs <- list()
  for (t in 1:T) {
    obs[[t]] <- matrix(
      c(p[t], 1 - p[t],
           0,        1), 
      nrow = 2, ncol = 2, byrow = TRUE)
  }
  
  # simulate capture histories
  y <- matrix(0, nrow = n_tot, ncol = T)
  for (i in 1:n_tot) {
    y[i, fc[i]] <- 1  # observed on first capture occasion
    for (t in (fc[i] + 1):T) {
      y[i,t] <- sample(1:2, 1, prob = obs[[t]][z[i,t],])
    }
  }
  
  # create alternative version of capture histories 
  # with non-detections after first capture coded as '0' rather than '2'
  y_zero <- y
  y_zero[y_zero == 2] <- 0
  
  return(list(T = T, n = n, phi = phi, p = p, z = z, y = y, y_zero = y_zero))
}
