betasw_ZHH2009 <- function(lambda, Tc, theta, S, delta = 0.039) {
  # Constants
  Na <- 6.0221417930e23   # Avogadro's constant
  Kbz <- 1.3806503e-23    # Boltzmann constant
  Tk <- Tc + 273.15       # Absolute temperature
  M0 <- 18e-3             # Molecular weight of water in kg/mol
  
  if (!is.numeric(Tc) || length(Tc) != 1 || !is.numeric(S) || length(S) != 1) {
    stop("Both Tc and S need to be scalar variables.")
  }
  
  lambda <- as.numeric(lambda)  # Ensure lambda is numeric
  rad <- theta * pi / 180       # Convert angle to radians
  
  # Helper functions
  RInw <- function(lambda, Tc, S) {
    n_air <- 1.0 + (5792105.0 / (238.0185 - 1 / (lambda / 1000)^2) +
                      167917.0 / (57.362 - 1 / (lambda / 1000)^2)) / 1e8
    n0 <- 1.31405
    n1 <- 1.779e-4
    n2 <- -1.05e-6
    n3 <- 1.6e-8
    n4 <- -2.02e-6
    n5 <- 15.868
    n6 <- 0.01155
    n7 <- -0.00423
    n8 <- -4382
    n9 <- 1.1455e6
    nsw <- n0 + (n1 + n2 * Tc + n3 * Tc^2) * S + n4 * Tc^2 +
      (n5 + n6 * S + n7 * Tc) / lambda + n8 / lambda^2 + n9 / lambda^3
    nsw <- nsw * n_air
    dnswds <- (n1 + n2 * Tc + n3 * Tc^2 + n6 / lambda) * n_air
    return(list(nsw = nsw, dnswds = dnswds))
  }
  
  BetaT <- function(Tc, S) {
    kw <- 19652.21 + 148.4206 * Tc - 2.327105 * Tc^2 +
      1.360477e-2 * Tc^3 - 5.155288e-5 * Tc^4
    a0 <- 54.6746 - 0.603459 * Tc + 1.09987e-2 * Tc^2 - 6.167e-5 * Tc^3
    b0 <- 7.944e-2 + 1.6483e-2 * Tc - 5.3009e-4 * Tc^2
    Ks <- kw + a0 * S + b0 * S^1.5
    return(1 / Ks * 1e-5)  # Unit is Pa
  }
  
  rhou_sw <- function(Tc, S) {
    a0 <- 8.24493e-1
    a1 <- -4.0899e-3
    a2 <- 7.6438e-5
    a3 <- -8.2467e-7
    a4 <- 5.3875e-9
    a5 <- -5.72466e-3
    a6 <- 1.0227e-4
    a7 <- -1.6546e-6
    a8 <- 4.8314e-4
    b0 <- 999.842594
    b1 <- 6.793952e-2
    b2 <- -9.09529e-3
    b3 <- 1.001685e-4
    b4 <- -1.120083e-6
    b5 <- 6.536332e-9
    
    density_w <- b0 + b1 * Tc + b2 * Tc^2 + b3 * Tc^3 + b4 * Tc^4 + b5 * Tc^5
    density_sw <- density_w + ((a0 + a1 * Tc + a2 * Tc^2 + a3 * Tc^3 + a4 * Tc^4) * S +
                                 (a5 + a6 * Tc + a7 * Tc^2) * S^1.5 + a8 * S^2)
    return(density_sw)
  }
  
  dlnasw_ds <- function(Tc, S) {
    (-5.58651e-4 + 2.40452e-7 * Tc - 3.12165e-9 * Tc^2 + 2.40808e-11 * Tc^3) +
      1.5 * (1.79613e-5 - 9.9422e-8 * Tc + 2.08919e-9 * Tc^2 - 1.39872e-11 * Tc^3) * S^0.5 +
      2 * (-2.31065e-6 - 1.37674e-9 * Tc - 1.93316e-11 * Tc^2) * S
  }
  
  PMH <- function(nsw) {
    nsw2 <- nsw^2
    (nsw2 - 1) * (1 + 2 / 3 * (nsw2 + 2) * (nsw / 3 - 1 / (3 * nsw))^2)
  }
  
  # Refractive index of seawater and partial derivative w.r.t salinity
  ref_data <- RInw(lambda, Tc, S)
  nsw <- ref_data$nsw
  dnds <- ref_data$dnswds
  
  # Isothermal compressibility from Lepple & Millero
  IsoComp <- BetaT(Tc, S)
  
  # Density of seawater
  density_sw <- rhou_sw(Tc, S)
  
  # Water activity derivative
  dlnawds <- dlnasw_ds(Tc, S)
  
  # PMH model for density derivative
  DFRI <- PMH(nsw)
  
  # Volume scattering at 90 degrees due to density fluctuation
  beta_df <- pi^2 / 2 * (lambda * 1e-9)^(-4) * Kbz * Tk * IsoComp * DFRI^2 * (6 + 6 * delta) / (6 - 7 * delta)
  
  # Volume scattering at 90 degrees due to concentration fluctuation
  flu_con <- S * M0 * dnds^2 / density_sw / (-dlnawds) / Na
  beta_cf <- 2 * pi^2 * (lambda * 1e-9)^(-4) * nsw^2 * flu_con * (6 + 6 * delta) / (6 - 7 * delta)
  
  # Total volume scattering at 90 degrees
  beta90sw <- beta_df + beta_cf
  
  # Total scattering coefficient
  bsw <- 8 * pi / 3 * beta90sw * (2 + delta) / (1 + delta)
  
  # Volume scattering at angles defined by theta
  betasw <- sapply(lambda, function(l) {
    beta90sw * (1 + (cos(rad)^2) * (1 - delta) / (1 + delta))
  })
  
  return(list(betasw = betasw, beta90sw = beta90sw, bsw = bsw))
}
