#'Finds wave periods, thermocline depth, and metalimnion boundaries 
#'
#'Uses the temperature profile (z,T) and computes profile of buoyancy
#'frequency (zN,N) to solve for the V1 and V2 seiche modal structure. The
#'thermocline depth (zt) is defined as the depth of maximum displacement
#'during a V1 seiche and is the peak in the V1 modal structure. The
#'boundaries of the metalimnion (zm1, zm2) are the depths of maximum
#'displacement during a V2 seiche and are identified from the peaks in the
#'V2 modal structure. From this and the basin length L, periods of the V1H1
#'(T1) and V2H1 (T2) seiches (in s) are computed. This is based on Munnich et al
#'1992 and solves for the modes using many layers as opposed to a two or
#'three layer approximation.
#'Returns a list
#'
#'@param wtr numeric; vector of water temperatures in degC
#'@param depths numeric; vector of depths in m, positive downward
#'@param basin_length numeric; basin length in m
#'@param rLA_method logical; use rLakeAnalyzer's buoyancy frequency calculation,
#' defaults to FALSE
#'
#'@importFrom data.table nafill
#'@importFrom QZ qz.geigen
#'
#'@examples find_V1V2_modes(c(14, 13, 9, 8), c(0, 2, 4, 6), 7400)
#'
#'@export

find_V1V2_modes = function(wtr, depths, basin_length, rLA_method = FALSE){
  
  buoy_freq = buoyancy_freq(wtr, depths, dz = 0.1, rLA_method = rLA_method)
  
  # Potentially fill NAs
  buoy_freq$N = nafill(buoy_freq$N, type = "nocb")
  buoy_freq$N = nafill(buoy_freq$N, type = "locf")
  
  # Normalize N^2/mean(N)^2
  Nmean = mean(buoy_freq$N)
  gamma2 = (buoy_freq$N / Nmean)^2
  
  # Normalize z by max depth
  H = max(buoy_freq$z)
  D = buoy_freq$z / H
  
  # Solve the eigenvalue problem
  nseg = length(buoy_freq$N) + 1
  A = toeplitz(c(-2, 1, rep(0, nseg - 3)))
  delta = 1 / nseg
  B = diag(x = -delta^2 * gamma2)
  
  # The lambda calculation looks weird, but this seems to reproduce the
  # matlab code [phi, lambda] = eig(A, B)
  phi = qz.geigen(A = A, B = B)$VR
  lambda_vect = ((A %*% phi) / (B %*% phi))[1,]
  lambda = diag(lambda_vect)
  
  # Find two smallest eigenvalues, which correspond to V1 and V2
  lambda_vect[lambda_vect < 0] = NA # Not sure if this should be done?
  lmin = min(lambda_vect, na.rm = TRUE)
  index = which(lambda_vect == lmin)
  lambda1 = lmin # 1st vertical mode eigenvalue
  lambda_vect[index] = NA
  
  lmin = min(lambda_vect, na.rm = TRUE)
  index[2] = which(lambda_vect == lmin)
  lambda2 = lmin # 2nd vertical mode eigenvalue
  
  # wave speeds (m/s) from eigenvalues
  c1 = sqrt(1 / lambda1) * Nmean * H
  c2 = sqrt(1 / lambda2) * Nmean * H
  
  # period of the v1h1 mode (s)
  T1 = 2 * basin_length / c1
  # period of the v2h1 mode (s)
  T2 = 2 * basin_length / c2
  
  # phi contains information on the modal structure. Finding the minimum two
  # eigenvalues above, we have found which columns contain the modal
  # structure for those two modes. The peak in the V1 profile corresponds
  # with the thermocline and the peaks in the V2 profile correspond to the
  # boundaries of the metalimnion from an internal seiche energetics
  # standpoint as there are the depths of the maximum seiche amplitude. 
  
  phi1 = abs(phi[, index[1]])
  indx_zT = which.max(phi1)
  zt = buoy_freq$z[indx_zT]
  
  phi2 = (phi[, index[2]])
  indx_zm1 = which.max(phi2)
  temp1 = buoy_freq$z[indx_zm1]
  indx_zm2 = which.min(phi2)
  temp2 = buoy_freq$z[indx_zm2]
  zm1 = min(c(temp1, temp2))
  zm2 = max(c(temp1, temp2))
  
  return(list(N = buoy_freq$N,
              zN = buoy_freq$z,
              T1 = T1,
              T2 = T2,
              zt = zt,
              zm1 = zm1,
              zm2 = zm2))
}
