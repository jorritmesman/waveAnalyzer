#'Calculate buoyancy frequency
#'
#'Calcululation of the buoyancy frequency, in rad/s
#'
#'@param wtr numeric; vector of water temperatures in degC
#'@param depths numeric; vector of depths in m, positive downward
#'@param dz numeric; depth spacing to interpolate to, defaults to 0.1 m
#'
#'@importFrom rLakeAnalyzer water.density
#'
#'@examples buoyancy_freq2(c(14, 13, 9, 8), c(0, 2, 4, 6))
#'
#'@export

buoyancy_freq2 = function(wtr, depths, dz = 0.1){
  g = 9.81
  
  # Interpolate onto regular grid of dz (linear interpolation)
  grid = seq(min(depths), max(depths), by = dz)
  wtr_grid = approx(x = depths, y = wtr, xout = grid)$y
  
  rho = water.density(wtr_grid)
  
  # Matlab:
  # % drhodz calculated as a center difference
  # drhodz = nan*ones(length(rho),1);
  # for i = 2:length(rho)-1
  #   drhodz(i) = (rho(i+1)-rho(i-1))/(zN(i+1)-zN(i-1));
  # end
  drhodz = diff(rho, lag = 2L) / diff(grid, lag = 2L)
  
  # Because z is positive downwards, drhodz should always be positive
  drhodz[drhodz < 0] = 0
  
  # return
  sqrt(g / mean(rho) * c(NA, drhodz, NA))
}
