#'Calculate buoyancy frequency
#'
#'Calculation of the buoyancy frequency "N", in rad/s, returns data.frame
#'
#'@param wtr numeric; vector of water temperatures in degC
#'@param depths numeric; vector of depths in m, positive downward
#'@param dz numeric; depth spacing to interpolate to, defaults to 0.1 m
#'@param rLA_method logical; use rLakeAnalyzer method to calculate z, defaults to FALSE
#'
#'@importFrom rLakeAnalyzer water.density buoyancy.freq
#'
#'@examples buoyancy_freq(c(14, 13, 9, 8), c(0, 2, 4, 6))
#'
#'@export

buoyancy_freq = function(wtr, depths, dz = 0.1, rLA_method = FALSE){
  
  if(rLA_method){
    returned_obj = buoyancy.freq(wtr, depths)
    return(data.frame(z = attr(returned_obj, "depths"),
                      N = sqrt(returned_obj)))
  }
  
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
  data.frame(z = grid,
             N =sqrt(g / mean(rho) * c(NA, drhodz, NA)))
}
