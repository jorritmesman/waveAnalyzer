#'Average wind vectors
#'
#'Averages wind vectors but retains average absolute wind speed
#'
#'@param u_wind numeric; vector of values
#'@param v_wind numeric; vector of values
#'@param na.rm logical; remove NAs or not
#'
#'@examples average_wind_vectors(u_wind = c(2, -2), v_wind = c(0, 4))
#'
#'@export

# Returned direction can be NA if vectors cancel each other out completely
# e.g. average_wind_vectors(u_wind = c(2, -2), v_wind = c(-2, 2))
average_wind_vectors = function(u_wind, v_wind, na.rm = FALSE){
  if(length(u_wind) != length(v_wind)){
    stop("u_wind and v_wind do not have the same length!")
  }
  
  if(na.rm){
    ind_na_u = is.na(u_wind)
    ind_na_v = is.na(v_wind)
    ind_to_use = !ind_na_u & !ind_na_v
    u_wind = u_wind[ind_to_use]
    v_wind = v_wind[ind_to_use]
  }else{
    if(any(is.na(c(u_wind, v_wind)))){
      stop("The data contains NA values!")
    }
  }
  
  abs_wind = sqrt(u_wind^2 + v_wind^2)
  
  av_u = mean(u_wind)
  av_v = mean(v_wind)
  
  speed_and_dir = wind_vectors_to_speed(av_u, av_v)
  
  return(c(speed = mean(abs_wind), direction = speed_and_dir[["direction"]]))
}
