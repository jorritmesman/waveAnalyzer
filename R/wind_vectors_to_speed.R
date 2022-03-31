#'Wind vectors to speed and dir
#'
#'Function to convert wind vectors into speed and direction in degrees, returns vector
#' Returned direction: 0 is north, 90 is east, etc
#'
#'@param u_wind numeric;
#'@param v_wind numeric;
#'
#'@examples wind_vectors_to_speed(-2, 3)
#'
#'@export

# Returns NA when two zeroes are given
wind_vectors_to_speed = function(u_wind, v_wind){
  windspeed = sqrt(u_wind^2 + v_wind^2)
  wind_dir_rad = atan2(v_wind, u_wind) 
  wind_dir_degrees = wind_dir_rad * 180/pi
  speed_and_dir = c(speed = windspeed, direction = 270 - wind_dir_degrees)
  if(speed_and_dir["direction"] >= 360){
    speed_and_dir["direction"] = speed_and_dir["direction"] - 360
  }
  if(u_wind == 0 & v_wind == 0){
    speed_and_dir["direction"] = NA
  }
  
  return(speed_and_dir)
}
