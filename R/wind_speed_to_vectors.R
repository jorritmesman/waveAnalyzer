#'Wind speed and direction into vectors
#'
#'Function to convert wind speed and direction into vectors, returns vector
#'
#'@param speed numeric;
#'@param direction numeric; wind direction in degrees, 0 = north, 90 = east, etc
#'
#'@examples wind_speed_to_vectors(5, 45)
#'
#'@export

wind_speed_to_vectors = function(speed, direction){
  direction = 270 - direction # Converting the wind direction to the "math" direction
  rads = direction / 180 * pi
  xcomp = speed * cos(rads)
  ycomp = speed * sin(rads)
  wsvectors = c(u_wind = xcomp, v_wind = ycomp)
  return(wsvectors)
}
