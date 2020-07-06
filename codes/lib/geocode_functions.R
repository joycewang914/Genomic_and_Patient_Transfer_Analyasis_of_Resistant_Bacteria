# V_long_lat - 
# Input: 
#   1) google_geocode object (as a matrix)
# Returns the longitude and latitude coordinates of an address

V_long_lat = function(address){
  temp_loc = apply(address, 2, FUN = function(x){ x[['results']][,"geometry"][grep("location", colnames(x$results[,"geometry"]))][1][1,]} )
  loc_df = do.call(rbind, temp_loc)
  loc_df = loc_df[,2:1]
  return(loc_df)
}
# End V_long_lat



