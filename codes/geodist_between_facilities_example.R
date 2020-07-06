# Demonstrates how distance between nursing facilities is calculated

# Start example

# Geocode of address of interest can be retrieved with Google API
# e.g. some hospitals
UM_address = "1500 E Medical Center Dr, Ann Arbor, MI 48109, United States"
Rush_address = "1620 W Harrison St, Chicago, IL 60612, United States"
Mayo_address = "200 First St. SW Rochester, MN 55905, United States"
Cleveland_address = "11100 Euclid Ave, Cleveland, OH 44106, United States"

example_hospitals = structure(c(UM_address, Rush_address, Mayo_address, Cleveland_address), names = c("Michigan", "Rush", "Mayo", "Cleveland"))

api_key = "" # Obtained from Google Cloud Platform
hospital_geocode = sapply(example_hospitals, FUN = function(x){print(x)
  google_geocode(address = x, key = api_key)
})

# code above generate this result:
hospital_geocode = readRDS("hospital_geocode.RDS")
hospital_loc = V_long_lat(hospital_geocode)

hospital_dist = 0.001*outer(1:nrow(hospital_loc), 1:nrow(hospital_loc),
                            Vectorize(FUN = function(x,y){
                              distm(c(hospital_loc[x,"lng"], hospital_loc[x,"lat"]),
                                    c(hospital_loc[y,"lng"], hospital_loc[y,"lat"]), fun = distHaversine)})) #meters -> km
rownames(hospital_dist) = rownames(hospital_loc)
colnames(hospital_dist) = rownames(hospital_loc)

### End example