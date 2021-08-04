#read in ncss download and set projection.
test <- raster::stack("W:/GIS/Daymet/raw/1840_daymet_v4_daily_na_tmax_2019.nc") 
raster::projection(test) <- "+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"

#create bounding box using exact same coordinates in request URL:
df <- data.frame(north_lat = 64.145, south_lat = 56.0596, west_lng = -163.7336, east_lng = -141.9427)
poly_df <- matrix(c(df[1, 'west_lng'], df[1, 'north_lat'], 
                    df[1, 'east_lng'], df[1, 'north_lat'], 
                    df[1, 'east_lng'], df[1, 'south_lat'], 
                    df[1, 'west_lng'], df[1, 'south_lat'],
                    df[1, 'west_lng'], df[1, 'north_lat'])  ## need to close the polygon
                  , ncol =2, byrow = T) 
poly <- st_polygon(list(poly_df))
poly <- st_sfc(poly)
st_crs(poly) <- 4326

#convert to same lcc projection as ncss download.
poly_llc <- st_transform(poly, crs = st_crs(test))

#plot the two together.
plot(test[[1]])
plot(poly_llc, add = TRUE, color = "red")