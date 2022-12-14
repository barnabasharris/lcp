
writeLines("ncols        11
nrows        11
xllcorner    173550.000000000000
yllcorner    26250.000000000000
cellsize     50.000000000000
NODATA_value  nan
 nan 7.1999998092651367188 10 15 20 25 30 33.40000152587890625 33.90000152587890625 33.90000152587890625 nan
 14.699999809265136719 4.8000001907348632812 -1.6000000238418579102 -1.6000000238418579102 nan nan nan 26.399999618530273438 27.200000762939453125 26.600000381469726562 26.600000381469726562
 5.0999999046325683594 1.7999999523162841797 -1.6000000238418579102 -1.6000000238418579102 nan nan nan 16.700000762939453125 18.200000762939453125 18.100000381469726562 18
 3.2000000476837158203 1.7999999523162841797 0.5 -1.6000000238418579102 nan nan nan nan nan 3.2999999523162841797 3
 8.3999996185302734375 6.6999998092651367188 2.7000000476837158203 -1.6000000238418579102 -1.6000000238418579102 nan nan nan nan nan nan
 16.5 12.300000190734863281 7.5 1.7999999523162841797 -1.6000000238418579102 -1.6000000238418579102 -1.6000000238418579102 nan nan nan nan
 23.899999618530273438 19.299999237060546875 15 10.899999618530273438 4.1999998092651367188 -1.6000000238418579102 -1.6000000238418579102 -1.6000000238418579102 -1.6000000238418579102 nan nan
 29.799999237060546875 26 21.700000762939453125 20 17 14.5 12.300000190734863281 11 9.6999998092651367188 7.1999998092651367188 3.4000000953674316406
 29.299999237060546875 26.100000381469726562 23.200000762939453125 21.100000381469726562 19.299999237060546875 17.5 15.600000381469726562 15 14 12.899999618530273438 12.399999618530273438
 24.700000762939453125 20.200000762939453125 16.5 16.5 14.800000190734863281 12.699999809265136719 9.3000001907348632812 6.4000000953674316406 3.5 8.1000003814697265625 9.6000003814697265625
 nan 3.5 -1.6000000238418579102 -1.6000000238418579102 -1.6000000238418579102 -1.6000000238418579102 -1.6000000238418579102 -1.6000000238418579102 -1.6000000238418579102 -1.6000000238418579102 nan
",'r.asc')

origin <- 
  vect(data.frame(x = 173775,
                  y = 26565),
       geom=c('x','y'),
       crs='EPSG:27700') %>% 
  st_as_sf

# origin <- 
#   vect(data.frame(x = 178557,
#                   y = 36776),
#        geom=c('x','y'),
#        crs='EPSG:27700') %>% 
#   st_as_sf

dest <- 
  vect(data.frame(x = 173985,
                  y = 26720),
       geom=c('x','y'),
       crs='EPSG:27700') %>% 
  st_as_sf

# dest <- 
#   vect(data.frame(x = 190725,
#                   y = 46344),
#        geom=c('x','y'),
#        crs='EPSG:27700') %>% 
#   st_as_sf



r <- hex.r
r <- rast('r.asc')
# r <- rast('test_dem.tif')
neighbours32 <- neighbourhood(neighbours = 32)

cells <- which(!is.na(terra::values(r)))
na_cells <- which(is.na(terra::values(r)))

adj <- terra::adjacent(x = r, cells = cells, 
                       directions = neighbours32, pairs = TRUE)
adj <- adj[!adj[,2] %in% na_cells,]


# visualize
r.all <- x <- r
r.all[] <- 1
r.v <- as.polygons(r.all,dissolve=F,values=F)
r.v$ID <- 1:ncell(x)
plot(x)
plot(r.v,add=T)
plot(r.v[dangerLinks$ID,],add=T,col='red')
text(r.v,label=r.v$ID)
r <- rast('r.asc')

x <- r

gb <- rast('bigdata/os_50m_masked.tif')
riv <- rast('bigdata/navigable_rivers.tif')

gb[riv==1] <- NaN
writeRaster(gb,'bigdata/os_50m_masked_w_navigable_rivers.tif')
plot(gb)

x <- rast('r.asc')

cs <- create_slope_cs2(x,neighbours = 32)
lcp <- create_lcp(cs, origin, dest)
x.bin <- cs$x
x.bin[is.nan(x.bin)] <- 9999
plot(x.bin)
plot(lcp$geometry,add=T,col='blue')
st_write(lcp,'lcp3.gpkg')
writeRaster(x.bin,'xbin.tif')






writeLines("ncols        5
nrows        5
xllcorner    173550.0
yllcorner    26250.0
cellsize     50.0
NODATA_value  nan
 1.0 2.0 3.0 4.0 5.0 6.0 7.0	nan	9.0	10.0 11.0	12.0	13.0	14.0	15.0 16.0	17.0	18.0	19.0	20.0 21.0	22.0	23.0	24.0	25.0
",'r.test.asc')

r <- rast('r.test.asc')
plot(r)
neighbours32 <- neighbourhood(neighbours = 32)

sum(neighbours32)
cells <- which(!is.na(terra::values(r)))
na_cells <- which(is.na(terra::values(r)))

adj <- terra::adjacent(x = r, cells = 13, 
                       directions = 16)

adj <- adj[!adj[,2] %in% na_cells,]


# visualize
r.all <- x
r.all[] <- 1
r.v <- as.polygons(r.all,dissolve=F,values=F)
r.v$ID <- 1:ncell(x)
plot(x)
plot(v,add=T)
plot(v[dangerLinks$ID,],add=T,col='red')
text(r.v,label=r.v$ID)
r <- rast('r.asc')
