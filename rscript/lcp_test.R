remotes::install_github('josephlewis/leastcostpath')
library(sf)
library(elevatr)
library(dplyr)
library(terra)
library(leastcostpath)
library(purrr)
set.seed(123)

wgsLoc <- data.frame(x=-2.00715,y=51.38306) %>% st_as_sf(coords=c('x','y'),crs=4326)
origin <- st_transform(wgsLoc,27700)

dem <- 
  elevatr::get_elev_raster(wgsLoc,z = 8) %>% 
  as(.,'SpatRaster') %>% 
  terra::project(.,y = "epsg:27700")

destinations.point <- origin %>% 
  st_buffer(20000) %>% 
  st_bbox %>% 
  st_as_sfc %>% 
  st_sample(10) %>% 
  st_sf

destinations.multipoint <- 
  destinations.point %>% 
  summarise() %>% 
  st_cast('MULTIPOINT')

plot(dem)
plot(origin,add=T,col='red')
plot(destinations.point,add=T)

# destinations as single point df return lcp for first point
cs <- create_slope_cs(dem)
lcp <- create_lcp(cs, origin, destination=destinations.point)
plot(lcp$geometry,add=T)

# destinations as multipoint df return lcp for first point also
cs <- create_slope_cs(dem)
lcp <- create_lcp(cs, origin, destination=destinations.multipoint)
plot(lcp$geometry,add=T)

# can iterate over rows of single point df
system.time({
  lcps.orig <- 1:nrow(destinations.point) %>% 
    map(~create_lcp(cs, origin, destination=destinations.point[.x,])) %>% 
    bind_rows()
})

plot(lcps$geometry,add=T)


# from create_lcp()
system.time({
  x <- cs
  cs_rast <- terra::rast(nrow = x$nrow, ncol = x$ncol, xmin = x$extent[1], xmax = x$extent[2], ymin = x$extent[3], ymax = x$extent[4],crs = x$crs)
  
  from_coords <- sf::st_coordinates(origin)[1, 1:2, drop = FALSE]
  # to_coords <- sf::st_coordinates(destination)[1, 1:2, drop = FALSE]      # allow all rows of points df
  to_coords <- sf::st_coordinates(destinations.point)[, 1:2, drop = FALSE]
  
  from_cell <- terra::cellFromXY(cs_rast, from_coords)
  to_cell <- terra::cellFromXY(cs_rast, to_coords)
  
  cm_graph <- igraph::graph_from_adjacency_matrix(x$conductanceMatrix, mode = "directed", weighted = TRUE)
  
  igraph::E(cm_graph)$weight <- (1/igraph::E(cm_graph)$weight)
  
  lcp_graph <- igraph::shortest_paths(cm_graph, from = from_cell, to = to_cell, mode = "out")
  # lcp_cells <- unlist(lcp_graph$vpath) # need to iterate over multiple paths
  
  lcps.new <- 
    lcp_graph$vpath %>% 
    map(~terra::xyFromCell(cs_rast, as.integer(.x))) %>% 
    map(~sf::st_sf(geometry = sf::st_sfc(sf::st_linestring(.x)), crs = x$crs)) %>% 
    bind_rows()
})

lcps.new$geometry==lcps.orig$geometry



