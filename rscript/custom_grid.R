os50k <- vect('/media/mal/git/OS-British-National-Grids/os_bng_grids.gpkg',layer='50km_grid')
gbvect <- vect('/media/mal/git/lcp/bigdata/os_50m_polygon_outline.gpkg')

os50k.gb <- terra::mask(os50k, gbvect)
plot(os50k.gb)
text(os50k.gb, labels=os50k.gb$tile_name,cex=0.6)

1:nrow(os50k.gb) %>% 
  walk(~writeVector(os50k.gb[.x,],
                   glue::glue('/media/mal/git/OS-British-National-Grids/50km_grid/os_bng_grids_50km_{os50k.gb[.x,]$tile_name}.kml')))
library(terra)
sqrt(expanse(os50k.gb[1,])*2)
70750^2
ngrid <- gbvect %>% 
  st_as_sf() %>% 
  st_make_grid(cellsize=70700) %>% 
  vect() %>% 
  mask(gbvect)

1:nrow(ngrid) %>% 
  walk(~writeVector(ngrid[.x,],
                    glue::glue('/media/mal/git/lcp/bigdata/custom_grid/cg_{.x}.kml')),
       overwrite=T)

dir.create()
writeVector(ngrid[2,],'bigdata/custom_grid/cg1.kml',overwrite=T)


plot(gbvect)
plot(customG,add=T)
text(customG, labels=1:nrow(customG),cex=0.6)
customG <- vect('bigdata/sub5.shp')

1:nrow(customG.c) %>% 
  walk(~writeVector(customG.c[.x,],
                    glue::glue('/media/mal/git/lcp/bigdata/custom_grid/cg_{.x}.kml'),overwrite=T))


cg2 <- vect('bigdata/custom_grid/cg_2.kml')


