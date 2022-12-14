

os50m.v <- st_read('bigdata/os_50m_polygon_outline.gpkg')

polgrid <- st_make_grid(os50m.v, cellsize = 5000) %>% 
  vect()

polgrid.osgb <- mask(polgrid,vect(os50m.v))
polgrid.osgb.ag <- aggregate(polgrid.osgb,dissolve=T)
writeVector(polgrid.osgb.ag,'bigdata/gbgrid.gpkg')

# import polygon maximising tiles from Edina
poldiv.sf <- st_read('~/temp_stuff/gbgrid_49_divider.shp')
poldiv <- poldiv.sf %>% 
  arrange(POINTY) %>% 
  mutate(area = as.numeric(st_area(.))) %>% 
  filter(area > 1000000000) %>% 
  mutate(fid = 1:nrow(.)) %>% 
  vect()
poldiv$AREA <- NULL

plot(poldiv)
text(poldiv, poldiv$fid)

poldiv$fid %>% 
  map(~writeVector(poldiv[poldiv$fid==.x],glue::glue('bigdata/sw{.x}.kml'),
                   overwrite=T))
writeVector(poldiv,'bigdata/custom_grid/allpol.gpkg')
list.files('bigdata',pattern='^sw*',full.names = T)

# use to download from Edina

####

# process mm water network downloads
wat.v <- 
  list.files('/media/mal/os_master_map_water_network',
             full.names = T,pattern='*gz$',
             recursive = T) %>% 
  map(~tryCatch({vect(.x,layer='WatercourseLink')},
                error = function(cond) return(NULL)))

wat.v.all <- vect(purrr::compact(wat.v))

dem <- rast('bigdata/os_50m_masked.tif')

# first, remove any links with widths as NaN. This means rivers with widths less than 5m, 3m, or 1m will
# be removed depending on their location. With worst accuracy in moors and mountains.
wat.v.havewidth <- wat.v.all[wat.v.all$width != 'nan',]
wat.v.havewidth$width <- as.numeric(wat.v.havewidth$width)
# second, remove any waterway that isn't permanent
wat.v.havewidth.perm <- wat.v.havewidth[wat.v.havewidth$permanence,]
# third remove any water way less than 1m across
wat.v.havewidth.perm.over1m <- wat.v.havewidth.perm[wat.v.havewidth.perm$width > 1,]
# saved as 'bigdata/sw_mm_water_processed.gpkg'
wat.v.havewidth.perm.over1m <- vect('bigdata/sw_mm_water_processed.gpkg')
# fourth, remove attributes and merge lines that are touching one another (takes a long time)
lines <- st_as_sf(wat.v.havewidth.perm.over1m)
lines <- lines[,!names(lines) %in% names(st_drop_geometry(lines))] # remove any attributes

nrow_start <- nrow(lines)
nrow_run <- nrow_start+1
ncol_start <- ncol(lines)

while (nrow_run > dim(lines)[1]) { # iterate through until complete
  
  nrow_run <- nrow(lines)
  
  lines$intersect_group <- unlist(map(st_intersects(lines),1)) # create grouping based on first line of intersection
  lines <- aggregate(lines, by=list(lines$intersect_group),FUN=first,do_union=TRUE)[,3] # aggregate data based on grouping
  
  print(paste0(nrow_start-nrow(lines)," features merged"))
  
}

# remove stretches of river less than 50m in length
lines.over50m <- lines %>% 
  mutate(length = as.numeric(st_length(.))) %>% 
  filter(length > 50)

# save output as 'bigdata/sw_mm_water_processed_merged.gpkg'
# st_write(lines.over50m,'bigdata/sw_mm_water_processed_merged.gpkg',delete_dsn = T)
rivers <- st_read('bigdata/sw_mm_water_processed_merged.gpkg')

# now 'punch out' known fords
# first, use highways fords
highways.fords <- st_read('bigdata/os_mm_highways_fords.gpkg')
highways.fords.buff <- highways.fords %>% 
  st_buffer(5) # small buffer required as point match rivers closely

highways.fords.buff.onrivers <- highways.fords.buff %>% 
  st_filter(lines.over50m)

highways.fords.onrivers <- highways.fords %>% 
  filter(gml_id %in% highways.fords.buff.onrivers$gml_id)

st_write(highways.fords.onrivers,'bigdata/os_mm_highways_fords_on_rivers.gpkg')
# just over 1000 fords on rivers

# convert lines to raster, and punch out fords later

# get raster template from os dem
dem <- rast('bigdata/os_50m_masked.tif')
rivers.ext <- ext(st_bbox(rivers))
rivers.ext.aligned <- align(rivers.ext,dem)

rivers.ext.aligned.c <- 
  c(rivers.ext.aligned[1],rivers.ext.aligned[3],
    rivers.ext.aligned[2],rivers.ext.aligned[4])

`%+%` <- function(x, y) paste0(x,y)

system(glue::glue(
  "gdal_rasterize -te {paste0(rivers.ext.aligned.c,collapse=' ')} " %+% 
    "-tr 50 50 -ot GTiff " %+% 
    "-burn NaN -init 1 -add -at {file.path(getwd(),'bigdata/sw_mm_water_processed_merged.gpkg')} " %+% 
    "{file.path(getwd(),'bigdata/sw_mm_water_processed_merged.tif')}"
)
)

# load raster
x <- rast('bigdata/sw_mm_water_processed_merged.tif')

# remove islands
x.sieved <- sieve(x, 300)
# writeRaster(x.sieved, 'bigdata/sw_mm_water_processed_merged_sieved.tif',overwrite=T)

# load fords on rivers and buffer by 75m to ensure at least 4 cells
fords <- buffer(vect('bigdata/os_mm_highways_fords_on_rivers.gpkg'),75)

ford.cells <- terra::extract(x.sieved,fords,cell=T)

x.sieved[ford.cells$cell] <- 1

# now add footpath fords - LINE DATA
path.fords <- st_read('bigdata/os_mm_paths_fords.gpkg')

path.fords.touchingrivers <- 
  path.fords %>% 
  st_filter(lines.over50m)

path.fords.touchingrivers.p <-
  path.fords.touchingrivers %>% 
  st_intersection(lines.over50m)

# load fords on rivers and buffer by 75m to ensure at least 4 cells
fords <- buffer(vect(path.fords.touchingrivers.p),75)

ford.cells <- terra::extract(x.sieved,fords,cell=T)

x.sieved[ford.cells$cell] <- 1

writeRaster(x.sieved,'bigdata/rivers_with_all_fords.tif',overwrite=T)

gbdem <- rast('bigdata/os_50m_masked.tif')

gbdem.crop <- crop(gbdem,x.sieved)

crs(x.sieved) <- crs(gbdem.crop)

gbdem.crop[x.sieved == 0] <- NaN

plot(gbdem.crop)

writeVector(wat.v.buff100.ag,'bigdata/sw_mm_water_processed_buffer100.gpkg')
plot(wat.v.buff100.ag)


wat.v.havewidth.perm.over1m.ag

l1 <- wat.v.havewidth.perm.over1m[wat.v.havewidth.perm.over1m$gml_id=='osgb5000005149300171',]
l2 <- wat.v.havewidth.perm.over1m[wat.v.havewidth.perm.over1m$gml_id=='osgb5000005146756014',]

l.all <- merge(l1,l2)


plot(wat.v.havewidth.perm.over1m[20:30,])

wat.v.havewidth.perm.over1m.ag <- aggregate(wat.v.havewidth.perm.over1m,dissolve=T)
plot(wat.v.havewidth.perm.over1m.ag)
terra::intersect()
wat.v.all.buff <- buffer(wat.v.all,100)
wat.v.all.buff.ag <- aggregate(wat.v.all.buff,dissolve=T)
wat.ext <- ext(wat.v.all)
wat.ext.aligned <- align(wat.ext,dem,snap='near')


plot(vect('bigdata/os_50m_polygon_outline.gpkg'),add=T)

wat.n <- vect(purrr::compact(wat.v))
wat.all <- rbind(wat.n,wat.tmp)
plot(wat.all)



g <- vect('bigdata/custom_grid/allpol.gpkg')
g$fid <- 1:nrow(g)
plot(g)
text(g, g$fid)
plot(wat.v.bind)

wat.all <- do.call('bind_rows',wat.c)
wat.all.v <- vect(wat.all)

wat.all.bbox <- wat.all %>% 
  st_bbox %>% st_as_sfc %>% vect()

gbdem <- rast('os_50m_masked.tif')
crsObj <- crs(gbdem)
originObj <- origin(gbdem)
extObj <- align(ext(wat.all.bbox),gbdem,snap='out')

r.clip <- rast(crs=crsObj,extent=extObj,resolution=c(50,50))
r.clip[] <- 1

wat.all.v$width
library(dplyr)
library(purrr)
library(terra)
paths <- 
  list.files('/media/mal/os_master_map_paths',full.names = T) %>% 
  map(vect)

future::plan('multisession',workers=6)

paths.fords <- paths %>% 
  map(.f = function(v) {
    v[v$formOfWay=='Path With Ford',]
  })

paths.fords.v <- vect(paths.fords)
writeVector(paths.fords.v,'bigdata/os_mm_paths_fords.gpkg')

highways <- vect('/media/mal/os_master_map_highways/Highways_Rrami_Hazard_FULL_001.gml.gz')
highways.fords <- highways[highways$hazard=='Ford',]
library(sf)
highways.fords.sf <- st_as_sf(highways.fords)
st_write(highways.fords.sf,'bigdata/os_mm_highways_fords.gpkg',delete_dsn = T)
writeVector(highways.fords, 'bigdata/os_mm_highways_fords.gpkg',overwrite=T)


lines <- st_read('bigdata/sw_mm_water_processed.gpkg',stringsAsFactors=FALSE)

# lines <- st_read(paste0(wd,input_file),stringsAsFactors=FALSE) # read in lines
lines <- lines[,!names(lines) %in% names(st_drop_geometry(lines))] # remove any attributes

nrow_start <- nrow(lines)
nrow_run <- nrow_start+1
ncol_start <- ncol(lines)

while (nrow_run > dim(lines)[1]) { # iterate through until complete
  
  nrow_run <- nrow(lines)
  
  lines$intersect_group <- unlist(map(st_intersects(lines),1)) # create grouping based on first line of intersection
  lines <- aggregate(lines, by=list(lines$intersect_group),FUN=first,do_union=TRUE)[,3] # aggregate data based on grouping
  
  print(paste0(nrow_start-nrow(lines)," features merged"))
  
}


