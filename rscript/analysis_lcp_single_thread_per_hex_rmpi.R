library(sf)
library(dplyr)
library(terra)
library(here)
library(purrr)
library(snow)
library(vapour)
library(glue)
library(stringr)
library(leastcostpath)
library(rgdal)
library(rgeos)

# library(sp)
# library(raster)
library(Matrix)
# library(gdistance)
library(cppRouting)
library(stringr)

# set up environment ----
wd <- here()

if (dir.exists('tmp')) {
  unlink('tmp',recursive=T)
}

dir.create('tmp')

# if (dir.exists('outputs')) {
#   unlink('outputs',recursive=T)
# }
# 
# dir.create('outputs')

if (dir.exists('logs')) {
  unlink('logs',recursive=T)
}

dir.create('logs')

# set system tmp
Sys.setenv(TMPDIR=file.path(here(),'tmp'))

print('head node system tempdir is...')
sysTmpDir <- Sys.getenv("TMPDIR")
print(sysTmpDir)

print('main R tempdir is...')
print(tempdir())

# source functions -----
source('rscript/lcp2.R')

# set parameters  -----
RcppParallel::setThreadOptions(numThreads = 1)
# RcppParallel::defaultNumThreads()
params <- list()
params$gridDist <- 1000 # dist between / density of path terminus points
params$cellsize <- 100000 # processing hex cell size, max dist 
params$edgeBuffer <- 2000 # distance to allow around edge of points
params$seed <- 1 # random seed
params$lcpNeigh <- 32 # lcp neighbours
params$demFile <- 'bigdata/os_50m_rivers_and_fords_sw.tif'
set.seed(params$seed)

vectFile <- 'bigdata/os_50m_polygon_outline.gpkg'
ewVectFile <- 'bigdata/ew.gpkg'

# run script ----

##•┣ load / process data ----
# load terrain raster
gbdem <- rast(params$demFile)
plot(gbdem)

# extract dem info
params$demCrs <- terra::crs(gbdem)
params$demExtent <- terra::ext(gbdem)
params$demRes <- terra::res(gbdem)

# load terrain vector outline
gbvect <- terra::crop(vect(vectFile),params$demExtent)

# dissolve to remove any errors
gbvect.ag <- terra::aggregate(gbvect)

##•┣ create hex grid ----

gbvect.ext <- ext(gbvect.ag)
gridOrigin <- c(gbvect.ext$xmin,gbvect.ext$ymin)
gridOrigin.offset <- c(gbvect.ext$xmin + (params$cellsize/1),
                       gbvect.ext$ymin + (params$cellsize/1)/2)
names(gridOrigin.offset) <- NULL

gbdem.hex <- gbvect.ag %>% 
  st_as_sf() %>% 
  st_make_grid(cellsize = params$cellsize,
               what = 'polygons',
               square = F,
               flat_topped = T) %>% 
  st_sf %>% 
  st_buffer(params$gridDist) %>%  # creates overlap between hex
  dplyr::mutate(hex_id = 1:nrow(.)) %>% 
  st_filter(st_as_sf(gbvect.ag),join=st_within)

gbdem.hex.offset <- gbvect.ag %>% 
  st_as_sf() %>% 
  st_make_grid(cellsize = params$cellsize,
               offset = gridOrigin.offset,
               what = 'polygons',
               square = F,
               flat_topped = T) %>% 
  st_sf %>% 
  st_buffer(params$gridDist) %>%  # creates overlap between hex
  dplyr::mutate(hex_id = 1:nrow(.)) %>% 
  dplyr::filter(hex_id %in% gbdem.hex$hex_id) %>% 
  st_filter(st_as_sf(gbvect.ag),join=st_touches)

gbdem.points <- gbvect.ag %>% 
  st_as_sf() %>% 
  st_make_grid(cellsize = params$gridDist,
               what = 'centers') %>% 
  st_sf %>% 
  st_join(gbdem.hex) %>% 
  st_filter(st_as_sf(gbvect.ag))

whichHex <- gbdem.points %>% 
    st_drop_geometry() %>% 
    group_by(hex_id) %>% 
    summarise(hex_num = n()) %>% 
  filter(hex_num > 20)
      
gbdem.hex.filtered <- 
  gbdem.hex %>% 
  dplyr::filter(hex_id %in% whichHex$hex_id) %>% 
  st_buffer(params$edgeBuffer) %>% # ADD edge without points
  vect()

# split into wrapped vectors for clusters
gbdem.hex.filtered.w <- 
  1:nrow(gbdem.hex.filtered) %>% 
  map(~wrap(gbdem.hex.filtered[.x,]))


lapply(gbdem.hex.filtered.w,
       calcLCP,
       points.sf=gbdem.points,
       params,
       maxPathsPerCore=100000,
       justCS=F
       )
##•┣ create cluster ----

# kathleen / non-interactive
# print('getting MPI cluster...')
# cl <- snow::getMPIcluster()
# print('done!')

# myriad / interactive
cl <- snow::makeSOCKcluster(13)


# Display info about each process in the cluster
print(clusterCall(cl, function() Sys.info()))

##•┣ fire! ----
print(Sys.time())
print('running analysis...')

snow::clusterEvalQ(cl,source(file.path(getwd(),'rscript','lcp2.R')))

l <- snow::clusterApplyLB(cl,
                     gbdem.hex.filtered.w[1:13],
                     calcLCP,
                     points.sf=gbdem.points,
                     params,
                     maxPathsPerCore=200000,
                     refreshCS=F)

gbdem.hex.filtered.w %>% 
  map(unwrap)

source('rscript/lcp2.R')

calcLCP(gbdem.hex.filtered.w[[10]],
                          points.sf=gbdem.points,
                          params,
                          maxPathsPerCore=200000,
                          refreshCS=F)


rams
datOut <- 
  snow::clusterApplyLB(cl,
                       gbdem.hex.filtered.w[19:22],
                       calcLCP,
                       points.sf=gbdem.points,
                       params,
                       maxPathsPerCore=300000)

# Clean up the cluster and release the relevant resources.
stopCluster(cl)
mpi.quit()

getwd()
r <- rast('/media/mal/git/lcp/outputs/hex_4_tmp.tif')
r1 <- rast('outputs/hex_4.tif')
r2 <- r + r1
r2 + r
plot(r2)
system('gdalwarp -r sum /media/mal/git/lcp/outputs/hex_4_tmp.tif /media/mal/git/lcp/outputs/hex_4.tif /media/mal/git/lcp/outputs/hex_4_warp.tif')

r3 <- rast('outputs/hex_4_warp.tif')
plot(r3)
