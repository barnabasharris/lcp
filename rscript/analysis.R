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
library(sp)
library(raster)
library(Matrix)
library(gdistance)


if (str_detect(here(),'lustre')) {
  env <- 'KATHLEEN'
} else env <- 'LOCAL'

wd <- here()

if (dir.exists('tmp')) {
  unlink('tmp',recursive=T)
}

dir.create('tmp')


if (dir.exists('outputs')) {
  unlink('outputs',recursive=T)
}

dir.create('outputs')

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

# set permissions for python script to allow executable
# system(glue('chmod -R +x {wd}/python/'))
# system(glue('chmod -R +x {wd}/bash/'))

# pre-process date for analysis -----
# gridRes <- 6500 # decided by maximum job length allowed by kathleen, see  below
params <- list()
params$gridDist <- 10000

demFile <- 'bigdata/os_50m_masked.tif'
vectFile <- 'bigdata/os_50m_polygon_outline.gpkg'
# load terrain raster
gbdem <- rast(demFile)
v <- vect(vectFile)

# copy original r to template
template <- terra::rast(gbdem)
# set new res of template
res(template) <- params$gridDist
# add arbitrary values
values(template) <- 1:ncell(template)
# convert to polygons
template.points <- as.points(template)
# clip grid by outline
template.points.m <- terra::mask(template.points,v)

# check coverage
# pols <- 1:1000 %>% 
#   map(.f = function(x) {
#     d <- symdif(buffer(template.points.m[x,],48000),
#            buffer(template.points.m[x,],32000))
#   }) %>% do.call('rbind',.)
# plot(pols, alpha=0.1,col='red',border=NA)

params$pathRange <- c(32000,48000) # a day's walk (should be larger than gridDist)
params$maxPathDist <- 48000 # a day's walk
params$numPaths <- 100 # number of paths from central point
params$edgeBuffer <- 5000
params$seed <- 1

set.seed(params$seed)
# terra::spatSample very slow
v.samp <- sf::st_sample(sf::st_as_sf(v),10000,type='regular') %>% 
  terra::vect()
# select mid point in grid
mid <- length(v.samp)/2
# find surrounding points
n <- as.numeric(nearby(v.samp[mid,],v.samp,k=9)) %>% 
  .[3:length(.)]
# calc distance from centre to surrounding point(s)
dists <- distance(v.samp[mid,],v.samp[n,])
meanDist <- mean(dists)
# approximate number of neighbours to achieve max dist
params$maxPathDist




centre <- vect(data.frame(lon=442098,lat=379548),crs=crs(gbdem)) %>% 
  buffer(50000)

d <- distance(mask(v.samp,centre))
hist(d)


plot(gbdem)
sf::st_distance(v.samp[1:1000,])


# loop through each central point
# x <- 200
template.points.m.wrapped <- wrap(template.points.m)
demres <- res(gbdem)
demcrs <- crs(gbdem)

cl <- parallel::makePSOCKcluster(6)
parallel::clusterExport(cl, 
                        varlist = c('template.points.m.wrapped',
                                    'demres',
                                    'demcrs',
                                    'demFile',
                                    'params'))


# x <- 199
lcp.tiles <-
  parallel::clusterApplyLB(cl, 1:200, function(x) {
    library(terra)
    print(x)
    # unwrap
    template.points.m <- vect(template.points.m.wrapped)
    p <- template.points.m[x,]
    
    # total width is 2 * (max pathRange + edgeBuffer)
    p.buff.wedge <- buffer(p,max(params$pathRange) + params$edgeBuffer)
    dm <- (max(params$pathRange) + params$edgeBuffer) / demres
    crsObj <- demcrs
    exSpat <- ext(p.buff.wedge)
    exNum <- c(exSpat$xmin,exSpat$xmax,exSpat$ymin,exSpat$ymax)
    names(exNum) <- NULL
    
    vals <- vapour::vapour_warp_raster(demFile, extent = exNum, dimension = dm*2, projection = crsObj)
    # can't set CRS here, as for some reason Kathleen throws an error when unwrapping the rast
    r <- setValues(rast(extent=exSpat, nrows = dm[2]*2, ncols = dm[1]*2), vals[[1]])
    
    p.buff <- buffer(p,max(params$pathRange))
    p.minbuff <- buffer(p, min(params$pathRange))
    pathzone <- symdif(p.buff,p.minbuff)
    r.pathzone <- mask(r, pathzone)
    
    # check if all not nan
    if (length(r.pathzone[is.na(r.pathzone)]) == terra::ncell(r.pathzone)) {
      print('no landscape cells in path zone...')
      return(NULL)
    }
    
    # distribute points 
    destinations <- spatSample(pathzone, params$numPaths)
    
    # extract vals from raster
    destinations.r <- terra::extract(r.pathzone, destinations)
    # keep only non-NaN points
    idx <- which(!is.nan(destinations.r[,2]))
    destinations.ok <- destinations[idx,]
    
    allp <- buffer(rbind(destinations.ok,p),params$edgeBuffer)
    processingExtent <- ext(allp)
    processingExtent.pol <- as.polygons(processingExtent)
    terra::crs(processingExtent.pol) <- demcrs
    r.crop <- crop(r,buffer(processingExtent.pol,params$edgeBuffer))
    plot(r)
    plot(destinations.ok,add=T)
    plot(p, add=T, col='red')
    
    r.dem <- wrap(r.crop)
    p$type <- 'origin'
    destinations.ok$type <- 'destination'
    r.points <- wrap(rbind(p,destinations.ok))
    
    
    return(list(r.dem = r.dem,
                r.points = r.points))
  })

paramString <- params %>% map2(.y = names(params), ~paste0(.y, '_',.x,collapse='_')) %>% 
  paste0(collapse = '_')

saveRDS(lcp.tiles, paste0('bigdata/lcp_tiles_test_',paramString,'.RDS'))

lcp.tile <- lcp.tiles[[150]]
lcp <- function(lcp.tile, neigh = 32) {
  library(terra)
  library(sf)
  r <- rast(lcp.tile$r.dem)
  tile.points <- vect(lcp.tile$r.points)
  crs(r) <- crs(tile.points)
  
  tile.points <- 
    st_as_sf(tile.points)
  
  slope_cs <- 
    create_slope_cs(r,cost_function = "tobler", neighbours = neigh)
  
  # rivers_cs <- 
  #   create_barrier_cs(raster = .rivers.r.tile,
  #                     barrier = .rivers.r.tile,
  #                     neighbours = neigh,
  #                     field=0)
  # 
  # slope_rivers_cs <- slope_cs * rivers_cs
  st <- Sys.time()
  lcp <- create_lcp(slope_cs, 
                    origin = tile.points %>% dplyr::filter(type == 'origin'), 
                    destination = tile.points %>% dplyr::filter(type == 'destination') %>% 
                      dplyr::slice(n = 50))
  Sys.time() - st
  
  plot(lcp$geometry,add=T)
}

((((length(template.points.m) * 100) * 34) / 60) / 60) / 250

# get rivers
rivers <- st_read(con, Id(schema='ext',table='open_rivers_processed_simplified_50m'))

t.p <- as.polygons(ext(r)) %>% 
  st_as_sf

rivers.tile.buff <- 
  st_crop(rivers,t.p) %>% 
  st_buffer(100) %>% 
  as('SpatVector')

rivers.r.tile <- 
  terra::rasterize(rivers.tile.buff, rast(r))

# rivers.r.all <- rast('bigdata/open_rivers_50m.tif')
# rivers.r.tile <- terra::crop(rivers.r.all,r)
rivers.r.tile[rivers.r.tile==1] <- 0
rivers.r.tile[rivers.r.tile!=0] <- NA
.crs <- raster::crs(terra::crs(rivers.r.tile))

.r <- as(r,'Raster')
raster::crs(.r) <- .crs

.rivers.r.tile <- as(rivers.r.tile,'Raster')
raster::crs(.rivers.r.tile) <- .crs

slope_cs <- create_slope_cs(dem = .r, cost_function = "tobler", neighbours = neigh)

rivers_cs <- 
  create_barrier_cs(raster = .rivers.r.tile,
                    barrier = .rivers.r.tile,
                    neighbours = neigh,
                    field=0)

slope_rivers_cs <- slope_cs * rivers_cs

dem_extent <- as(raster::extent(.r), 'SpatialPolygons')

lcp_A <- sp::SpatialPoints(cbind(233722,764293))
lcp_B <- sp::SpatialPoints(cbind(245253,764438))

plot(raster(slope_rivers_cs))
plot(lcp_A,add=T)
plot(lcp_B,add=T)

lcp <- create_lcp(cost_surface = slope_rivers_cs, 
           origin = lcp_A, 
           destination = lcp_B, 
           directional = FALSE)

plot(lcp,add=T,col='black')
# if parallel is TRUE then need to specify number of cores via ncores argument 

fete_lcp <- create_FETE_lcps(cost_surface = slope_rivers_cs, 
                             locations = rbind(lcp_A,lcp_B), 
                             cost_distance = FALSE, 
                             parallel = FALSE)
fete_lcp.v <- as(fete_lcp,'SpatVector')
writeVector(fete_lcp.v,'lcp.gpkg')
raster::writeRaster(raster(slope_rivers_cs),'lcp_cs.tif')


plot(fete_lcp,add=T)
# # visualize
# if (env == 'LOCAL') {
#   plot(r)
#   plot(template.pols.m, add=T)
# }
# #



r.tiles.length <- r.tiles %>%
  map_int(.f = function(x) {
    length(as.points(rast(x$r.views)))
  })

r.tiles <- readRDS('bigdata/rtiles_25000m.RDS')
# multi-thread analysis ----

##•┣ define function ----
# x <- 1
# x <- 17
# tile <- r.tiles[[17]]


viewpointAnalysis <- function(tile) {
  st <- Sys.time()
  library(terra)
  library(glue)
  setwd(wd)
  
  # extract tile number
  x <- tile$r.num
  
  sink(glue('logs/analysis_sinkout_{x}.txt'))
  
  print('compute node system tempdir is...')
  sysTmpDir <- Sys.getenv("TMPDIR")
  print(sysTmpDir)
  
  # extract tile
  print('extracting raster tile...')
  r <- rast(tile$r.views)
  crs(r) <- 'EPSG:27700'
  
  # make into points
  print('converting to points...')
  r.points <- as.points(r)
  r.points.length <- length(r.points)
  
  # write points to disk
  print('writing points to disk...')
  pointsLoc <- glue('{sysTmpDir}/viewpoints_{x}.gpkg')
  if (file.exists(pointsLoc)) file.remove(pointsLoc)
  writeVector(r.points,pointsLoc,overwrite=T)
  
  # write section of dem to disk
  r <- rast(tile$r.dem)
  crs(r) <- 'EPSG:27700'
  demLoc <- glue('{sysTmpDir}/viewdem_{x}.tif')
  if (file.exists(demLoc)) file.remove(demLoc)
  writeRaster(r,demLoc,overwrite=T)
  
  viewobserver <-  1.75
  viewtarget <- 0
  pointsloc <- pointsLoc
  output_loc <- file.path(wd,'outputs')
  output_it <- x
  
  print('node preparation complete')
  diff <- Sys.time() -  st
  print(diff)
  # execute bash script
  st <- Sys.time()
  
  visOutputLocs <- glue('{sysTmpDir}/vislocs_{x}.gpkg')
  visOutputCum <- glue('outputs/cum_view_{x}.tif')
  # call viewshed plugin through bash script
  system2(glue("{getwd()}/bash/visibility_viewpoint_locations.sh"),
          glue('{demLoc} {pointsLoc} {viewdist} {viewobserver} {viewtarget} {visOutputLocs}'))
  
  system2(glue("{getwd()}/bash/visibility_viewshed.sh"),
          glue('{demLoc} {visOutputLocs} {visOutputCum}'),
          stderr = paste0(getwd(),'/logs/viz_e',output_it,'.txt'),
          stdout = paste0(getwd(),'/logs/viz_o',output_it,'.txt'))
  
  print(glue('viewshed analysis of {r.points.length} points complete'))
  diff <- Sys.time() -  st
  print(diff)
  
  print(glue('Time per point is: {round(difftime(Sys.time(),st,units="secs") / r.points.length,5)} seconds'))
  
  # remove the grass mapset
  # unlink(grassMapset,recursive = T)
  sink()
  return(print(glue('node with job {x} finished')))
}

varsToExport <- c('sysTmpDir','wd','viewdist')

# remove already processed tiles
alreadyProcessed <- list.files('outputs') %>% 
  map(~stringr::str_extract(.x,'[0-9]+')) %>% 
  unlist() %>% 
  as.numeric()

toProcess <- setdiff(1:length(r.tiles),alreadyProcessed)


##•┣ {parallel} version for testing ----
if (env == 'LOCAL') {
  cl <- parallel::makePSOCKcluster(6)
  parallel::clusterExport(cl, varlist = varsToExport)
  datOut <- parallel::clusterApplyLB(cl, r.tiles[toProcess], viewpointAnalysis)
}


##•┣ {snow} version for kathleen ----
if (env == 'KATHLEEN') {
  #!! if offSet = T the use pd$tiles$pol !!
  print('getting MPI cluster...')
  cl <- snow::getMPIcluster()
  print('done!')
  
  # Display info about each process in the cluster
  print(clusterCall(cl, function() Sys.info()))
  print('exporting vars to nodes...')
  snow::clusterExport(cl, varsToExport)
  print(Sys.time())
  print('running analysis...')
  datOut <- snow::clusterApply(cl, c(17L, 25L, 84L, 96L, 97L, 99L, 108L, 115L, 116L, 117L, 118L, 
                                     119L, 120L, 121L, 122L, 123L, 128L, 133L, 136L, 138L, 139L, 140L, 
                                     141L, 142L, 143L, 144L, 145L, 146L, 147L, 148L, 155L, 156L, 157L, 
                                     163L, 164L, 165L, 166L, 167L, 168L, 169L, 170L, 171L, 172L, 173L, 
                                     174L, 175L, 176L, 177L, 178L, 184L, 185L, 186L, 192L, 193L, 194L, 
                                     195L, 196L, 197L, 198L, 199L, 200L, 201L, 202L, 203L, 204L, 205L, 
                                     206L, 207L, 214L, 215L, 226L, 227L, 228L, 229L, 230L, 231L, 232L, 
                                     233L, 234L, 235L, 236L, 237L, 242L, 244L, 245L, 251L, 252L, 255L, 
                                     256L, 257L, 258L, 259L, 260L, 261L, 262L, 263L, 264L, 265L, 266L, 
                                     271L, 279L, 280L, 281L, 282L, 283L, 284L, 285L, 286L, 287L, 288L, 
                                     289L, 290L, 291L, 298L, 305L, 306L, 307L, 308L, 309L, 310L, 311L, 
                                     312L, 313L, 314L, 315L, 316L, 317L, 330L, 331L, 332L, 333L, 334L, 
                                     335L, 336L, 337L, 338L, 339L, 340L, 341L, 358L, 359L, 360L, 361L, 
                                     362L, 363L, 364L, 365L, 366L, 386L, 387L, 388L, 389L, 390L, 393L, 
                                     409L, 410L, 411L, 412L, 413L, 414L),
                               viewpointAnalysis, cva=F)
  # datOut <- snow::clusterApplyLB(cl, 1:50, viewpointAnalysis)
}

# remove grass location
unlink(grassloc,recursive = T)
print(Sys.time())
print('done!')
# Clean up the cluster and release the relevant resources.
stopCluster(cl)
mpi.quit()




