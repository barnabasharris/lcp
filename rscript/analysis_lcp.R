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
library(sp)
library(raster)
library(Matrix)
library(gdistance)
library(cppRouting)

# set up environment ----
if (str_detect(here(),'lustre')) {
  env <- 'KATHLEEN'
} else env <- 'LOCAL'

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

# define functions -----



# hex.w <- gbdem.hex.filtered.w[[3]]
# points.sf = gbdem.points

calcLCP <- 
  function(hex.w,
           points.sf,
           params) {
    
    library(sf)
    library(terra)
    library(purrr)
    library(glue)
    library(stringr)
    library(leastcostpath)
    library(Matrix)
    library(gdistance)
    library(cppRouting)
    library(vapour)
    library(data.table)
    
    # helper functions
    `%+%` <- function(x, y) paste0(x,y)
    
    
    source('rscript/lcp2.R')
    st <- Sys.time()
    # setwd(wd)
    print('compute node system tempdir is...')
    sysTmpDir <- Sys.getenv("TMPDIR")
    print(sysTmpDir)
    
    # nodeBatchNumber <- 1
    getLCPsLines <- function(nodeBatchNumber, gpkgOut) {
      
      idx <- nodesMaxList[[nodeBatchNumber]]
      
      print(glue::glue('calculating shortest paths for hex {hex.id}, ' %+% 
                         'node batch {nodeBatchNumber}...'))
      
      lcp_cells_list <- cppRouting::get_multi_paths(cs$cpp_cm_graph,
                                                    from=point.nodes.clean[idx],
                                                    to=point.nodes.clean,
                                                    long=F)
      print('done!')
      print('extracting routes and converting to line geometries...')
      # f <- lcp_cells_list[[1]]
      lcp.all <- lcp_cells_list %>% 
        map(.f = function(f) {

          # data.table / terra approach
          nodes.df = dplyr::bind_rows(
            lapply(f, function(x) data.frame(x)),
            .id = 'to'
          ) %>% 
            dplyr::rename('node' = x)
          
          # get coords
          coords <- 
            terra::xyFromCell(hex.r, as.integer(nodes.df$node)) %>% 
            as.data.frame()
                              
          nodes.df.coords <- 
            nodes.df %>% 
            dplyr::bind_cols(coords) %>% 
            dplyr::group_by(to) %>% 
            dplyr::mutate(object = dplyr::cur_group_id()) %>% 
            dplyr::ungroup() %>% 
            dplyr::mutate(part = object,
                          hole = 0) %>% 
            dplyr::select(c(object,part,x,y,hole)) %>% 
            as.matrix
          
          lcps <- terra::vect(x=nodes.df.coords,
                              crs=crs(hex.r)
                              ,type='lines'
          ) %>% 
            st_as_sf
          return(lcps)
        }) %>% do.call('bind_rows',.)
      print('done!')
      
      print('writing line geometries to geopackage...')
      if (!file.exists(gpkgOut)) {
        # writeVector(allLCPs, file.path(getwd(),'outputs',glue::glue("lcps_all.gpkg")),
        #             layer='lcps',overwrite=T)
        st_write(lcp.all, gpkgOut,
                 layer='lcps',delete_dsn= T,
                 quiet=T)
      } else {
        # writeVector(allLCPs, file.path(getwd(),'outputs',glue::glue("lcps_all.gpkg")),
        #             insert=T, layer='lcps')
        st_write(lcp.all, gpkgOut,
                 layer='lcps',append= T,
                 quiet=T)
      }
      # SF to export to single GPKG
      
      print(glue::glue('finished calculating shortest paths for hex {hex.id}, ' %+% 
                         'node batch {nodeBatchNumber}...'))
      return(gpkgOut)
    }
    
    # first extract the hex polygon
    hex.pol <- terra::unwrap(hex.w)
    hex.id <- hex.pol$hex_id
    
    sink(glue('logs/calc_lcp_sinkout_hex_{hex.id}.txt'))
    
    # subset points and convert to spatvector
    hex.points <- points.sf %>% 
      dplyr::filter(hex_id == hex.id) %>% 
      vect()
    
    # load dem and mask / trim
    print('trimming main raster to hex area...')
    gb.dem <- rast(params$demFile)
    hex.r <- terra::trim(terra::mask(gb.dem,hex.pol))
    
    # create cost surface
    print('creating cost surface...')
    cs <- 
      create_slope_cs2(hex.r,
                       cost_function = 
                         "campbell 2019", 
                       neighbours = params$lcpNeigh)
    # TODO: might be worth experimenting with {dodgr} to see if any speed gains
  
    # get cell numbers of points
    point.nodes <- terra::extract(cs$x, hex.points, cells=T)
    # check if any points fall on NaN cells, and remove if so
    point.nodes.clean <- point.nodes[!is.na(point.nodes[,2]),]$cell
    
    # split all terminus points into correct sized list
    maxPointsPerCore <- floor(params$maxPathsPerCore/length(point.nodes.clean))
    nodesIdx <- 1:length(point.nodes.clean)
    nodesMaxList <- split(nodesIdx, 
                          ceiling(seq_along(nodesIdx) / 
                                    maxPointsPerCore))
    
    # lines approach
    gpkgOut <- file.path(getwd(),'outputs',glue::glue("hex_{hex.id}.gpkg"))
    
    # process in chunks
    print('finding least cost paths for node batches...')
    1:length(nodesMaxList) %>% 
      walk(getLCPsLines,gpkgOut)
    
    print(glue::glue('finished calculating shortest paths for hex {hex.id}...'))
    
    tifOut <- str_replace(gpkgOut,'.gpkg','.tif')
    lcpsExt.v <- c(cs$extent[1],cs$extent[3],cs$extent[2],cs$extent[4])
    
    system(glue::glue(
      "gdal_rasterize -te {paste0(lcpsExt.v,collapse=' ')} " %+% 
        "-tr {paste0(params$demRes,collapse=' ')} -ot GTiff -l lcps " %+% 
        "-burn 1 -init 0 -add -at {gpkgOut} {tifOut}"
    )
    )
    
    file.remove(gpkgOut)
    sink()
    return(difftime(Sys.time(),st))
  }

# set parameters  -----
params <- list()
params$gridDist <- 1000 # dist between / density of path terminus points
params$cellsize <- 100000 # processing hex cell size, max dist 
params$edgeBuffer <- 2000 # distance to allow around edge of points, for paths that go beyond terminus
params$seed <- 1 # random seed
params$lcpNeigh <- 32 # lcp neighbours
params$maxPathsPerCore <- 300000 # limit of paths to process at one time (higher increases RAM)
params$demFile <- 'bigdata/os_50m_masked_w_navigable_rivers.tif'
set.seed(params$seed)

vectFile <- 'bigdata/os_50m_polygon_outline.gpkg'
ewVectFile <- 'bigdata/ew.gpkg'

# run script ----

##•┣ load / process data ----
# load terrain raster
gbdem <- rast(params$demFile)
# extract dem info
params$demCrs <- terra::crs(gbdem)
params$demExtent <- terra::ext(gbdem)
params$demRes <- terra::res(gbdem)
# load terrain vector outline
gbvect <- vect(vectFile)
# dissolve to remove any errors
gbvect.ag <- aggregate(gbvect)

##•┣ create hex grid ----
gbdem.hex <- gbvect.ag %>% 
  st_as_sf() %>% 
  st_make_grid(cellsize = params$cellsize,
               what = 'polygons',
               square = F,
               flat_topped = T) %>% 
  st_sf %>% 
  st_buffer(params$gridDist) %>%  # creates overlap between hex
  st_filter(st_as_sf(gbvect.ag),join=st_within) %>% 
  dplyr::mutate(hex_id = 1:nrow(.))

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
  st_buffer(params$edgeBuffer) %>% # ADD edge without points onto hex 
  vect()

# split into wrapped vectors for clusters
gbdem.hex.filtered.w <- 
  1:nrow(gbdem.hex.filtered) %>% 
  map(~wrap(gbdem.hex.filtered[.x,]))

# plot(gbdem)
# plot(gbdem.hex.filtered,add=T)
# text(gbdem.hex.filtered, gbdem.hex.filtered$hex_id, halo=FALSE, inside=FALSE)
lapply(gbdem.hex.filtered.w[1:4], function(x) {
  calcLCP(x, 
          points.sf = gbdem.points,
          params
  )
})

alarm()
##•┣ create cluster ----

print('getting MPI cluster...')
# cl <- snow::getMPIcluster() #
cl <- parallel::makePSOCKcluster(3) # testing
print('done!')

varsToExport <- 
  c("sysTmpDir")

# Display info about each process in the cluster
# print(clusterCall(cl, function() Sys.info()))
print('exporting vars to nodes...')

# snow::clusterExport(cl, varsToExport)
parallel::clusterExport(cl, varsToExport)
print(Sys.time())
print('running analysis...')



parallel::clusterApply(cl, 
                       gbdem.hex.filtered.w[1:3],
                       calcLCP,
                       points.sf = gbdem.points,
                       params)
                       

clusterApply(cl, 1:6, get("+"), 3)

hex2 <- 
  calcLCP(x = 2,
          hex.spt = gbdem.hex.filtered,
          points.sf = gbdem.points,
          dem = gbdem,
          params)

sessionTag <- 'east_wansdyke'

outTimes <- lapply(point_nodes[1:10], FUN=calcLCP, point_nodes, gbdem.ew, cpp_cm_graph,
                   landscapeTag=sessionTag)



plot(rast(glue::glue('outputs/{sessionTag}_lcps.tif')))
varsToExport <- c('cpp_cm_graph','nodes','cs_rast.w','params','wd')

lcpsExt <- ext(v)
lcpsExt.v <- c(lcpsExt$xmin,lcpsExt$ymin,lcpsExt$xmax,lcpsExt$ymax)
names(lcpsExt.v) <- NULL
system(glue::glue(
  "gdal_rasterize -te {paste0(lcpsExt.v,collapse=' ')} -tr 50 50 -ot GTiff -burn 1 -init 0 -add /media/mal/git/lcp/outputs/v.gpkg /media/mal/git/lcp/outputs/vt_new.tif"
)
)

r <- rast('outputs/vt_new.tif')
r.m <- mask(r,ewBuff.less.ag)
writeRaster(r.m,'5000m_buff_eastwans_1200p.tif')

# cl <- parallel::makePSOCKcluster(6)
# parallel::clusterExport(cl, 
#                         varlist = varsToExport)
# 
outTimes <- parallel::clusterApplyLB(cl, nodes[1:10], fun=calcLCP)

print('getting MPI cluster...')
cl <- snow::getMPIcluster()
print('done!')

# Display info about each process in the cluster
print(clusterCall(cl, function() Sys.info()))
print('exporting vars to nodes...')
snow::clusterExport(cl, varsToExport)
print(Sys.time())
print('running analysis...')
datOut <- snow::clusterApply(cl, nodes, fun=calcLCP)


# Clean up the cluster and release the relevant resources.
stopCluster(cl)
mpi.quit()




create cost surfaces and lcps
# x <- 5
# terminus generation
generateTerminus <- 
  function(gbdem.ew, ewvect, ewBuff, params,
           avoidEdgeBy=500, method='grid') {
    
    ewBuff.less <- buffer(ewvect, params$ewBuff-avoidEdgeBy)
    ewBuff.less.ag <- aggregate(ewBuff.less)
    
    # random band approach
    if (method == 'random_band') {
      ewBuff.diff <- terra::symdif(ewBuff.ag,ewBuff.less.ag)
      template.points.m <- terra::spatSample(ewBuff.diff, 1200, method='random')  
    }
    
    # perimeter line sample approach
    if (method == 'perimeter_line') {
      ewBuff.less.sf <- st_as_sf(ewBuff.less.ag) %>%
        st_cast('LINESTRING')
      p <- st_line_sample(ewBuff.less.sf, density=0.001) %>%
        st_sf %>%
        st_cast('POINT')
      template.points.m <- vect(p)
    }
    
    # grid approach
    if (method == 'grid') {
      # # # copy original r to template
      template <- terra::rast(gbdem.ew)
      # set new res of template
      res(template) <- params$gridDist
      # add arbitrary values
      values(template) <- 1:ncell(template)
      # convert to polygons
      template.points <- as.points(template)
      # clip grid by outline
      template.points.m <- terra::mask(template.points,ewBuff.ag)
    }
    return(template.points.m)
  }
