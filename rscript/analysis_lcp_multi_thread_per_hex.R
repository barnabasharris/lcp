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

if (str_detect(sysNodeName,'lustre')) {
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



# hex.w <- gbdem.hex.filtered.w[[15]]
# points.sf = gbdem.points
# maxPathsPerCore=300000
calcLCP <- 
  function(hex.w,
           cl,
           points.sf,
           params,
           maxPathsPerCore=300000) {
    
    # helper functions
    `%+%` <- function(x, y) paste0(x,y)
    `%!in%` <- function(x,y)!('%in%'(x,y))
    
    source('rscript/lcp2.R')
    
    st <- Sys.time()
    # setwd(wd)
    print('compute node system tempdir is...')
    sysTmpDir <- Sys.getenv("TMPDIR")
    print(sysTmpDir)
    
    ##•┣ define func ----
    
    # idx <- 1:10
    # idx <- nodesMaxList[[2]]
    # cpp_cm_graph <- cs$cpp_cm_graph
    getLCPsLines <- function(idx, hex.id, hex.r.w,
                             cpp_cm_graph,
                             point.nodes.clean,
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
      `%+%` <- function(x, y) paste0(x,y)
      `%!in%` <- function(x,y)!('%in%'(x,y))
      
      nodeBatchNumber <- paste0(idx[1],'-',idx[length(idx)])
      
      sink(glue('logs/calc_lcp_sinkout_hex_{hex.id}_node_batch_{nodeBatchNumber}.txt'))
      
      # unwrap rast
      hex.r <- rast(hex.r.w)
      
      # lines approach
      gpkgOut <- file.path(getwd(),'outputs',glue::glue("hex_{hex.id}_batch_{nodeBatchNumber}.gpkg"))
      
      print(glue::glue('calculating shortest paths for hex {hex.id}, ' %+% 
                         'node batch {nodeBatchNumber}...'))
      
      lcp_cells_list <- cppRouting::get_multi_paths(cpp_cm_graph,
                                                    from=point.nodes.clean[idx],
                                                    to=point.nodes.clean,
                                                    long=F)
      
      print('done!')
      print('extracting routes and converting to line geometries...')
      # f <- lcp_cells_list[[50]]
      # names(lcp_cells_list)[50]

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
          )
          return(lcps)
        })
      print('done!')
      
      print('selecting only points with at least one successful path...')
      lcp.all.c <- lcp.all[unlist(map(lcp.all, terra::geomtype))=='lines']
    
      print('binding paths into a single vect object..')
      lcp.all.v <- vect(lcp.all.c)
      
      print('writing line geometries to geopackage...')
      writeVector(lcp.all.v, gpkgOut, overwrite=T)
      print('done!')
      
      # gpkgOut <- '/media/mal/git/lcp/outputs/hex_4_batch_65-128.gpkg'
      tifOut <- str_replace(gpkgOut,'.gpkg','.tif')
      
      hex.exent <- ext(hex.r)
      
      lcpsExt.v <- c(hex.exent[1],hex.exent[3],hex.exent[2],hex.exent[4])
      
      system(glue::glue(
        "gdal_rasterize -te {paste0(lcpsExt.v,collapse=' ')} " %+% 
          "-tr {paste0(params$demRes,collapse=' ')} -ot GTiff " %+% 
          "-burn 1 -init 0 -add -at {gpkgOut} {tifOut}"
      )
      )
      
      file.remove(gpkgOut)
      
      
      print(glue::glue('finished calculating shortest paths for hex {hex.id}, ' %+% 
                         'node batch {nodeBatchNumber}...'))
      sink()
      return(gpkgOut)
    }
    
    ##•┣ run script ----
    # first extract the hex polygon
    hex.pol <- terra::vect(hex.w)
    hex.id <- hex.pol$hex_id
    
    # subset points and convert to spatvector
    hex.points <- points.sf %>% 
      dplyr::filter(hex_id == hex.id) %>% 
      vect()
    
    # load dem and mask / trim
    print(glue::glue('trimming main raster to hex {hex.id} area...'))
    gb.dem <- rast(params$demFile)
    hex.r <- terra::trim(terra::mask(gb.dem,hex.pol))
    # wrap it up
    hex.r.w <- wrap(hex.r)
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
    maxPointsPerCore <- floor(maxPathsPerCore/length(point.nodes.clean))
    nodesIdx <- 1:length(point.nodes.clean)
    # nodesMaxList <- split(nodesIdx, sort(nodesIdx%%maxPointsPerCore))
    nodesMaxList <- split(nodesIdx, 
                          ceiling(seq_along(nodesIdx) / 
                                    maxPointsPerCore))
    
    # process in chunks
    
    print(glue::glue('finding least cost paths for {length(nodesMaxList)} node batches...'))
    # peakRAM::peakRAM({
    # getLCPsLines(nodesMaxList[[2]],
    #              hex.id,
    #              hex.r.w,
    #              cs$cpp_cm_graph,
    #              point.nodes.clean,
    #              params)
    # })
    
    # lapply(nodesMaxList,
    #       FUN=getLCPsLines,
    #      hex.id, 
    #      hex.r.w,
    #      cs$cpp_cm_graph,
    #      point.nodes.clean,
    #      params)
                             
    ##•┣ call func ----
    
    # check for pre-existing hex and node batches
    hexBatchOutputs <- nodesMaxList %>% 
      map_chr(.f = function(idx) {
        nodeBatchNumber <- paste0(idx[1],'-',idx[length(idx)])
        tifOut <- file.path(getwd(),'outputs',glue::glue("hex_{hex.id}_batch_{nodeBatchNumber}.tif"))
      })
  
    toProcess <- 
      which(hexBatchOutputs %!in% 
              list.files(file.path(getwd(),'outputs'),full.names = T))
    
    snow::clusterApplyLB(cl,
                       nodesMaxList[toProcess],
                       getLCPsLines,
                       hex.id, 
                       hex.r.w,
                       cs$cpp_cm_graph,
                       point.nodes.clean,
                       params)
    
    print(glue::glue('finished calculating shortest paths for hex {hex.id}...'))
    
    # sum rasters
    hexRasOut <- glue::glue('outputs/hex_{hex.id}.tif')
    hexRasOutTmp <- glue::glue('outputs/hex_{hex.id}_tmp.tif')
    hex.r.zero <- hex.r
    hex.r.zero[!is.nan(hex.r.zero)] <- 0
    writeRaster(hex.r.zero,hexRasOut,overwrite=T)
    
    lcpRas <- file.path(getwd(),list.files('outputs',full.names = T,
               pattern = glue::glue('^hex_{hex.id}_batch.*tif'))) %>% 
      map(.f = function(r) {
        system(glue::glue('gdal_calc.py -A {hexRasOut} -B {r} --outfile={hexRasOutTmp} --calc="A+B"'))
        file.copy(hexRasOutTmp, hexRasOut, overwrite=T)
      })
    
    file.remove(hexRasOutTmp)
    file.remove(file.path(getwd(),list.files('outputs',full.names = T,
                                             pattern = glue::glue('^hex_{hex.id}_batch.*tif'))))
    
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
  st_buffer(params$edgeBuffer) %>% # ADD edge without points onto hex 
  vect()

# split into wrapped vectors for clusters
gbdem.hex.filtered.w <- 
  1:nrow(gbdem.hex.filtered) %>% 
  map(~wrap(gbdem.hex.filtered[.x,]))

##•┣ create cluster ----

# kathleen / non-interactive
print('getting MPI cluster...')
cl <- snow::getMPIcluster()
print('done!')

# myriad / interactive
cl <- snow::makeSOCKcluster(30)

# Display info about each process in the cluster
print(clusterCall(cl, function() Sys.info()))

##•┣ fire! ----
print(Sys.time())
print('running analysis...')
datOut <- lapply(gbdem.hex.filtered.w[19:22],
                 FUN=calcLCP, cl, 
                 points.sf=gbdem.points,
                 params,
                 maxPathsPerCore=300000)

# Clean up the cluster and release the relevant resources.
stopCluster(cl)
mpi.quit()


# create cost surfaces and lcps
# x <- 5


hex4 <- rast('outputs/hex_4.tif')
hex4.f10 <- focal(hex4,9,fun='mean')
plot(hex4.f10)


os50m.v <- st_read('bigdata/os_50m_polygon_outline.gpkg')

polgrid <- st_make_grid(os50m.v, cellsize = 5000) %>% 
  vect()

polgrid.osgb <- mask(polgrid,vect(os50m.v))
polgrid.osgb.ag <- aggregate(polgrid.osgb,dissolve=T)
writeVector(polgrid.osgb.ag,'bigdata/gbgrid.gpkg')

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

oswater <- list.files('~/Downloads',pattern='21479*',full.names = T)
file.copy(oswater, paste0('/media/mal/os_master_map_water_network/',basename(oswater)))

list.files('/media/mal/os_master_map_water_network/')
unzip(paste0('/media/mal/os_master_map_water_network/',basename(oswater))[1])

expanse(poldiv.sw[5,]) / (5000^2)


plot(poldiv.sw[1,])
plot(mask(polgrid,poldiv.sw[1,]),add=T)

poldiv.sw.1sf <- poldiv.sw[1,] %>% 
  st_as_sf() 

poldiv.sw.1sf.grid <- poldiv.sw.1sf %>% 
  st_make_grid(cellsize = 5000)
expanse()
plot(poldiv.sw.1sf.grid[poldiv.sw.1sf])






