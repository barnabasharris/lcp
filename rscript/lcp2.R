

cf <- function(x){
  (terms$term_c*(1/(pi*terms$term_b*(1+((slope2deg(x)-terms$term_a)/terms$term_b)^2)))+terms$term_d)+terms$term_e*slope2deg(x)
}

neighbourhood <- function(neighbours) { 
  
  neighbours_32 <- matrix(c(0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 
                            1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0), nrow = 7, ncol = 7, byrow = TRUE)
  
  neighbours_48 <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 
                            1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 
                            1, 0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 9, ncol = 9, byrow = TRUE)
  
  if (inherits(neighbours, "matrix")) { 
    neighbours <- neighbours 
  } else if (neighbours == 4) {  
    neighbours <- 4
  } else if (neighbours == 8) {
    neighbours <- 8      
  } else if (neighbours == 16) { 
    neighbours <- 16
  } else if (neighbours == 32) { 
    neighbours <- neighbours_32
  } else if (neighbours == 48) { 
    neighbours <- neighbours_48
  } else (
    stop(paste0("neighbours argument invalid. Expecting 4, 8, 16, 32, 48, or a matrix object"))
  )
  
  return(neighbours)
}



# list.files('leastcostpath/R',full.names = T) %>% 
#   map(source)
# x = gbdem.ew
# cost_function = "campbell 2019"
# crit_slope = 12
# max_slope = NULL
# percentile = 0.5
# exaggeration = FALSE
# neighbours <- 32
# neighbours = params$lcpNeigh
# x <- rast('test_dem.tif')
# x <- rast('r.asc')
create_slope_cs2 <- function(x, 
                             cost_function = "tobler", 
                             neighbours = 16, crit_slope = 12, max_slope = NULL, percentile = 0.5, exaggeration = FALSE) {
  library(leastcostpath)
  
  `%!in%` <- function(x,y)!('%in%'(x,y))
  
  list.files('leastcostpath/R',full.names = T) %>%
    map(source)
  
  if(terra::is.lonlat(x)) { 
    stop("supplied digital elevation model (DEM) is invalid. x argument expects DEM with a projected coordinate system")
  }
  # neighbours <- 32
  neighbours <- neighbourhood(neighbours = neighbours)
  
  # duplicate matrix but with all cells positive
  # 
  # message('removing low lying land and resulting islands...')
  # # remove low-lying land
  # x[x <= 0] <- NaN
  # # remove islands 8
  # freqNan <- function(x) sum(is.nan(x))
  # x.foc <- terra::focal(x,w=3,freqNan,na.policy='omit')
  # 
  # # count nan
  # numNanIsland <- 7
  # islandCells <- sum(terra::values(x.foc) >= numNanIsland,na.rm=T)
  # 
  # while (islandCells > 0) {
  #   x.foc[x.foc >= numNanIsland] <- NaN
  #   x.foc.tmp <- terra::focal(x.foc,w=3,freqNan,na.policy='omit')
  #   islandCells <- sum(terra::values(x.foc.tmp) >= numNanIsland,na.rm=T)
  #   x.foc <- x.foc.tmp
  # }
  # # remove cells identified as islands
  # x[is.nan(x.foc)] <- NaN
  
  # generate neighbours matrix using all neighbours and then subset to ensure
  # cells are consistent across comparisons
  
  x.values <- terra::values(x)
  cells <- which(!is.na(x.values))
  na_cells <- which(is.na(x.values))
  
  adj <- terra::adjacent(x = x, cells = cells, directions = 
                           neighbours, pairs = TRUE)
  adj <- adj[!adj[,2] %in% na_cells,]
  
  message('detecting and removing links crossing NaNs...')
  # identify danger cells - those with 1 more more nan neighbours
  freqNan <- function(x) sum(is.nan(x))
  x.foc <- terra::focal(x,w=3,freqNan,na.policy='omit')
  
  # get danger cell numbers
  dangerCells.all <- which(terra::values(x.foc) >= 1) 
  
  # option 1: remove danger cell from / to links from adj matrix
  {
    adj <- adj[adj[,1] %!in% dangerCells.all,]
    adj <- adj[adj[,2] %!in% dangerCells.all,]

    # limit these links to just 8 neighbours
    adj.danger <- terra::adjacent(x = x, cells = dangerCells.all,
                                  directions = 'queen',
                                  pairs = TRUE)
    adj.danger <- adj.danger[!adj.danger[,2] %in% na_cells,]
    adj <- rbind(adj,adj.danger)
  }
  
  # option 2 - visualise linkages using polylines and remove
  # based on extracting raster values
  {
    # adj.danger <- adj[adj[,1] %in% dangerCells.all,]
    # 
    # fromXY <- xyFromCell(x, adj.danger[,1]) %>%
    #   cbind(.,object = 1:nrow(.))
    # 
    # toXY <- xyFromCell(x, adj.danger[,2]) %>%
    #   cbind(.,object = 1:nrow(.))
    # 
    # pointMat <- rbind(fromXY, toXY) %>%
    #   cbind(.,part = 1) %>%
    #   .[,c(3,4,1,2)] %>%
    #   .[order(.[,1]),]
    # # create vector lines showing path of links
    # v <- vect(pointMat,'lines',crs="EPSG:27700")
    # # add to and from cell numbers
    # v$from <- adj.danger[,1]
    # v$to <- adj.danger[,2]
    # # extract dem values along links lines
    # ras.ex <- extract(x,v)
    # 
    # # count number of NAs
    # ras.ex$nan <- as.integer(as.logical(is.na(ras.ex[,2])))
    # ras.ex.sum <- rowsum(ras.ex[,3], ras.ex[,1],na.rm=T)
    # ras.ex.sum.id <- as.data.frame(cbind(ras.ex.sum,unique(ras.ex[,1])))
    # names(ras.ex.sum.id) <- c('nan','ID')
    # dangerLinks <- ras.ex.sum.id[ras.ex.sum.id$nan >= 3,]
    # 
    # # visualize
    # r.all <- x
    # r.all[] <- 1
    # r.v <- as.polygons(r.all,dissolve=F,values=F)
    # r.v$ID <- 1:ncell(x)
    # plot(x)
    # plot(v,add=T)
    # plot(v[dangerLinks$ID,],add=T,col='red')
    # text(r.v,label=r.v$ID)
    # 
    # adj.sub1 <- adj[!(adj[,1] %in% v[dangerLinks$ID,]$from &
    #                     adj[,2] %in% v[dangerLinks$ID,]$to),]
    # 
    # adj.sub2 <- adj.sub1[!(adj.sub1[,2] %in% v[dangerLinks$ID,]$from &
    #                          adj.sub1[,1] %in% v[dangerLinks$ID,]$to),]
    # 
    # adj <- adj.sub2
  }
  
  # option 3 - use neighbours to work out which links cross 
  {
    # neighboursData <- createMovementRules(neighbours)
    # required.cells <- neighboursData$required.cells
    # neighbours.all <- neighboursData$neighbours.all
    # library(data.table)
    # adj.danger <- as.data.table(adj[adj[,1] %in% dangerCells.all,])
    # # idea is to whittle away DT, leaving only those links that are bad
    # # first, remove any links to immediate neighbour cells of the from cells, as these are OK
    # danger.neigh.queen <- as.data.table(terra::adjacent(x, unique(adj.danger$from),directions='queen',
    #                 pairs=T))
    # adj.danger.c <- adj.danger[!danger.neigh.queen, on=c("from", "to")]
    # setnames(adj.danger.c,c('from','to'),c('from_raster','to_raster'))
    # 
    # # now check 81 neighbours of remaining from cells
    # # need to translate original cell values of raster to relative cell values, i.e. cell
    # # number of 81 in relation to the 'from' cel
    # # col names relate to relative cell numbers i.e. those in required cells object below
    # x.dt <- data.table(cell = 1:ncell(x),
    #                    values = x.values)
    # 
    # danger.neigh <- as.data.table(terra::adjacent(x, unique(adj.danger.c$from),
    #                                               directions=neighbours.all,
    #                                               pairs=F)) %>% 
    #   .[,from := unique(adj.danger.c$from)] %>% 
    #   melt(., id.vars='from',
    #        variable.name = 'to_relative',
    #        value.name = 'to') %>% 
    #   .[, to_relative := as.integer(str_replace(to_relative,'V',''))]
    # 
    # setkey(danger.neigh, 'to')
    # setkey(x.dt, 'cell')
    # 
    # danger.neigh.vals <- x.dt[danger.neigh]
    # danger.neigh.vals[is.na(values.r), values.r := NaN]
    # setnames(danger.neigh.vals, c('cell','values.r','from'),
    #          c('to_raster','value_raster','from_raster'))
    # danger.neigh.vals <- danger.neigh.vals[,c(3,4,1,2)]
    # 
    # adj.danger.j <- danger.neigh.vals[adj.danger.c, on=c('from_raster','to_raster')]
    # adj.danger.j.required <- required.cells[adj.danger.j,on='to_relative',allow.cartesian=TRUE]
    # adj.danger.j.required$value_raster <- NULL
    # adj.danger.j.required$to_relative <- NULL
    # setnames(adj.danger.j.required,'required_relative','to_relative')
    # 
    # required.vals <- adj.danger.j.required[danger.neigh.vals,on=c('from_raster','to_relative')] %>% 
    #   .[!is.na(to_raster),]
    # 
    # required.vals.s <- 
    #   required.vals[, .(required_cell_nan_count = sum(is.nan(value_raster)),
    #                     num_required = .N),
    #                 by=list(from_raster,to_raster)] %>% 
    #   .[,required_cell_nan_per := required_cell_nan_count / num_required]
    # 
    # badlinks <- required.vals.s[required_cell_nan_per >= 0.25,]
    # 
    # # show links
    # fromXY <- xyFromCell(x, badlinks$from) %>%
    #   cbind(.,object = 1:nrow(.))
    # 
    # toXY <- xyFromCell(x, badlinks$to) %>%
    #   cbind(.,object = 1:nrow(.))
    # 
    # pointMat <- rbind(fromXY, toXY) %>%
    #   cbind(.,part = 1) %>%
    #   .[,c(3,4,1,2)] %>%
    #   .[order(.[,1]),]
    # # create vector lines showing path of links
    # v <- vect(pointMat,'lines',crs="EPSG:27700")
    # 
    # r.all <- x
    # r.all[] <- 1
    # r.v <- as.polygons(r.all,dissolve=F,values=F)
    # r.v$ID <- 1:ncell(x)
    # plot(x)
    # plot(r.v,add=T)
    # text(r.v,label=r.v$ID)
    # plot(v,add=T,col='red')
    # 
    # plot(mask(v,r.v[r.v$ID==37,]),add=T)
    # 
    # # need to build from / to data.table with a link expressed as indices of neighbours_81
    # # showing 'killer' cells, i.e. those cells which if NaN would mean the link is invalid
    # 
    # r.all <- x
    # r.all[] <- 1
    # r.v <- as.polygons(r.all,dissolve=F,values=F)
    # r.v$ID <- 1:ncell(x)
    # plot(x)
    # plot(r.v,add=T)
    # text(r.v,label=r.v$ID)
    
  }
  
  elev_values <- terra::values(x)[,1]
  
  message("calculating slope...")
  
  rise <- (elev_values[adj[,2]] - elev_values[adj[,1]])
  run <- calculate_distance(x = x, adj = adj)
  
  mathematical_slope <- rise/run
  
  if(exaggeration) { 
    mathematical_slope <- ifelse(mathematical_slope > 0, mathematical_slope * 1.99, mathematical_slope * 2.31)
  }
  
  ncells <- length(cells) + length(na_cells)
  
  cf <- cost(cost_function = cost_function, crit_slope = crit_slope, percentile = percentile)
  
  if(is.function(cost_function)) { 
    message(c("Applying ", deparse(body(cost_function)[[2]]), " cost function"))
  } else{ 
    message(c("Applying ", cost_function, " cost function"))
  }
  
  speed <- cf(mathematical_slope)
  
  conductance <- speed/run
  
  if(!is.null(max_slope)) {
    max_slope <- max_slope/100
    index <- abs(mathematical_slope) >= max_slope
    conductance[index] <- 0
  }
  
  cs_matrix <- Matrix::Matrix(data = 0, nrow = ncells, ncol = ncells, sparse = TRUE)
  cs_matrix[adj] <- conductance
  
  cs <- list("conductanceMatrix" = cs_matrix, 
             "costFunction" = cost_function,
             "maxSlope" = ifelse(!is.null(max_slope), paste0(max_slope*100, "%"), NA), 
             "exaggeration" = exaggeration,
             "criticalSlope" = ifelse(test = !is.function(cost_function), yes = ifelse(test = cost_function == "wheeled transport", yes = paste0(max_slope, "%"), no = NA), no = NA),
             "percentile" = ifelse(test = !is.function(cost_function), yes = ifelse(test = cost_function == "campbell 2019", yes = percentile, no = NA), no = NA),
             "neighbours" = sum(neighbours, na.rm = TRUE),
             "resolution" = terra::res(x), 
             "nrow" = terra::nrow(x), 
             "ncol" = terra::ncol(x), 
             "extent" = x@ptr$extent$vector, 
             "crs" = terra::crs(x, proj = TRUE))
  
  message('preparing igraphs...')
  cs_rast <- terra::rast(nrow = cs$nrow, 
                         ncol = cs$ncol, 
                         xmin = cs$extent[1], 
                         xmax = cs$extent[2], 
                         ymin = cs$extent[3], 
                         ymax = cs$extent[4],
                         crs = cs$crs)
  
  col_sum <- Matrix::colSums(cs$conductanceMatrix)
  row_sum <- Matrix::rowSums(cs$conductanceMatrix)
  logical_sm <- methods::as(cs$conductance, "lMatrix")
  ncols <- Matrix::colSums(logical_sm)
  nrows <- Matrix::rowSums(logical_sm)
  vals <- ((col_sum / ncols) + (row_sum / nrows)) / 2
  
  cs_rast.wvals <- terra::setValues(cs_rast, vals)
  cs_rast_df <- terra::as.data.frame(cs_rast.wvals,xy=T,cells=T)
  cs_rast_df[,4] <- NULL
  
  cm_graph <- igraph::graph_from_adjacency_matrix(cs$conductanceMatrix,
                                                  mode = "directed", weighted = TRUE)
  
  igraph::E(cm_graph)$weight <- (1/igraph::E(cm_graph)$weight)
  
  cm_graph_df <- igraph::as_data_frame(cm_graph)
  
  cm_graph_df[is.na(cm_graph_df$weight),]
  cs$all_nodes <- unique(c(cm_graph_df$from,cm_graph_df$to))
  cs$cpp_cm_graph <-  cppRouting::makegraph(cm_graph_df)
  cs$x <-  x
  rm(cm_graph_df)
  gc()
  
  class(cs) <- "conductanceMatrix"
  
  return(cs)
}


getLinkRequirements <- function(neighboursMat) {
  
  neighbours.r <- rast(neighboursMat, crs='EPSG:27700')
  
  neighbours_all <- matrix(rep(1,length(neighboursMat)), nrow = sqrt(length(neighboursMat)),
                          ncol = sqrt(length(neighboursMat)), byrow = TRUE)
  neighbours_all.r <- rast(neighbours_all, crs='EPSG:27700')
  neighbours_all.v <- as.polygons(neighbours_all.r,dissolve=F)
  neighbours_all.v$ID <- 1:length(neighbours_all)
  
  fromXY <- xyFromCell(neighbours_all.r, 
                       rep(ceiling(length(neighbours_all)/2),
                           length(neighbours_all))) %>%
    cbind(.,object = 1:nrow(.))
  
  toXY <- xyFromCell(neighbours_all.r, 1:length(neighbours_all)) %>%
    cbind(.,object = 1:nrow(.))
  
  pointMat <- rbind(fromXY, toXY) %>%
    cbind(.,part = 1) %>%
    .[,c(3,4,1,2)] %>%
    .[order(.[,1]),]
  
  # create vector lines showing path of links
  v <- vect(pointMat,'lines',crs="EPSG:27700")
  v$ID <- 1:nrow(v)
  v$from <- ceiling(length(neighbours_all)/2)
  v$to <- 1:length(neighbours_all)
  plot(neighbours_all.r)
  neighbours.r[25] <- 3
  plot(neighbours.r)
  plot(v,add=T)
  text(neighbours_all.v,neighbours_all.v$ID)
  
  # extract cell numbers for each link
  v.link.cells <- extract(neighbours_all.r, v, cells=T)
  v.link.cells[,2] <- NULL
  # remove reference to centre cell and self
  v.link.required.cells <- 
    v.link.cells[v.link.cells$ID != v.link.cells$cell & 
                   v.link.cells$cell != ceiling(length(neighbours_all)/2),]
  
  v.link.required.cells <- as.data.table(v.link.required.cells)
  names(v.link.required.cells) <- c('to_relative','required_relative')
  # if any of the listed cells per destination are NaN then the link is invalid.
  l <- list(required.cells = v.link.required.cells,
            neighbours.all = neighbours_all)
  return(l)
}

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


# hex.w <- gbdem.hex.filtered.w[[10]]
# points.sf = gbdem.points
# maxPathsPerCore=200000
# refreshCS = F
calcLCP <-
  function(hex.w,
           points.sf,
           params,
           maxPathsPerCore=300000,
           refreshCS = F) {
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
    st <- Sys.time()
    `%+%` <- function(x, y) paste0(x,y)
    `%!in%` <- function(x,y)!('%in%'(x,y))

    st <- Sys.time()
    # setwd(wd)
    print('compute node system tempdir is...')
    sysTmpDir <- Sys.getenv("TMPDIR")
    print(sysTmpDir)

    ##•┣ define func ----

    # idx <- nodesMaxList$`59`
    getLCPsLines <- function(idx, hex.id, hex.r,
                             cpp_cm_graph,
                             point.nodes.clean,
                             params) {

      st <- Sys.time()
      `%+%` <- function(x, y) paste0(x,y)
      `%!in%` <- function(x,y)!('%in%'(x,y))

      nodeBatchNumber <- paste0(idx[1],'-',idx[length(idx)])

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
                              crs=terra::crs(hex.r)
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

      print(glue::glue('finished calculating shortest paths for hex {hex.id}, ' %+%
                         'node batch {nodeBatchNumber}...'))

      print('writing raster...')
      lcpOut <- str_replace(gpkgOut,'.gpkg','.tif')
      hex.exent <- ext(hex.r)
      lcpsExt.v <- c(hex.exent[1],hex.exent[3],hex.exent[2],hex.exent[4])

      print('rasterizing lcps...')
      system(glue::glue(
        "gdal_rasterize -te {paste0(lcpsExt.v,collapse=' ')} " %+%
          "-tr {paste0(params$demRes,collapse=' ')} -ot GTiff " %+%
          "-burn 1 -init 0 -add -at {gpkgOut} {lcpOut}"
      )
      )
      

      print('summing hex raster with new lcp raster...')
      hexRasOut <- file.path(getwd(),'outputs',glue::glue('hex_{hex.id}.tif'))
      hexRasOutTmp <- file.path(getwd(),'outputs',glue::glue('hex_{hex.id}_tmp.tif'))
      
      system(glue::glue('gdal_calc.py -A {hexRasOut} -B {lcpOut} --outfile={hexRasOutTmp} --calc="A+B"'))
      file.copy(hexRasOutTmp, hexRasOut, overwrite=T)
      
      # tmp.r <- rast(hexRasOut)
      # lcp.r <- rast(lcpOut)
      # terra::crs(lcp.r) <- terra::crs(tmp.r)
      # sum.r <- tmp.r + lcp.r
      # terra::writeRaster(sum.r, filename=hexRasOut, overwrite=T)
      
      file.remove(gpkgOut)
      file.remove(lcpOut)

      diff <- difftime(Sys.time(),st)
      print(glue::glue('finished all calculations for hex {hex.id}, ' %+%
                         'node batch {nodeBatchNumber}'))
      print(diff)
      return(diff)
    }

    ##•┣ run script ----
    # hex.w <- gbdem.hex.filtered.w[[30]]
    # first extract the hex polygon
    hex.pol <- terra::vect(hex.w)
    hex.id <- hex.pol$hex_id


    sink(glue('logs/calc_lcp_sinkout_hex_{hex.id}.txt'))

    # subset points and convert to spatvector
    hex.points <- points.sf %>%
      dplyr::filter(hex_id == hex.id) %>%
      vect()

    # load dem and mask / trim
    if (refreshCS==T | !file.exists(glue::glue('bigdata/hex_{hex.id}_cs.RDS'))) {
      print(glue::glue('trimming main raster to hex {hex.id} area...'))
      gb.dem <- rast(params$demFile)

      hex.r <- terra::trim(terra::mask(gb.dem,hex.pol))

      # create cost surface
      print('creating cost surface...')
      cs <-
        create_slope_cs2(hex.r,
                         cost_function =
                           "campbell 2019",
                         neighbours = params$lcpNeigh)
      # cs.orig <- cs
      cs$x <- wrap(cs$x)
      saveRDS(cs,glue::glue('bigdata/hex_{hex.id}_cs.RDS'))
      print(glue::glue('cost surface saved to bigdata/hex_{hex.id}_cs.RDS'))
    }

    print('loading cost surface...')
    cs <- readRDS(glue::glue('bigdata/hex_{hex.id}_cs.RDS'))
    hex.r <- rast(cs$x)
    cpp_cm_graph <- cs$cpp_cm_graph

    # TODO: might be worth experimenting with {dodgr} to see if any speed gains

    print(glue::glue('writing base raster for hex {hex.id}...'))
    
    # required for gdal procedure only
    hexRasOut <- file.path(getwd(),'outputs',glue::glue('hex_{hex.id}.tif'))
    
    hex.r.zero <- hex.r
    hex.r.zero[!is.nan(hex.r.zero)] <- 0
    terra::writeRaster(hex.r.zero, filename=hexRasOut,
                       overwrite=T)

    # set up lcp data...
    # get cell numbers of points
    point.nodes <- terra::extract(hex.r, hex.points, cells=T)
    # check if any points fall on NaN cells, and remove if so
    point.nodes.clean <- point.nodes[!is.na(point.nodes[,2]),]$cell

    # split all terminus points into correct sized list
    maxPointsPerCore <- floor(maxPathsPerCore/length(point.nodes.clean))
    nodesIdx <- 1:length(point.nodes.clean)
    nodesMaxList <- split(nodesIdx,
                          ceiling(seq_along(nodesIdx) /
                                    maxPointsPerCore))

    # # check for pre-existing hex and node batches
    hexBatchOutputs <- nodesMaxList %>%
      map_chr(.f = function(idx) {
        nodeBatchNumber <- paste0(idx[1],'-',idx[length(idx)])
        tifOut <- file.path(getwd(),'outputs',glue::glue("hex_{hex.id}_batch_{nodeBatchNumber}.tif"))
      })

    toProcess <-
      which(hexBatchOutputs %!in%
              list.files(file.path(getwd(),'outputs'),full.names = T))
    
    print(glue::glue('calculating {length(nodesMaxList)} lcp node batches for hex {hex.id}...'))
    ##•┣ call func ----
    lapply(nodesMaxList[toProcess],
           getLCPsLines,
           hex.id,
           hex.r,
           cpp_cm_graph,
           point.nodes.clean,
           params)
    
    # file.remove(hexRasOutTmp)

    diff <- difftime(Sys.time(),st)
    print(glue::glue('finished calculating least cost paths for hex {hex.id} in {diff}...'))
    sink()
    return(diff)
  }



