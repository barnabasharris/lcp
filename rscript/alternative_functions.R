# quickly counts occurences of values from vector
countOccurences = function(x) {
  data.table::data.table(x)[, .N, keyby = x]
}

# this calculates the routes, extracts the raster cell numbers and counts
# and directly applies these to a pre-existing raster for speed. Problem
# is that some cells involved in the path are not identified (i.e. those
# that are only clipped by the putative path line). Partly addressed by
# buffering round the cells but produces noisier dataset overall.
getLCPsNoLines <- function(nodeBatchNumber, nodesMaxList) {
  
  idx <- nodesMaxList[[nodeBatchNumber]]
  
  print(glue::glue('calculating shortest paths for hex {hex.pol$hex_id}, node batch {nodeBatchNumber}...'))
  
  lcp_cells_list <- cppRouting::get_multi_paths(cs$cpp_cm_graph,
                                                from=point.nodes.clean[idx],
                                                to=point.nodes.clean,
                                                long=F)
  pathCellsList <- lcp_cells_list %>% 
    map(~.Internal(unlist(.x, FALSE, FALSE)))
  
  pathCells <- 
    as.integer(.Internal(unlist(pathCellsList, FALSE, FALSE)))
  
  occurences <- countOccurences(pathCells)
  setnames(occurences, 'x', 'cell')
  
  # buffer occurences
  adj <- adjacent(hex.r,occurences$cell,direction='8') %>%
    as.data.table() %>%
    .[,cell := occurences$cell] %>%
    .[occurences, on='cell'] %>%
    data.table::melt(id.vars = 'N',
                     measure.vars = 1:8) %>%
    .[,c('value','N')] %>%
    setnames(., 'value','cell')
  
  # bind to path cells dt
  occurences.adj <- rbind(occurences,adj)
  
  # sum any repeats
  occurences.adj.s <- occurences.adj[, .(N = sum(N)), by=cell]
  
  return(occurences.adj.s)
}


# gen data.table with all cells and xys
rasDt <- data.table(cell = 1:ncell(hex.r),
                    allXy = terra::xyFromCell(hex.r,1:ncell(hex.r))
)

# potentially faster way to create lines
dt <- rbindlist(lapply(f,as.data.table), idcol = 'from')
xy <- xyFromCell(hex.r, as.numeric(dt$V1))
dt$x <- xy[,1]
dt$y <- xy[,2]

setnames(dt,'from','object')
setnames(dt,'V1','part')
dt.m <- matrix(as.numeric(unlist(dt)),nrow=nrow(dt))
lcps <- terra::vect(x=dt.m,
                    crs=crs(hex.r)
                    ,type='lines'
)
writeVector(lcps,'test.gpkg')


gb <- rast(params$demFile)
gb[is.nan(gb)] <- 99999
gb.nan <- gb[gb == 99999] 
gb.nan.buff <- buffer(gb.nan,150)

