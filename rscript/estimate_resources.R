# plot(hex.r)
# plot(hex.points[!is.na(point.nodes[,2])],add=T)
mp <- 50000
mp <- 100000
peakRams <- seq(50000,500000,by=25000) %>% 
  map(.f = function(mp) {
    
    # maxPointsPerCore <- floor(params$maxPathsPerCore/length(point.nodes.clean))
    maxPointsPerCore <- floor(mp/length(point.nodes.clean))
    
    nodesIdx <- 1:length(point.nodes.clean)
    nodesMaxList <- split(nodesIdx, 
                          ceiling(seq_along(nodesIdx) / 
                                    maxPointsPerCore))
    
    # fast approach 
    
    ram <- peakRAM::peakRAM({
      allOccurences <- 1:length(nodesMaxList[1:3]) %>% 
        map(~getLCPsNoLines(.x,nodesMaxList))
    })
    return(ram)
  })

mpTimes <- map2_df(.x =  seq(50000,500000,by=25000), .y = peakRams, function(mp, y) {
  maxPointsPerCore <- floor(params$maxPathsPerCore/length(point.nodes.clean))
  maxPointsPerCore <- floor(mp/length(point.nodes.clean))
  
  nodesIdx <- 1:length(point.nodes.clean)
  nodesMaxList <- split(nodesIdx, 
                        ceiling(seq_along(nodesIdx) / 
                                  maxPointsPerCore))
  
  data.frame(mp = mp, secsperhex = (y$Elapsed_Time_sec/3) * length(nodesMaxList),
             peakram = y$Peak_RAM_Used_MiB)
})

load('bigdata/mpTimes.RData')

secsPerHex <- (300/5)*length(nodesMaxList) # maxNodes = 250000
secsPerHex <- (210/2)*length(nodesMaxList)
minsPerHex <- secsPerHex/60
hoursPerHex <- minsPerHex/60
totalProcessHours <- hoursPerHex*nrow(gbdem.hex.filtered)
HourPerCore <- totalProcessHours / 40 