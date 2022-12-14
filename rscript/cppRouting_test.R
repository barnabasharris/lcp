library(cppRouting)
library(dplyr)
library(sf)
library(ggplot2)
library(reshape2)


#Reading french road data

roads <- read.csv("/media/mal/git/cppRouting/data_readme/roads.csv",
                  colClasses = c("character","character","numeric"))
#Shapefile data of communes (polygons)
com<-read_sf("/media/mal/git/cppRouting/data_readme/com_simplified_geom.shp")

#General practitioners locations
med<-read.csv("/media/mal/git/cppRouting/data_readme/doctor.csv",colClasses = c("character","numeric","character","numeric"))
#Import materinty ward locations
maternity<-read.csv("/media/mal/git/cppRouting/data_readme/maternity.csv",colClasses = c("character","numeric"))
#Commuting data from national census
load("/media/mal/git/cppRouting/data_readme/commuting.Rds")
#Import nodes coordinates (projected in EPSG : 2154)
coord<-read.csv("/media/mal/git/cppRouting/data_readme/coordinates.csv",colClasses = c("character","numeric","numeric"))

graph <- makegraph(roads,directed = T,coords = coord)

#Compute shortest path
trajet<-get_path_pair(cpp_cm_graph,from="205793",to="212490")

#Compute shortest path
distance<-get_distance_pair(graph,from="205793",to="212490")

#Compute detour time of 25 and 45 minutes
det25<-get_detour(graph,from="205793",to="212490",extra=25)
det45<-get_detour(graph,from="205793",to="212490",extra=45)



#Create sf object of nodes
pts<-st_as_sf(coord,coords=c("X","Y"),crs=2154)
pts<-st_transform(pts,crs=4326)
pts$time<-ifelse(pts$ID %in% unlist(det45),"45","0")
pts$time<-ifelse(pts$ID %in% unlist(det25),"25",pts$time)
pts$time<-ifelse(pts$ID %in% unlist(trajet),"Shortest Path",pts$time)
pts$time<-factor(pts$time,levels = c("25","45","Shortest Path","0"))

#Plot
dijon=get_map(location=c(lon=5.041140,lat=46.48),zoom=8, source="google",maptype = "toner-2010")

p <- ggplot() +
  geom_sf(data=pts[pts$time!="0",],aes(color=time),inherit.aes = FALSE)+
  ggtitle(paste0("Detours around Dijon-lyon path - ",round(distance,digits = 2)," minutes"))+
  labs(color="Minutes")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y=element_blank(),axis.title.x=element_blank())
p

##

#Data describing edges of the graph
edges<-data.frame(from_vertex=c(0,0,1,1,2,2,3,4,4), 
                  to_vertex=c(1,3,2,4,4,5,1,3,5), 
                  cost=c(9,2,11,3,5,12,4,1,6))

#Get all nodes
nodes<-unique(c(edges$from_vertex,edges$to_vertex))



#Construct directed and undirected graph 
directed_graph<-makegraph(edges,directed=TRUE)
non_directed<-makegraph(edges,directed=FALSE)

#Sampling origin and destination nodes
origin<-sample(nodes,10,replace=TRUE)
destination<-sample(nodes,10,replace=TRUE)

#Get distance between origin and destination in the two graphs
dir_dist<-get_distance_pair(Graph=directed_graph, from=origin, to=destination, allcores=FALSE)
non_dir_dist<-get_distance_pair(Graph=non_directed, from=origin, to=destination, allcores=FALSE)
print(dir_dist)
print(non_dir_dist)

dir_path <- get_path_pair(Graph=directed_graph, from=0, to=5)

tcoords <- data.frame(X = c(1,2,3,1,2,3),
                      Y = c(1,2,1,2,1,2),
                      node = nodes)

pts <- st_as_sf(tcoords,coords=c("X","Y"),crs=27700)

pts$time<-ifelse(pts$node %in% unlist(dir_path),"Shortest Path",NA)
pts$time<-factor(pts$time,levels = c("25","45","Shortest Path","0"))


##

















x <- slope_cs
cs_rast <- terra::rast(nrow = x$nrow, ncol = x$ncol, xmin = x$extent[1], xmax = x$extent[2], ymin = x$extent[3], ymax = x$extent[4],crs = x$crs)
col_sum <- Matrix::colSums(x$conductanceMatrix)
row_sum <- Matrix::rowSums(x$conductanceMatrix)
logical_sm <- methods::as(x$conductance, "lMatrix")
ncols <- Matrix::colSums(logical_sm)
nrows <- Matrix::rowSums(logical_sm)
vals <- ((col_sum / ncols) + (row_sum / nrows)) / 2
cs_rast <- terra::setValues(cs_rast, vals)
library(igraph)

img <- terra::as.data.frame(cs_rast,xy=T)

# make a 2D lattice graph, with same dimensions as "img"
l <- make_lattice(dimvector = c(length(unique(img$y)), 
                                length(unique(img$x))), directed=F, circular=F)
summary(l)
plot(l)
# set vertex attributes
V(l)$x = img$x
V(l)$y = img$y
V(l)$v = img$lyr.1


# compute edge weights as a function of attributes of the two connected vertices
el <- get.edgelist(l)

# "weight" is a known edge attribute, and is used in shortest_path()
# I was confused about weights... lower weights are better, Inf weights will be avoided.
# also note from help: "if all weights are positive, then Dijkstra's algorithm is used."
E(l)$weight <- 1/(pmax(V(l)[el[, 1]]$v, V(l)[el[, 2]]$v))

# E(l)$color = grey.colors(length(unique(E(l)$weight)))[E(l)$weight+1]

idf <- igraph::as_data_frame(l)
df$ID <- 1:nrow(df)
nCoords <- df %>% 
  dplyr::select(ID,x,y) %>% 
  dplyr::rename('X' = 2,
                'Y' = 3)

# "color" is a
graph <- makegraph(idf,directed = F,coords = nCoords)








df <- data.frame(x=seq(0,12,by=.01), 
                 y=sapply(seq(0,12,by=.01), FUN = function(i) 10*sin(i)+rnorm(1)))

# convert to "pixels" of raster data
# assumption: image color is greyscale, only need one numeric value, v
img <- data.frame(table(round(df$y,0), round(df$x,1)))
names(img) <- c("y","x","v")
img$y <- as.numeric(as.character(img$y))
img$x <- as.numeric(as.character(img$x))

# make a 2D lattice graph, with same dimensions as "img"
l <- make_lattice(dimvector = c(length(unique(img$y)), 
                                length(unique(img$x))), directed=F, circular=F)
summary(l)

# set vertex attributes
V(l)$x = img$x
V(l)$y = img$y
V(l)$v = img$v

# compute edge weights as a function of attributes of the two connected vertices
el <- get.edgelist(l)

# "weight" is a known edge attribute, and is used in shortest_path()
# I was confused about weights... lower weights are better, Inf weights will be avoided.
# also note from help: "if all weights are positive, then Dijkstra's algorithm is used."
E(l)$weight <- 1/(pmax(V(l)[el[, 1]]$v, V(l)[el[, 2]]$v))
E(l)$color = grey.colors(length(unique(E(l)$weight)))[E(l)$weight+1]
# find the start/end vertices
start = V(l)[V(l)$x == 0 & V(l)$y == 0]
end = V(l)[V(l)$x == 12 & V(l)$y == -5] 

# get the shortest path, returning "both" (vertices and edges)...
result <- shortest_paths(graph = l, from = start, to = end,  output = "both")

# color the edges that were part of the shortest path green
V(l)$color = ifelse(V(l) %in% result$vpath[[1]], "green", V(l)$color)
E(l)$color = ifelse(E(l) %in% result$epath[[1]], "green", E(l)$color)

# color the start and end vertices red
V(l)$color = ifelse(V(l) %in% c(start,end), "red", V(l)$color)

plot(l, vertex.shape = "square", vertex.size=2, vertex.frame.color=NA, vertex.label=NA, curved=T)



l.df <- igraph::as_data_frame(l)

# l.df$weight[!is.finite(l.df$weight)] <- 1.1
ids <- unique(c(l.df$from,l.df$to))
ncoords <- data.frame(ID = ids, X = img$x, Y=img$y)
l.df$color <- NULL
nl <- makegraph(l.df, ncoords,directed = T)
#Contraction of input graph
graph3<-cpp_contract(nl,silent=F)
?get_distance_pair(nl,1,3000,algorithm = "Dijkstra")


#Data describing edges of the graph
edges<-data.frame(from_vertex=c(0,0,1,1,2,2,3,4,4), 
                  to_vertex=c(1,3,2,4,4,5,1,3,5), 
                  cost=c(9,2,11,3,5,12,4,1,6))

#Get all nodes
nodes<-unique(c(edges$from_vertex,edges$to_vertex))

#Construct directed and undirected graph 
directed_graph<-makegraph(edges,directed=TRUE)
non_directed<-makegraph(edges,directed=FALSE)

#Sampling origin and destination nodes
origin<-sample(nodes,10,replace=TRUE)
destination<-sample(nodes,10,replace=TRUE)

#Get distance between origin and destination in the two graphs
dir_dist <- get_distance_pair(Graph=directed_graph, from=origin, to=destination, allcores=FALSE)

non_dir_dist<-get_distance_pair(Graph=non_directed, from=origin, to=destination, allcores=FALSE)
print(dir_dist)









