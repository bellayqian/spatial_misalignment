# Load required package
library(sp)
library(sf)
library(tidyverse)
library(ggplot2)

# Define the resolution of the first grid (coarser resolution)
resolution1 <- c(5, 5)  # 5x5 grid

# Define the resolution of the second grid (finer resolution)
resolution2 <- c(10, 10)  # 10x10 grid

# Create the first spatial grid (coarser grid)
grid1 <- GridTopology(cellcentre.offset = c(0.5, 0.5), 
                      cellsize = c(2, 2),  # Cell size adjusted to match extent and resolution
                      cells.dim = resolution1)

sp_grid1 <- SpatialGrid(grid1)

# Create the second spatial grid (finer grid)
grid2 <- GridTopology(cellcentre.offset = c(0, 0), 
                      cellsize = c(1, 1),  # Smaller cell size for finer resolution
                      cells.dim = resolution2)

sp_grid2 <- SpatialGrid(grid2)

# Plot both grids to visualize the differences
plot(sp_grid2, col = "lightblue", main = "Spatial Grids with Different Resolutions")
plot(sp_grid1, col = 'red', add = TRUE)  # Second grid with finer resolution shown in red border

## now union together cells of the finer grid to create non-nested misalignment (there's probably a simpler/cleaner way to do this)
sp_grid2_poly<-st_as_sf(as(sp_grid2,'SpatialPolygons'))
sp_grid1_poly<-st_as_sf(as(sp_grid1,'SpatialPolygons'))

## add ids and plot to check ordering of grid cells in the dataset
## add cell ids
# sp_grid2_poly$ID<-1:nrow(sp_grid2_poly)
# ggplot() + 
#   geom_sf(data = sp_grid2_poly, aes(fill = ID))

## pick cells to cross-union (cell to union that cross boundaries of coarser grid)
## 2,3 or 4,5 or 6,7 or 8,9
## pick one of these every other row
set.seed(1)
pair_list<-list(c(2,3),c(4,5),c(6,7),c(8,9))
union_ids<-matrix(0,10,10)
id_num<-1
for (i in seq(1,10,2)){
  union_ids[i,pair_list[[sample(1:4,size=1)]]]<-id_num
  id_num<-id_num+1
}

sp_grid2_poly$union_ids<-c(t(union_ids))

store_xunions<-NULL
for (i in 1:max(sp_grid2_poly$union_ids)){
  temp<-sp_grid2_poly[sp_grid2_poly$union_ids==i, ] %>% 
    st_union() %>% # unite to a single geometry object
    st_sf() # make the geometry a data frame object
  
  store_xunions<-rbind(store_xunions,temp)
}

## now merge together remaining cells of finer grid to get as close as possible to grid 1
grid2_nounion<-subset(sp_grid2_poly,union_ids==0)
## find cells that intersect with each cell of grid 1 and union them
store_iunions<-NULL
for (i in 1:nrow(sp_grid1_poly)){
  temp<-grid2_nounion[c(st_contains(sp_grid1_poly[i,],grid2_nounion,sparse=F))==T,] %>% 
    st_union() %>% # unite to a single geometry object
    st_sf() # make the geometry a data frame object
  
  store_iunions<-rbind(store_iunions,temp)
}

grid2_final<-rbind(store_iunions,store_xunions)

plot(st_geometry(grid2_final),border='red')
plot(st_geometry(sp_grid1_poly),add=T)

## rename and save grids
gridx<-grid2_final
gridy<-sp_grid1_poly

save(gridx,gridy,
  file='sim_grids.RData')