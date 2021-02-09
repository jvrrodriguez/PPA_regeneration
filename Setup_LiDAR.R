#library(rlas)
#library(ForestTools) 
library(lidR)
library(rLiDAR)

library(sp)
library(rgdal)
library(shapefiles)
library(spatstat)
library(maptools)
library(raster)

#https://cran.r-project.org/web/packages/lidR/readme/README.html
#https://github.com/Jean-Romain/lidR/wiki/Rasterizing-perfect-canopy-height-models
#https://github.com/gisma/uavRst/wiki/Building-a-Canopy-Height-Model-(CHM)-using-lidR
#https://cran.r-project.org/web/packages/ForestTools/vignettes/treetopAnalysis.html
#https://rdrr.io/cran/ForestTools/f/vignettes/treetopAnalysis.Rmd

all.plot = list.files("~/Documentos/Datos NO publicados/BioIntForest/Data/Plots/",pattern = "^parcela.*\\poli.shp$")
all.plot = gsub(".shp","" , all.plot ,ignore.case = TRUE)
setwd("~/Documentos/Datos NO publicados/BioIntForest/Data/Plots")
#quadrats <- readOGR(dsn=".",layer="~/Documentos/Datos NO publicados/BioIntForest/Data/Plots/parcelas.shp")

all.las = list.files("~/Documentos/Datos NO publicados/BioIntForest/Data/LiDAR/",pattern = "2011.las")

Plot.LiDAR <- FALSE
LAD.LiDAR <- TRUE
Trees.LiDAR <- FALSE


for (j in 1:length(all.plot)) {
  
  ETRS89 <- CRS("+proj=utm +zone=30 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
  WGS84 <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  
  
  #=======================================================================#
  # Importing LAS file:
  #=======================================================================#
  
  setwd("~/Documentos/Datos NO publicados/BioIntForest/Data/LiDAR")
  dataLas <- rLiDAR::readLAS(paste0("PuntosParcela", j, "-2011.las"), short=TRUE)
  
  setwd("~/Documentos/Datos NO publicados/BioIntForest/Data/Plots")
  plot   <- readOGR(dsn=".",layer=all.plot[j])
  
  if (projection(plot) != as.character(ETRS89)) plot <- spTransform(plot, CRS = ETRS89) #2, 3, 5, 7, 8, 9 son WGS84
  
  # Summary of the LAS file
  summary(dataLas)
  
  
  #=======================================================================#
  # LAS file visualization:
  #=======================================================================#
  
  if (Plot.LiDAR == TRUE) {
    
    # 01 Set a single color
    col<-"forestgreen"
    
    # plot 2D
    plot(dataLas[,1],dataLas[,2], col=col,xlab="UTM.Easting", ylab="UTM.Northing", main="Single color")
    
    # plot 3D
    library(rgl)
    points3d(dataLas[,1:3], col=col, axes=FALSE,xlab="", ylab="", zlab="")
    axes3d(c("x+", "y-", "z-"))                 # axes
    grid3d(side=c('x+','y-','z'), col="gray")   # grid
    title3d(xlab = "UTM.Easting", ylab = "UTM.Northing",zlab = "Height(m)", col="red") # title
    planes3d(0, 0, -1, 0.001, col="gray", alpha=0.7)   # terrain
    
    # 02 Set a color by height# color ramp
    
    myColorRamp <- function(colors, values) {
      v <- (values - min(values))/diff(range(values))
      x <- colorRamp(colors)(v)
      rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
    }
    
    # Color by height
    col <- myColorRamp(c("blue","green","yellow","red"), dataLas[,3])
    
    # plot 2D
    plot(dataLas[,1], dataLas[,2], col=col, xlab="UTM.Easting", ylab="UTM.Northing", main="Color by height")
    
    # plot 3D
    points3d(dataLas[,1:3], col=col, axes=FALSE, xlab="", ylab="", zlab="")
    axes3d(c("x+", "y-", "z-"))               # axes
    grid3d(side=c('x+','y-','z'), col="gray") # grid
    title3d(xlab = "UTM.Easting", ylab = "UTM.Northing",zlab = "Height(m)", col="red") # title
    planes3d(0, 0, -1, 0.001, col="gray",alpha=0.7) # terrain
    
  }
    
  
  #=======================================================================#
  # Within-plot canopy structure
  #=======================================================================#
  
  if (LAD.LiDAR == TRUE) {
  
    dataLas.sp <- SpatialPointsDataFrame(coordinates(dataLas[,1:2]), data.frame(z=dataLas[,3]) )
    crs(dataLas.sp) <- ETRS89
    
    coords = matrix(c(0,0, 0,2*20, 2*14.99,2*20, 2*14.99,0), # antes 14.5 pero haciendo analisis falta un cacho 
                    ncol = 2, byrow = TRUE) #excepto la A1 y A3, que van horizontales
    P1 = Polygon(coords)
    Plot = SpatialPolygons(list(Polygons(list(P1), ID = "a")))
    
    CuadCode <- paste0(as.vector(t(replicate(15, LETTERS[1:15]))), sprintf("%02d",rep(1:20, 15)))
    Points <- SpatialPointsDataFrame(makegrid(Plot, cellsize=c(2,2)), data.frame(CuadCode))
    
    setwd("~/Documentos/Datos NO publicados/BioIntForest/Analysis/")
    source("function rotated.grid.R")
    
    plot.grid <- make_grid_rot(plot, cellsize = c(2,2), rotation = 1)
    if (length(plot.grid) != length(Points)) plot.grid <- plot.grid[1:length(Points),]
    plot(plot, col = "blue")
    plot(plot.grid, add = TRUE)
    
    dataLas.cell <- sp::over(plot.grid, dataLas.sp, returnList = TRUE)
  
    lad <- matrix(0, length(dataLas.cell), 3)
    colnames(lad) <- c("lad.var", "lad.cv", "lad.shan") #"lad.mean"
    
    for (i in 1:length(dataLas.cell)) {
      
      lad.tmp <- dataLas.cell[[i]][[1]]
      if (length(lad.tmp)>0) {
        
        lad[i,1] <- var(LAD(lad.tmp, z0 = min(dataLas.sp$z))$lad, na.rm = T)
        lad[i,2] <- sd(LAD(lad.tmp, z0 = min(dataLas.sp$z))$lad, na.rm = T) / mean(LAD(lad.tmp, z0 = min(dataLas.sp$z))$lad, na.rm = T)
        p_i <- LAD(lad.tmp, 2)$lad[-1]/sum((LAD(lad.tmp, 2)$lad[-1]))
        lad[i,3] <- -sum(p_i[p_i!=0] * log(p_i[p_i!=0], base = exp(1)))
        
      }
      lad[i,] <- ifelse(!is.finite(lad[i,]), 0, lad[i,])
    }
    
    Points@data <-  cbind(Points@data, lad)
    
    library(automap)
    
    # This is the new grid for predicting
    Plot.grid <- as(raster(ncol = 15*2, nrow = 20*2, crs=NULL, ext = extent(Points@bbox)), "SpatialPolygons")
    Plot.grid <- SpatialPolygonsDataFrame(Plot.grid, data.frame(n=1:length(Plot.grid)))
    
    Plot.vars <- colnames(lad)
    
    for (i in 1:length(Plot.vars)) {
    
      sph.fit <- autoKrige(formula(paste(Plot.vars[i], "~ 1")), Points, Plot.grid) #, nmax = 1 for binary data
      plot(sph.fit)
      
      colnames(sph.fit$krige_output@data)[3:5] <- paste0(Plot.vars[i], c(".pred", ".var", ".stdev"))
      Plot.grid@data <- cbind(Plot.grid@data, sph.fit$krige_output@data[,3:5])
    
    }
    
    source("~/Documentos/Datos NO publicados/BioIntForest/Analysis/function Raster_extractGrid.R")
    ras.lad <- Raster_extractGrid(Plot.grid, colnames(Plot.grid@data)[grep(".pred", colnames(Plot.grid@data))], 0.5, plot = F)
    
    png(filename = paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/Environ_Maps/Lad_A",j,".png"), width = 800, height = 1000)
    plot(ras.lad)
    dev.off()
    
    source("~/Documentos/Datos NO publicados/BioIntForest/Analysis/function raster.as.im.R")
    im.lad <- listof( lad.var=raster.as.im(ras.lad[[1]]), lad.cv=raster.as.im(ras.lad[[2]]), lad.shan=raster.as.im(ras.lad[[3]]))
    
    save(ras.lad, im.lad, file = paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/Environ_Maps/Lad_A",j,".RData"))
    
  }
  
  #=======================================================================#
  # Individual tree detection using K-means cluster
  #=======================================================================#
  
  if (Trees.LiDAR == TRUE) {
    
    # Setting the xyz coordinates and subsetting the data
    xyz<-subset(rLAS[,1:3],rLAS[,3] >= 1.37)
    
    # Finding clusters (trees)
    clLAS<-kmeans(xyz, 50)
    
    # Set the tree id vector
    id<-as.factor(clLAS$cluster) #PLEASE NAME AS id!!!!
    
    # Combining xyzi and tree id
    xyId<-cbind(xyz[,1:2],id)
    
    # Compute the lidar convex hull of the clusters
    chullTrees<-chullLiDAR2D(xyId)
    
    # Plotting the lidar convex hull
    library(sp)
    
    plot(SpatialPoints(xyId[,1:2]),cex=0.5,col=xyId[,3], axes=T)
    plot(chullTrees$chullPolygon,add=TRUE, border="green")
    
    # Get the ground-projected area of lidar convex hull
    chullList< - chullTrees$chullArea
    summary(chullList) # summary
    
    
    #=======================================================================#
    #  Computing individual tree LiDAR metrics
    #=======================================================================#
    
    TreesMetrics <- CrownMetrics(xyziId)
    head(TreesMetrics)
    
    # Compute the lidar convex hull of the clusters 
    chullTrees <- chullLiDAR2D(xyziId)
    
    
    #===============================================================#
    # Smoothing the CHM using a Gaussian filter
    #==================================================================#
    
    # Smoothing CHM
    sCHM<-CHMsmoothing(chm, "Gaussian", 3, 0.6)
    plot(sCHM)
    
    setwd("~/Documentos/Datos NO publicados/BioIntForest/Data/LiDAR")
    las <- lidR::readLAS(paste0("PuntosVegetacionParcela", j, ".las"))
    
    # Basic triangulation and rasterization of first returns
    chm <- grid_canopy(las, res = 0.25, pitfree(max_edge = c(3, 1.5)))
    plot(chm, col = height.colors(50))
    
    #https://gis.stackexchange.com/questions/183175/rotating-90-using-two-point-equidistant-projection-with-proj4
    # rotate <- function(x, angle=0, resolution=res(x)) {
    #   y <- x; crs(y) <- "+proj=aeqd +ellps=sphere +lat_0=90 +lon_0=0"
    #   projectRaster(y, res=resolution, 
    #                 crs=paste0("+proj=aeqd +ellps=sphere +lat_0=90 +lon_0=", -angle))
    # }
    
    plot.corners <- plot@polygons[[1]]@Polygons[[1]]@coords
    
    plot.angle = 90 - (180/pi) * atan((plot.corners[2,1] - plot.corners[1,1]) / (plot.corners[2,2] - plot.corners[1,2]))
    
    #chm.rot <- rotate(chm, rectangle_angle - 90, resolution = 0.25) #aprox ese angulo
    
    #https://cran.r-project.org/web/packages/tangles/vignettes/deidentification.html
    
    source("~/Documentos/Datos NO publicados/BioIntForest/Analysis/function raster.as.im.R")
    im.canopyLidar <- raster.as.im(chm)
    im.canopyLidar <- shift(im.canopyLidar, origin="bottomleft")
    
    im.canopyLidar <- rotate.im(im.canopyLidar, angle = plot.angle, centre="centroid")
    plot(im.canopyLidar)
    
    save(im.canopyLidar, file = paste0("~/Documentos/Datos NO publicados/BioIntForest/Analysis/LiDAR_A",j,"_Canopy.RData"))
    
    #las <- lastrees(las, li2012())
    #col <- random.colors(200)
    #plot(las, color = "treeID", colorPalette = col)
    
    #lidR::plot(las,color="Z",colorPalette = pal(100),backend="pcv")
    
  }
  
}