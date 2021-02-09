library(sp)
library(rgdal)
library(tidyr)
library(raster)
library(gstat)
library(automap)
library(raster)
library(spatstat)
library(maptools)

# Setup data --------------------------------------------------------------

setwd("~/Documentos/Datos NO publicados/BioIntForest/Data/Inventarios_floristicos")
Data <- data.frame(read.table(file="Cover.csv", header=T, sep=";",dec=","))
Data <- Data  %>% unite("plot", c(Plot, Thinning), remove = FALSE)
# unique(Data$Plot)
# unique(Data$plot)
# unique(Data$CuadCode)

coords = matrix(c(0,0, 0,2*20, 2*14.99,2*20, 2*14.99,0), # antes 14.5 pero haciendo analisis falta un cacho 
                ncol = 2, byrow = TRUE) #excepto la A1 y A3, que van horizontales
P1 = Polygon(coords)
Plot = SpatialPolygons(list(Polygons(list(P1), ID = "a")))
# plot(Plot, axes = TRUE)

CuadCode <- paste0(as.vector(t(replicate(15, LETTERS[1:15]))), sprintf("%02d",rep(1:20, 15)))
Points <- SpatialPointsDataFrame(makegrid(Plot, cellsize = c(2,2)), data.frame(CuadCode))

# No se porque tengo que cambiarlos manualmente
Plot@bbox[,1] <- c(0,0)
Plot@bbox[,2] <- c(30,40)
Points@bbox[,1] <- c(0,0)
Points@bbox[,2] <- c(30,40) 

for (j in 1:length(unique(Data$Plot))) {
  
  Points.data <- subset(Data, Plot == unique(Data$Plot)[j]) %>% pivot_wider(names_from = c(Sp, Year), values_from = Cover)
  dim(Points.data)[1] == dim(Points)[1]
  #rep.points <- which(table(kk$CuadCode) > 1)
  #kk[kk$CuadCode == names(rep.points[2]),]
  
  Plot.tmp <- merge(Points, Points.data, by.x = "CuadCode", by.y = "CuadCode")
  Plot.vars <- apply(expand.grid(unique(Data$Sp), unique(Data$Year)), 1, paste, collapse="_")
  Plot.vars <- Plot.vars[6:10] # only 2008 data
  
  # This is the new grid for predicting
  Plot.grid <- as(raster(ncol = 15*2, nrow = 20*2, crs=NULL, ext = extent(Points@bbox)), "SpatialPolygons") #*4
  Plot.grid <- SpatialPolygonsDataFrame(Plot.grid, data.frame(n=1:length(Plot.grid)))
  
  #https://rstudio-pubs-static.s3.amazonaws.com/80464_9156596afb2e4dcda53e3650a68df82a.html
  
  
  # Now krigging variables --------------------------------------------------
  
  for (i in 1:length(Plot.vars)) {
    
    spplot(Plot.tmp, Plot.vars[i], main = Plot.vars[i], colorkey=TRUE)
    
    # Manually fitting variogram and kriging
    
    # Data.Var <- variogramformula(paste(Plot.vars[i], "~ 1")), Plot.tmp)
    # plot(Data.Var,pch=20,cex=1.5,col="black",
    #      ylab=expression("Semivariance ("*gamma*")"),
    #      xlab="Distance (m)", main = Plot.vars[i])
    # 
    # sph.model <- vgm(psill=15000, model="Sph", range=1000, nugget=4000)
    # sph.fit <- fit.variogram(object = Data.Var, model = sph.model)
    # # plot(Data.Var,pch=20,cex=1.5,col="black",
    # #     ylab=expression("Semivariance ("*gamma*")"),
    # #     xlab="Distance (m)", main = Plot.vars[i], model=sph.fit)
    # 
    # sph.pred <- krige(formula(paste(Plot.vars[i], "~ 1")), Plot.tmp, Plot.grid, model = sph.fit)
    # spplot(sph.pred, c("var1.pred"))
    # colnames(sph.pred@data)[3:4] <- paste0(Plot.vars[i], c(".pred", ".var"))
    # Plot.grid@data <- cbind(Plot.grid@data, sph.pred@data[,3:4])
    
    
    # Automatically fitting variogram and kriging
    
    # sph.fit <- autofitVariogram(formula(paste(Plot.vars[i], "~ 1")), Plot.tmp, model="Sph")
    # summary(sph.fit)
    # plot(sph.fit)
    
    #you should conver to im file for spatstat
    
    Plot.tmp <- Plot.tmp[is.finite(Plot.tmp[Plot.vars[i]][[1]]),] #remove NA values, if necessary
    
    sph.fit <- autoKrige(formula(paste(Plot.vars[i], "~ 1")), Plot.tmp, Plot.grid)
    plot(sph.fit)
    colnames(sph.fit$krige_output@data)[3:5] <- paste0(Plot.vars[i], c(".pred", ".var", ".stdev"))
    Plot.grid@data <- cbind(Plot.grid@data, sph.fit$krige_output@data[,3:5])
  }
  
  source("~/Documentos/Datos NO publicados/BioIntForest/Analysis/function Raster_extractGrid.R")
  ras.Understory <- Raster_extractGrid(Plot.grid, colnames(Plot.grid@data)[grep(".pred", colnames(Plot.grid@data))], 0.5, plot = F)
  
  source("~/Documentos/Datos NO publicados/BioIntForest/Analysis/function raster.as.im.R")
  
  png(filename = paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/Environ_Maps/Understory_",unique(Data$Plot)[j],".png"), width = 800, height = 1000)
  plot(ras.Understory)
  dev.off()
  
  Vars <- Plot.vars
  
  if (unique(Data$Plot)[j] != "A3") 
    im.Understory <- listof( Fs=raster.as.im(ras.Understory[[1]]), Hed=raster.as.im(ras.Understory[[2]]), Pter=raster.as.im(ras.Understory[[3]]), Rub=raster.as.im(ras.Understory[[4]]), Scl=raster.as.im(ras.Understory[[5]]))
  
  if (unique(Data$Plot)[j] == "A3") 
    im.Understory <-listof( Hed=raster.as.im(ras.Understory[[1]]), Pter=raster.as.im(ras.Understory[[2]]), Rub=raster.as.im(ras.Understory[[3]]), Scl=raster.as.im(ras.Understory[[4]]))
  
  save(ras.Understory, im.Understory, Vars, file = paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/Environ_Maps/Understory_",unique(Data$Plot)[j],".RData"))
  
}

writeRaster(raster, "test_output11", format = "GTiff")