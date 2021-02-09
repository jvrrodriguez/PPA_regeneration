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
  
  Data <- data.frame(read.table(file="~/Documentos/Datos NO publicados/BioIntForest/Data/Fotos_hemisfericas/Hemiphot_Plots_2004.csv", fileEncoding="utf-8", row.names = NULL, header=T, sep=",",dec="."))
  loc.tmp <- substring(as.character(Data$File), 3, 14)
  loc.tmp2 <- unlist(strsplit(loc.tmp, split='.JPG', fixed=TRUE))
  Data <- data.frame(Data, X = substr(loc.tmp2, 1, 1), Y = as.integer(gsub("[^0-9\\.]", "", loc.tmp2)), Exp = ifelse(grepl("luz",loc.tmp2),1,0))
  Data$X <- ifelse(as.character(Data$X) =="ñ", "o", as.character(Data$X))
  Data$CuadCode <- paste0(toupper(Data$X), sprintf("%02d",Data$Y))
  
  # unique(Data$Plot)
  # unique(Data$plot)
  # unique(Data$CuadCode)
  
  coords = matrix(c(0,0, 0,2*20, 2*14.99,2*20, 2*14.99,0), # antes 14.5 pero haciendo analisis falta un cacho 
                  ncol = 2, byrow = TRUE) #excepto la A1 y A3, que van horizontales
  P1 = Polygon(coords)
  Plot = SpatialPolygons(list(Polygons(list(P1), ID = "a")))
  # plot(Plot, axes = TRUE)
  
  CuadCode <- paste0(as.vector(t(replicate(15, LETTERS[1:15]))), sprintf("%02d",rep(1:20, 15)))
  Points <- SpatialPointsDataFrame(makegrid(Plot, cellsize=c(2,2)), data.frame(CuadCode))
  
  # No se porque tengo que cambiarlos manualmente
  Plot@bbox[,1] <- c(0,0)
  Plot@bbox[,2] <- c(30,40)
  Points@bbox[,1] <- c(0,0)
  Points@bbox[,2] <- c(30,40) 
  
  for (j in 1:length(unique(Data$Plot))) {
    
    Points.data <- subset(Data, Plot == unique(Data$Plot)[j])
    
    if (sum(table(Points.data$CuadCode) > 1) > 0 ) {  # Si hay repetidos
      sel.repeat <- names(which(table(Points.data$CuadCode) > 1))
      
      #quito los repetidos y sobreexpuestos
      Points.data <- Points.data[!(Points.data$CuadCode == sel.repeat & Points.data$Exp == 1),]
      #Points.data <- aggregate(Points.data, list(Data$X, Data$Y), mean)
    }
    
    #Points$CuadCode <- as.character(Points$CuadCode)
    #dim(Points.data)[1] == dim(Points)[1]
    #Points.data <- cbind(Data, CuadCode)
    
    Plot.tmp <- merge(Points, Points.data, by.x = "CuadCode", by.y = "CuadCode")
    #Plot.vars <- apply(expand.grid(unique(Data$Sp), unique(Data$Year)), 1, paste, collapse="_")
    Plot.vars <- c("CanOpen", "LAI", "DirectBelow", "DiffBelow", "DirectBelow.Yr", "DiffBelow.Yr", "N.Sunflecks", "Mdn.Sunflecks", "Max.Sunflecks") #colnames(Plot.tmp@data)
    ### IMPORTANT. You should to modify in the last lines depending on the number and name of variables
    
    spplot(Plot.tmp, Plot.vars)
    
    # This is the new grid for predicting
    Plot.grid <- as(raster(ncol = 15*2, nrow = 20*2, crs=NULL, ext = extent(Points@bbox)), "SpatialPolygons") #*4
    Plot.grid <- SpatialPolygonsDataFrame(Plot.grid, data.frame(n=1:length(Plot.grid)))
    
    #https://rstudio-pubs-static.s3.amazonaws.com/80464_9156596afb2e4dcda53e3650a68df82a.html
    #https://rpubs.com/jguelat/autocorr2
    
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
      
      Plot.tmp[is.infinite(Plot.tmp[Plot.vars[i]][[1]])] <- NA #Remove infinite values if necessary
      Plot.tmp <- Plot.tmp[is.finite(Plot.tmp[Plot.vars[i]][[1]]),] #remove NA values, if necessary
      
      if (j != 6) sph.fit <- autoKrige(formula(paste(Plot.vars[i], "~ 1")), Plot.tmp, Plot.grid)
      if (j == 6) sph.fit <- autoKrige(formula(paste(Plot.vars[i], "~ 1")), Plot.tmp, Plot.grid, fix.values = c(0,NA,NA))
      plot(sph.fit)
      
      colnames(sph.fit$krige_output@data)[3:5] <- paste0(Plot.vars[i], c(".pred", ".var", ".stdev"))
      Plot.grid@data <- cbind(Plot.grid@data, sph.fit$krige_output@data[,3:5])
    }
    
    source("~/Documentos/Datos NO publicados/BioIntForest/Analysis/function Raster_extractGrid.R")
    ras.Light <- Raster_extractGrid(Plot.grid, colnames(Plot.grid@data)[grep(".pred", colnames(Plot.grid@data))], 0.5, plot = F)
    
    source("~/Documentos/Datos NO publicados/BioIntForest/Analysis/function raster.as.im.R")
    
    png(filename = paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/Environ_Maps/Light_",unique(Data$Plot)[j],".png"), width = 800, height = 1000)
    plot(ras.Light)
    dev.off()
    
    im.Light <- listof(CanOpen=raster.as.im(ras.Light[[1]]), LAI=raster.as.im(ras.Light[[2]]), 
                      DirectBelow=raster.as.im(ras.Light[[3]]), DiffBelow=raster.as.im(ras.Light[[4]]),
                      DirectBelow.Yr=raster.as.im(ras.Light[[5]]), DiffBelow.Yr=raster.as.im(ras.Light[[6]]),
                      N.Sunflecks=raster.as.im(ras.Light[[7]]), Mdn.Sunflecks=raster.as.im(ras.Light[[8]]), Max.Sunflecks=raster.as.im(ras.Light[[9]]))
    Vars <- Plot.vars
    save(ras.Light, im.Light, Vars, file = paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/Environ_Maps/Light_",unique(Data$Plot)[j],".RData"))
    
}
  