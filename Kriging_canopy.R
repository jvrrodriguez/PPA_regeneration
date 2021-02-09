# Setup data --------------------------------------------------------------

setwd("~/Documentos/Datos NO publicados/BioIntForest/Data/Inventarios_floristicos")

data.frond <- data.frame(read.table(file="cruce_sombras_plantulas_2008_2011_JVR.csv", header=T, sep=";",dec=","))


# Now krig mixed vs pure canopy -------------------------------------------

for (i in unique(data.frond$PARCELA)) {

  data.plot <- subset(data.frond, PARCELA==i & Dosel!="<NA>")
  
  #stantardize relative positions only for decidious trees
  data.plot$MAPX <- data.plot$MAPX - min(data.plot$MAPX)
  data.plot$MAPY <- data.plot$MAPY - min(data.plot$MAPY)
  p.range <- c(round(min(data.plot$MAPX)), round(max(data.plot$MAPX)), round(min(data.plot$MAPY)), round(max(data.plot$MAPY)))
  
  coordinates(data.plot) <- c("MAPX", "MAPY")
  data.plot@bbox[,2] <- c(30,40)
  data.plot$Dosel <- as.numeric(ifelse(data.plot$Dosel == "mixto", 1, 0))
  
  Plot.grid <- as(raster(ncol = 15*2, nrow = 20*2, crs=NULL, ext = extent(data.plot@bbox)), "SpatialPolygons")
  Plot.grid <- SpatialPolygonsDataFrame(Plot.grid, data.frame(n=1:length(Plot.grid)))
  
  sph.fit <- autoKrige.cv(Dosel~ 1, data.plot, Plot.grid, nfold = 10) #, nmax = 1 for binary data
  #sph.fit <- idw(Dosel~ 1, data.plot, Plot.grid,)

  Plot.grid@data <- data.frame(Plot.grid@data, krig.pred = sph.fit$krige_output@data[,3])
  
  source("~/Documentos/Datos NO publicados/BioIntForest/Analysis/function Raster_extractGrid.R")
  ras.Canopy <- Raster_extractGrid(Plot.grid, "krig.pred", 0.5, plot = F)
  png(filename = paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/Environ_Maps/Canopy_A",i,".png"), width = 800, height = 1000)
  plot(ras.Canopy)
  plot(data.plot, col=as.integer(data.plot$Dosel), add=T, pch = 19, xlab = "x coords", ylab = "y coords")
  dev.off()
  
  source("~/Documentos/Datos NO publicados/BioIntForest/Analysis/function raster.as.im.R")
  im.Canopy <-listof( Canopy=raster.as.im(ras.Canopy[[1]]))
  
  save(ras.Canopy, im.Canopy, file = paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/Environ_Maps/Canopy_A",i,".RData"))
  
}