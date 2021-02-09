
library(reshape2)
library(ggplot2)
library(spatstat)
library(raster)

Year <- c(2008, 2009, 2010, 2011)
Plot <- c("A2", "A3", "A4", "A5", "A7", "A9")

data.plot <- hyperframe(Year = Year, A2 = listof(rep(NA, 4)), A3 = listof(rep(NA, 4)), A4 = listof(rep(NA, 4)), A5 = listof(rep(NA, 4)), A7 = listof(rep(NA, 4)), A9 = listof(rep(NA, 4)))

for (j in 1:length(Year)){

  load(file=paste0("~/Documentos/Datos NO publicados/BioIntForest/Results/Results.PPA_", Year[j],".RData"))
  
  data.sp <- data.sp.reg$pred.dens
  
  data.plot$A2[[j]] <- data.sp[[2]][[2]][[1]]
  data.plot$A3[[j]] <- data.sp[[3]][[1]]
  data.plot$A4[[j]] <- data.sp[[4]][[1]]
  data.plot$A5[[j]] <- data.sp[[5]][[1]]
  data.plot$A7[[j]] <- data.sp[[7]][[1]]
  data.plot$A9[[j]] <- data.sp[[9]][[1]]
  
}

#https://stackoverflow.com/questions/35202728/standard-deviation-function-throws-error-when-na-rm-true-while-using-calc-in-r
pop.sd=function(x, na.rm){
  v = var(x,na.rm = na.rm)
  l = if(na.rm) sum(!is.na(x)) else length(x)
  sqrt(v*(l-1)/l)
}

raster.summary <- stack()

for (j in 1:length(Plot)){
  
  #plot(data.plot[,j+1,1])
  raster.tmp <- stack(raster(data.plot[,j+1,1][[1]]), raster(data.plot[,j+1,1][[2]]), raster(data.plot[,j+1,1][[3]]), raster(data.plot[,j+1,1][[4]]))
  mean <- mean(raster.tmp)
  sd <- raster::calc(raster.tmp, pop.sd, na.rm=T)
  raster.tmp2 <- stack(mean, sd)
  names(raster.tmp2) <- c(paste(Plot[j],'mean'), paste(Plot[j],'sd'))
  raster.summary <- stack(raster.summary, raster.tmp2) 
  
}

plot(raster.summary)
