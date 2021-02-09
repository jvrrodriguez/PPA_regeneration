##########################################
####   ANALIZE POINT DATA           ######
##########################################

library(shapefiles)
library(spatstat)
library(spatstat.local)
library(spatstat.utils)
library(forestSAS)
library(raster)
library(rgdal)
library(maptools)

#library(spatialEco) #https://github.com/jeffreyevans/spatialEco
library(shar) #https://r-spatialecology.github.io/shar/
library(ecespa)
library(splines)
library(reshape2)
library(bipartite)


library(ggplot2)
library(gridExtra)
library(dplyr)
library(grid)
library(MASS)
library(cowplot)
library(rasterVis)


#datos con la informacion de posicion de las plantas
setwd("~/Documentos/Datos NO publicados/BioIntForest/Data/Inventarios_floristicos")
data.frond <- data.frame(read.table(file="plantulas_2008_2013.csv", header=T, sep=",",dec=".")) #data.frame(read.table(file="cruce_sombras_plantulas_2008_2011_JVR.csv", sep=";" ,dec=","))
data.pinus <- data.frame(read.table(file="PS09Medidas pino 2009 con coordenadas.csv", header=T, sep=";",dec=","))
data.pinus99 <- data.frame(read.table(file="Medidas pinos 1999.csv", header=T, sep=";",dec=","))

#Treatment plots
Treat <-  c("20%", "30%", "0%", "0%", "30%", "20%", "30%", "20%", "0%") #ONLY 0% and 30%
p.range.max_2013 <- c(NA, 7*2, 8*2, 6*2, 8*2, NA, 4*2, NA, 9*2) #para modificar el tamaño de la parcela (eje X) para el año 2013; he expandido un poco la parcela 9 para que incluya algun Fs adicional (sino el script no funciona)

data.frond$SP2 <- ifelse(data.frond$SP=="Qi" | data.frond$SP=="Qir" | data.frond$SP=="Qic" | data.frond$SP=="Qi?", "Qi",
                         ifelse(data.frond$SP=="Qh" | data.frond$SP=="Qhr", "Qh",
                                ifelse(data.frond$SP=="Fs" | data.frond$SP=="Fsr", "Fs",NA)))
data.frond$SP2<-as.factor(data.frond$SP2)
data.frond$RES <- ifelse(data.frond$SP=="Qir" | data.frond$SP=="Qhr" | data.frond$SP=="Fsr", 1, 0)

#convert to an hyperframe plots
Year <- c(2008, 2009, 2010, 2011, 2012, 2013)
R.rec <- 0.5 #radio del buffer para reclutas
R.under <- 0.5 #radio para understory
R.canopy <- 1 #radio para el dosel
save.output <- TRUE

netmetrics <- data.frame()

for (j in 1:(length(Year)-1)){
  
  cat("Analising year...", Year[j], "\n")
  
  if (save.output == T) pdf(paste0("~/Documentos/Datos NO publicados/BioIntForest/Results/SpatNet_", Year[j],".pdf"), width=10, height=7) #7, 5  
  
  #reset lists for each year
  Plots <- unique(data.frond$PARCELA)

  for (i in Plots) {
    
    cat("Setting up Plot...", i, "\n")
    
    # Setup recruit data --------------------------------------------------------------
    
    data.plot <- subset(data.frond, PARCELA==i & SP2!="<NA>") 
    
    #standardize relative positions only for deciduous trees
    data.plot$MAPX <- data.plot$MAPX - min(data.plot$MAPX)
    data.plot$MAPY <- data.plot$MAPY - min(data.plot$MAPY)
    
    p.range <- c(round(min(data.plot$MAPX)), round(max(data.plot$MAPX)), round(min(data.plot$MAPY)), round(max(data.plot$MAPY)))
    
    if (Year[j] == 2013) p.range[2] <- p.range.max_2013[i]
    
    ### if (j != length(Year)){
      
      data.plot$Size <- data.plot[paste0("ALT",Year[j])][[1]]
      data.sel <- data.plot[which(grepl("ALT", colnames(data.plot)))][, (j+1):length(Year)]
      
      if (!is.null(dim(data.sel))) {
        
        data.plot$Size <- ifelse(rowSums(data.sel== 201.0 | data.sel == -9999.0 | is.na(data.plot$Size), na.rm=T) > 0, 
                                 apply(data.sel, 1, FUN=min, na.rm=T), data.plot$Size)
        
      }
      
      data.plot$Size <- ifelse(is.infinite(data.plot$Size), NA, data.plot$Size)
      data.plot$Size <- ifelse(data.plot$Size == -9999.0, 201, data.plot$Size)
      
      if (!is.null(dim(data.sel))) {
        
        data.plot$Fate <- as.factor(ifelse(is.na(rowSums(data.sel)), 0, 1))
        
      } else {
        
        data.plot$Fate <- as.factor(ifelse(is.na(data.sel), 0, 1))
        
      }
      
      data.year <- data.plot[is.finite(data.plot$Size), ]
      data.year <- data.frame(x = data.year$MAPX, y = data.year$MAPY, SP = as.factor(data.year$SP2), SIZE = data.year$Size / 100, FATE = data.year$Fate, GROWTH = 0)
      data.year <- data.year[!is.na(data.year$SIZE),]
      data.year$SP <- as.factor(as.character(data.year$SP))
      
      data.year <- with(data.year, 
                        ppp(x = x, y = y, marks = data.frame(SP = as.factor(SP), SIZE = log(SIZE + 1), FATE = FATE, GROWTH = GROWTH), 
                            owin(xrange=c(p.range[1],p.range[2]), yrange=c(p.range[3],p.range[4]))))
      
      data_tmp <- data.frame(x = data.year$x, y = data.year$y, species = data.year$marks$SP, size = data.year$marks$SIZE, shape = ifelse(data.year$marks$FATE=="0", 19, 1))
      p <- ggplot(data_tmp, aes(x, y, colour = species, size = size)) + xlab("x coords") + ylab("y coords") + geom_point(shape = data_tmp$shape)
      name <- paste("p",i,sep="_")
      tmp <- list(p)
      #map.Plots[[name]] <- tmp

      layout(matrix(c(1,1,2,3), 2, 2, byrow = FALSE), heights = c(3, 3), widths = c(3, 3))
      par(mar = c(3, 3, 3, 4))
      plot(data.year$x,data.year$y,col=as.integer(data.year$marks$SP), 
           cex = data.year$marks$SIZE, pch = ifelse(data.year$marks$FATE=="0", 19, 1), xlab = "x coords", ylab = "y coords", main = NULL)
      title(paste0("Plot A", i, "; Yr: ", Year[j], "; Th: ",Treat[i], "; n: ", data.year$n, "; Surv : ", round(sum(data.year$marks$FATE == "1") / data.year$n,3)), line = -1, cex.main = 1, outer = TRUE)
      
      Q <- quadratcount(data.year, nx = p.range[2]/2, ny = p.range[4]/2)
      quadrat.test(data.year, nx = p.range[2]/2, ny = p.range[4]/2)
      plot(Q, add = TRUE, cex = 0.5, col="grey")
      
      par(mar = c(3, 0, 3, 3))
      hist(array(Q), xlab="",  main="density plants (2x2 m2)")
      par(mar = c(3, 0, 3, 3))
      hist(data.year$marks$SIZE, xlab="", main="Size plants - log(height+1)")
      layout(matrix(c(1), 1, 1, byrow = FALSE))
      
      # compute neighborhood
      
      SP.rec_neigh <- nnIndex(data.year,id=paste("T",1:data.year$n, sep=""), R = R.rec, smark="SP", buffer=TRUE)
      Rich.rec_year <- as.vector(data.year$n)
      for (k in 1:data.year$n) Rich.rec_year[k] <- length(unique(unlist(SP.rec_neigh$nnSP[k,2:ncol(SP.rec_neigh$nnSP)]))) - 1 #le quito el NA, que es una categoria
      
      l <- apply(SP.rec_neigh$nnSP,1, table)
      #l <- lapply(l, prop.table)  # percentage
      
      library(tidyverse)
      SP.rec_year <-  map_dfr(l, ~ .x %>%
                                 as.list %>%
                                 as_tibble)
      
      SP.rec_year <- as.matrix(SP.rec_year)
      SP.rec_year <- data.frame(ifelse(is.na(SP.rec_year), 0, SP.rec_year))
      
      
      # Setup understory data ---------------------------------------------------

      source("~/Documentos/Datos NO publicados/BioIntForest/Analysis/function raster.as.im.R")
      
      load(paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/Environ_Maps/Understory_A", i, ".RData"))
      ras.Understory <- normalize_stack(ras.Understory)
      
      if (is.null(ras.Understory$Fs_2008.pred)) {
        ras.tmp <- ras.Understory$Hed_2008.pred 
        ras.tmp$v <- ifelse(is.nan(values(ras.tmp)), 0, values(ras.tmp))
        ras.Understory <- c(list(Fs_2008.pred = im.tmp), ras.Understory)
      }
      
      shorter.names <- function(x) {
        
              j <- strsplit(x, "[.]")[[1]]
              k <-substr(j, 1, 1)[[1]] 
              l <- substr(j, 1, 3)[[2]]
              paste0(k,".",l)
              
      }
      
      names(ras.Understory) <- unlist(lapply(names(ras.Understory), function(x, ...) strsplit(x, "_")[[1]][1] ))

      SP.under_tmp <- raster::extract(ras.Understory, method='simple', as.SpatialPointsDataFrame.ppp(data.year), buffer = R.canopy)
      l <- SP.under_tmp
      SP.under_year <- data.frame(matrix(unlist(l), nrow=length(l), byrow=T))
      colnames(SP.under_year) <- names(ras.Understory)
      R.under.max <- max(apply(SP.under_year, 2, max, na.rm=TRUE))
      SP.under.max <- apply(SP.under_year,1,which.max)
      
      over.SP.under_year <- as.vector(SP.under.max)
      for (k in 1:length(SP.under.max)) over.SP.under_year[k] <- SP.under_year[k,SP.under.max[k]]/R.under.max
      
      
      # Setup canopy data ------------------------------------------------------

      load(paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/Environ_Maps/Canopy_A", i, ".RData"))
      ras.Canopy <- ras.Canopy
      ras.canopy_tmp <- lapply (names(ras.Canopy), shorter.names)
      names(ras.Canopy) <- unlist(ras.canopy_tmp)
      
      SP.canopy_tmp <- raster::extract(ras.Canopy, method='simple', as.SpatialPointsDataFrame.ppp(data.year), buffer = R.canopy)
      l <- lapply(SP.canopy_tmp, colSums)
      #l <- lapply(l, prop.table)  # percentage
      #l <- lapply(l, function(x, ...){ ifelse(is.nan(x), 0, x) })
      SP.canopy_year <- data.frame(matrix(unlist(l), nrow=length(l), byrow=T))
      colnames(SP.canopy_year) <- names(ras.Canopy)
      R.canopy.max <- max(apply(SP.canopy_year, 2, max, na.rm=TRUE))
      SP.canopy.max <- apply(SP.canopy_year,1,which.max)
      
      over.SP.canopy_year <- as.vector(SP.canopy.max)
      for (k in 1:length(SP.canopy.max)) over.SP.canopy_year[k] <- SP.canopy_year[k,SP.canopy.max[k]]/R.canopy.max
      
      #SP.max.canopy_year <- ifelse(over.SP.canopy_year == 0, "NA", SP.max.canopy_year)


      # Summarize into a data.frame -----------------------------------------------------

      data.nn_year <- data.frame( lower = SP.rec_neigh$nnSP[1][[1]], 
                            mid = names(ras.Understory)[apply(SP.under_year, 1, which.max)],    
                            higher = names(ras.Canopy)[apply(SP.canopy_year, 1, which.max)], 
                            webID = i, 
                            freq.mid = over.SP.under_year,
                            freq.higher = over.SP.canopy_year,
                            N.rec = applynbd(data.year, R = R.rec, function(Y, ...){npoints(Y)-1}) + 1,
                            Rich.rec = Rich.rec_year,
                            SIZE.rec = applynbd(subset.ppp(data.year, select = "SIZE"), R=R.rec, function(Y, ...) { mean(marks(Y))}),
                            SP1.rec = SP.rec_neigh$nnSP[2][[1]], 
                            Rich_under = rowSums( apply(SP.under_year, 2, function(x, ...){ ifelse(x > median(x, na.rm=T), 1, 0) }) ), 
                            Area_under = rowSums(SP.under_year)/R.under.max,
                            Rich_canopy = rowSums( ifelse(SP.canopy_year > 0, 1, 0) ), 
                            Area_canopy = rowSums(SP.canopy_year)/R.canopy.max  )
      
      data.int_year <- list(lower = SP.rec_year, mid = SP.under_year, higher = SP.canopy_year)
      
      #title(xlab="Service Providers", line=4, cex.lab=1)
      
      # to reorder matrix based on missing interactions
      
      net.Sp_full <- matrix(data = 0, nrow = dim(ras.Understory)[3], ncol = dim(ras.Canopy)[3], dimnames = list(names(ras.Understory), names(ras.Canopy)))
      
      sort.matrix <-  function (x, y) {
        
        if ( nrow(x) != nrow(y) ) {
          
          x <- x[match(rownames(y), rownames(x)), ]
          x <- ifelse(is.na(x), 0, x)
          rownames(x) <- rownames(y)
          
        }
        
        if ( ncol(x) != ncol(y) ) {
          
          x <- x[, match(colnames(y), colnames(x))]
          x <- ifelse(is.na(x), 0, x)
          colnames(x) <- colnames(y)
          
        }
        
        return(x)
        
      }
      
      
      net.Sp_year <- list()
      
      for (l in 1:3) {
        
        net.Sp_year[[l]] <- prop.table( frame2webs(data.frame(data.nn_year, freq = SP.rec_year[colnames(SP.rec_year)[l]][[1]]),  varnames = c("mid", "higher", "webID", "freq"), type.out="array")[,,1])
        net.Sp_year[[l]] <- sort.matrix(net.Sp_year[[l]], net.Sp_full)
        
        tmp <- array(which(colSums(net.Sp_year[[l]]) == 0))
        
        if (length(tmp) > 0) {
          
          #assign a small value to one of the rows whose column is zero
          for (ll in tmp) net.Sp_year[[l]][sample(1:nrow(net.Sp_year[[l]]),1), ll] <- 0.0001
          
        }
        
      }
      
      names(net.Sp_year) <- c("Qh", "Qi","Fs")
      
      net.higher_year <- frame2webs(data.nn_year,  varnames = c("lower", "higher", "webID", "freq.higher"), type.out="array")[,,1]
      
      #https://ibartomeus.github.io/hab-sp_ntw/demo.html
      
      library(RColorBrewer)
      numero_colores <- length(table(net.Sp_year$Qh))
      rgb_cols <- c("#FFFFFFFF", brewer.pal(numero_colores - 1, "YlOrRd"))

      par(mfrow=c(1,3))
      visweb(net.Sp_year$Qh, type="nested")#, square="defined", def.col=rgb_cols)
      visweb(net.Sp_year$Qi, type="nested")#, square="defined")#, def.col=rgb_cols)
      visweb(net.Sp_year$Fs, type="nested")#, square="defined")#, def.col=rgb_cols)
      
      #hist(nnRich.canopy_year)
      
      netmetrics_year <- data.frame(Year = Year[j], Plot = i, Treat = Treat[i], Data = "Obs",
                                    rbind(data.frame(Sp = "Qh", t(networklevel(net.Sp_year[[1]], index=c("connectance", "web asymmetry", "linkage density", "NODF", "ISA", "SA", "interaction evenness")))),
                                    data.frame(Sp = "Qi",t(networklevel(net.Sp_year[[2]], index=c("connectance", "web asymmetry", "linkage density", "NODF", "ISA", "SA", "interaction evenness")))),
                                    data.frame(Sp = "Fs",t(networklevel(net.Sp_year[[3]], index=c("connectance", "web asymmetry", "linkage density", "NODF", "ISA", "SA", "interaction evenness"))))))
      
      netmetrics <- rbind(netmetrics, netmetrics_year)

      
      # Simulate a Thomas process and extract vegetation ------------------------

      
      net.Qh_sim <- list(); net.Qi_sim <- list(); net.Fs_sim <- list()
      
      for (m in 1:5) {  
        
        data.Qh.year <- subset.ppp(subset.ppp(data.year, select = "SP"), marks== "Qh")
        data.Qi.year <- subset.ppp(subset.ppp(data.year, select = "SP"), marks== "Qi")
        data.Fs.year <- subset.ppp(subset.ppp(data.year, select = "SP"), marks== "Fs")
        
        fit.Qh.year <- kppm(unmark(data.Qh.year) ~ 1, "Thomas", statistic="pcf")
        data.Qh.fit <- simulate(fit.Qh.year)
        
        fit.Qi.year <- kppm(unmark(data.Qi.year) ~ 1, "Thomas", statistic="pcf")
        data.Qi.fit <- simulate(fit.Qi.year, 1)
        
        fit.Fs.year <- kppm(unmark(data.Fs.year) ~ 1, "Thomas", statistic="pcf")
        data.Fs.fit <- simulate(fit.Fs.year, 1)
        
        data.fit <- superimpose(Qh = data.Qh.fit[[1]], Qi = data.Qi.fit[[1]], Fs = data.Fs.fit[[1]])
        
        #re-hacer para que sea un data.frame
        data.fit <-  with(data.frame(x = data.fit$x, y = data.fit$y, SP = data.fit$marks), 
             ppp(x = x, y = y, marks = data.frame(SP = as.factor(SP), KK =1), 
                 owin(xrange=c(p.range[1],p.range[2]), yrange=c(p.range[3],p.range[4]))))
   
        SP.rec_sim = data.fit$marks$SP
        
        SP.rec_fit.neigh <- nnIndex(data.fit,id=paste("T",1:data.fit$n, sep=""), smark="SP", R = R.rec, buffer=TRUE)
        Rich.rec_sim <- as.vector(data.fit$n)
        for (k in 1:data.year$n) Rich.rec_sim[k] <- length(unique(unlist(SP.rec_fit.neigh$nnSP[k,2:ncol(SP.rec_fit.neigh$nnSP)]))) - 1 #le quito el NA, que es una categoria
        
        l <- apply(SP.rec_fit.neigh$nnSP,1, table)
        #l <- lapply(l, prop.table)  # percentage
        
        library(tidyverse)
        SP.rec_sim <-  map_dfr(l, ~ .x %>%
                                  as.list %>%
                                  as_tibble)
        SP.rec_sim <- as.matrix(SP.rec_sim)
        SP.rec_sim <- data.frame(ifelse(is.na(SP.rec_sim), 0, SP.rec_sim))
        
        ##
        
        SP.under_sim <- raster::extract(ras.Understory, method='simple', as.SpatialPoints.ppp(data.fit), buffer = R.canopy)
        l <- SP.under_sim
        SP.under_sim <- data.frame(matrix(unlist(l), nrow=length(l), byrow=T))
        colnames(SP.under_sim) <- names(ras.Understory)
        R.under.max <- max(apply(SP.under_sim, 2, max, na.rm=TRUE))
        SP.under.max <- apply(SP.under_sim,1,which.max)
        
        over.SP.under_sim <- as.vector(SP.under.max)
        for (k in 1:length(SP.under.max)) over.SP.under_sim[k] <- SP.under_sim[k,SP.under.max[k]]/R.under.max
        
        SP.under_sim <- names(ras.Understory)[apply(SP.under_sim, 1, which.max)]
        
        ##
        
        SP.canopy_sim <- raster::extract(ras.Canopy, method='simple', as.SpatialPoints.ppp(data.fit), buffer = R.canopy)
        
        l <- lapply(SP.canopy_sim, colSums)
        SP.canopy_sim <- data.frame(matrix(unlist(l), nrow=length(l), byrow=T))
        colnames(SP.canopy_sim) <- names(ras.Canopy)
        R.canopy.max <- max(apply(SP.canopy_sim, 2, max, na.rm=TRUE))
        SP.canopy.max <- apply(SP.canopy_sim,1,which.max)
        
        over.SP.canopy_sim <- as.vector(SP.canopy.max)
        for (k in 1:length(SP.canopy.max)) over.SP.canopy_sim[k] <- SP.canopy_sim[k,SP.canopy.max[k]]/R.canopy.max
        
        SP.canopy_sim <- names(ras.Canopy)[apply(SP.canopy_sim, 1, which.max)]
        
        ##
        
        data.nn_sim <- data.frame(lower = SP.rec_sim, mid = SP.under_sim, higher = SP.canopy_sim, webID = i)
        
        #par(mfrow=c(2,1), mar=c(8,2,2,1)+.1)
        #barplot(sort(table(paste(data.nn_year$higher, data.nn_year$mid, data.nn_year$lower, sep="-")), decreasing = T), las=2, border = 1, cex.lab=1, cex.axis=1, font=1, col.axis="black", main="Observed")
        #barplot(sort(table(paste(data.nn_year$higher, data.nn_sim$mid, data.fit$marks$SP, sep="-")), decreasing = T), las=2, border = 1, cex.lab=1, cex.axis=1, font=1, col.axis="black", main="Simulated")
        
        net.Qh_sim[[m]] <- prop.table( frame2webs(data.nn_sim,  varnames = c("mid", "higher", "webID", "lower.Qh"), type.out="array")[,,1])
        net.Qi_sim[[m]] <- prop.table( frame2webs(data.nn_sim,  varnames = c("mid", "higher", "webID", "lower.Qi"), type.out="array")[,,1])
        net.Fs_sim[[m]] <- prop.table( frame2webs(data.nn_sim,  varnames = c("mid", "higher", "webID", "lower.Fs"), type.out="array")[,,1])
        
        net.Qh_sim[[m]] <- sort.matrix(net.Qh_sim[[m]], net.Sp_full)
        net.Qi_sim[[m]] <- sort.matrix(net.Qi_sim[[m]], net.Sp_full)
        net.Fs_sim[[m]] <- sort.matrix(net.Fs_sim[[m]], net.Sp_full)
        
      }
      
      net.Sp_sim <- list()
      
      for (l in 1:3) {
        
        x <- get(paste0("net.",colnames(SP.rec_year)[l],"_sim"))
        net.Sp_sim[[l]] <- apply(array(unlist(x), c(dim(x[[1]]), length(x))), 1:2, mean, na.rm = TRUE)
        colnames(net.Sp_sim[[l]]) <- colnames(net.Sp_full)
        rownames(net.Sp_sim[[l]]) <- rownames(net.Sp_full)
        
        tmp <- array(which(colSums(net.Sp_sim[[l]]) == 0))
        
        if (length(tmp) > 0) {
          
          #assign a small value to one of the rows whose column is zero
          for (ll in tmp) net.Sp_sim[[l]][sample(1:nrow(net.Sp_sim[[l]]),1), ll] <- 0.0001
            
        }
        
      }
      
      names(net.Sp_sim) <- c("Qh", "Qi","Fs")
      
      
      par(mfrow=c(1,3))
      visweb(net.Sp_sim$Qh, type="nested")#, square="defined", def.col=rgb_cols)
      visweb(net.Sp_sim$Qi, type="nested")#, square="defined", def.col=rgb_cols)
      visweb(net.Sp_sim$Fs, type="nested")#, square="defined", def.col=rgb_cols)
      
      net.Qh_diff <- net.Sp_year$Qh/net.Sp_sim$Qh
      net.Qi_diff <- net.Sp_year$Qi/net.Sp_sim$Qi
      net.Fs_diff <- net.Sp_year$Fs/net.Sp_sim$Fs
      
      visweb(net.Qh_diff, type="nested", NA.col="grey")#, square="defined", def.col=rgb_cols)
      visweb(net.Qi_diff, type="nested", NA.col="grey")#, square="defined", def.col=rgb_cols)
      visweb(net.Fs_diff, type="nested", NA.col="grey")#, square="defined", def.col=rgb_cols)
      
      netmetrics_sim <- data.frame(Year = Year[j], Plot = i, Treat = Treat[i], Data = "Sim",
                                    rbind(data.frame(Sp = "Qh", t(networklevel(net.Sp_sim[[1]], index=c("connectance", "web asymmetry", "linkage density", "NODF", "ISA", "SA", "interaction evenness")))),
                                          data.frame(Sp = "Qi",t(networklevel(net.Sp_sim[[2]], index=c("connectance", "web asymmetry", "linkage density", "NODF", "ISA", "SA", "interaction evenness")))),
                                          data.frame(Sp = "Fs",t(networklevel(net.Sp_sim[[3]], index=c("connectance", "web asymmetry", "linkage density", "NODF", "ISA", "SA", "interaction evenness"))))))
      
      netmetrics <- rbind(netmetrics, netmetrics_sim)
    
  }
  
  dev.off()
  
  if (save.output == T) write.csv(netmetrics, paste0("~/Documentos/Datos NO publicados/BioIntForest/Results/netmetrics.csv"))
  

}

library(Rmisc)

index=c("connectance", "web.asymmetry", "NODF", "interaction.strength.asymmetry", "specialisation.asymmetry", "linkage.density", "interaction.evenness")

i=6

data.sum <- summarySE(netmetrics, measurevar=index[i], groupvars=c("Year","Data"))
colnames(data.sum)[4] <- "var"
data.sum$Year <- factor(data.sum$Year)
  
ggplot(data.sum, aes(x=Year, y=var, fill=Data)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=var-se, ymax=var+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))
