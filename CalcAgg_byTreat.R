##########################################
####   ANALIZE POINT DATA           ######
##########################################

library(raster)
library(shapefiles)

library(spatstat.core)
library(spatstat.local)
library(spatstat.utils)
#library(spatialEco) #https://github.com/jeffreyevans/spatialEco
library(shar) #https://r-spatialecology.github.io/shar/
library(ecespa)
library(onpoint)

library(MASS)
library(zoo)
library(splines)

library(ggplot2)
library(gridExtra)
library(dplyr)
library(reshape2)
library(grid)
library(cowplot)
library(corrplot)
library(rasterVis)

#Incluye funciones para resumir la información
source("~/Documentos/Datos NO publicados/BioIntForest/PPA_regeneration/FunsAgg.R")

#datos con la informacion de posicion de las plantas
setwd("~/Documentos/Datos NO publicados/BioIntForest/Data/Inventarios_floristicos")
data.frond <- data.frame(read.table(file = "plantulas_2008_2013.csv", header = T, sep = ",",dec = ".")) #data.frame(read.table(file="cruce_sombras_plantulas_2008_2011_JVR.csv", sep=";" ,dec=","))
data.pinus <- data.frame(read.table(file = "PS09Medidas pino 2009 con coordenadas.csv", header = T, sep = ";",dec = ","))
data.pinus99 <- data.frame(read.table(file = "Medidas pinos 1999.csv", header = T, sep = ";",dec = ","))

#Treatment plots
Treat <-  c("20%", "30%", "0%", "0%", "30%", "20%", "30%", "20%", "0%") #ONLY 0% and 30%
p.range.max_2013 <- c(NA, 7*2, 8*2, 6*2, 8*2, NA, 4*2, NA, 9*2) #para modificar el tamaño de la parcela (eje X) para el año 2013; he expandido un poco la parcela 9 para que incluya algun Fs adicional (sino el script no funciona)

data.frond$SP2 <- ifelse(data.frond$SP == "Qi" | data.frond$SP == "Qir" | data.frond$SP == "Qic" | data.frond$SP == "Qi?", "Qi",
                         ifelse(data.frond$SP == "Qh" | data.frond$SP == "Qhr", "Qh",
                                ifelse(data.frond$SP == "Fs" | data.frond$SP == "Fsr", "Fs",NA)))
data.frond$SP2 <- as.factor(data.frond$SP2)
data.frond$RES <- ifelse(data.frond$SP == "Qir" | data.frond$SP == "Qhr" | data.frond$SP == "Fsr", 1, 0)

#convert to an hyperframe plots
sizejuv <- 0.2 # umbral (m) para discriminar juvenil de plantula
sizead <- 2 # umbral (m) para discriminar adulto de juvenil
Year <- c(2008, 2009, 2010, 2011, 2012, 2013)[1:4]
nsim <- 199
fit.gam = FALSE
save.output <- TRUE
gof.win <- 3 # number of consecutive values out of the confidence intervals

summary.year <- data.frame(); summary.plot <- summary.year; 
summary.size.plot <- summary.year; summary.size <- summary.year; summary.rec.Q <- summary.year; summary.ad.Q <- summary.year


for (j in 1:length(Year)) {
  
  cat("Analising year...", Year[j], "\n")
  
  if (save.output == T) pdf(paste0("~/Documentos/Datos NO publicados/BioIntForest/PPA_Results/Reports/PPA_", Year[j],".pdf"), width = 10, height = 7) #7, 5  
  
  #reset lists for each year
  Plots <- unique(data.frond$PARCELA)
  map.Plots <- list()
  
  data.sp <- hyperframe(Plot = 1:max(Plots), Year = Year[j], Treat = Treat, ppp = listof(rep(NA, 9)), ppp = listof(rep(NA, 9)), covs = listof(rep(NA, 9)))
  data.size <- data.sp
  data.size.c <- data.sp
  data.sp.rec <- data.sp
  data.sp.rec.Fs <- data.sp
  data.sp.rec.Qh <- data.sp
  data.sp.rec.Qi <- data.sp
  data.size.rec <- data.sp
  data.size.c.rec <- data.sp
  data.size.c.ad <- data.sp
  data.fate <- data.sp
  data.fate.rec <- data.sp
  data.fate.rec.i <- data.sp
  data.fate.ad <- data.sp
  data.growth <- data.sp
  data.growth.ad <- data.sp
  
  dens.all <- hyperframe(Plot = 1:max(Plots), Year = Year[j], Treat = Treat, im = list(rep(NA, 9)))
  dens.sp <- dens.all
  dens.sp.rec <- dens.all
  dens.sp.rec.rich <- dens.all
  dens.sp.rec.shan <- dens.all
  dens.size <- dens.all
  dens.size.c <- dens.all
  dens.size.c.ad <- dens.all
  dens.size.c.rec <- dens.all
  
  data.sp.list <- list()
  
  matrix.env <- list() #para almacenar las matrices de correlaciones de las variables
  matrix.dens <- list()
  
  for (i in Plots) {
    
    cat("Setting up Plot...", i, "\n")
    
    data.plot <- subset(data.frond, PARCELA == i & SP2 != "<NA>") 
    
    #stantardize relative positions only for decidious trees
    data.plot$MAPX <- data.plot$MAPX - min(data.plot$MAPX)
    data.plot$MAPY <- data.plot$MAPY - min(data.plot$MAPY)
    
    p.range <- c(round(min(data.plot$MAPX)), round(max(data.plot$MAPX)), round(min(data.plot$MAPY)), round(max(data.plot$MAPY)))
    
    if (Year[j] == 2013) p.range[2] <- p.range.max_2013[i]
    
    data.pinus99.plot <- subset(data.pinus99, Parcela == paste0("A",i) & is.finite(DM))
    data.pinus99.plot$Size <- rowMeans(cbind(data.pinus99.plot$P.MeanHt, data.pinus99.plot$A.MeanHt), na.rm = T)
    x <- data.pinus99.plot$DM[!is.nan(data.pinus99.plot$Size)]
    y <- data.pinus99.plot$Size[!is.nan(data.pinus99.plot$Size)]
    
    # Modelo para calcular las relacion alometrica entre dbh y altura
    m <- nls(y~a*x/(b+x)) 
    cor(y, predict(m))
    
    data.pinus.plot <- subset(data.pinus, PARCELA == paste0("A",i))
    data.pinus.plot <- data.pinus.plot[!is.na(data.pinus.plot$Dn_NS_07) | !is.na(data.pinus.plot$Dn_EO_07),]
    
    #Load environmental maps
    
    source("~/Documentos/Datos NO publicados/BioIntForest/PPA_regeneration/function raster.as.im.R")
    
    #load(paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/Environ_Maps/MDS-MDT_A", i, ".RData"))

    load(paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/Environ_Maps/Light_A", i, ".RData"))
    im.Light <- normalize_stack(im.Light)
    im.Light$N.Sunflecks$v <- ifelse(is.nan(im.Light$N.Sunflecks$v), 1,  im.Light$N.Sunflecks$v)
    
    load(paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/Environ_Maps/Understory_A", i, ".RData"))
    im.Understory <- normalize_stack(im.Understory)
    # ras.Understory <- normalize_stack(stack(sum(normalize_stack(ras.Understory))))
    # names(ras.Understory) <- "Understory"
    # im.Understory <- listof( Canopy=raster.as.im(ras.Understory))
    
    if (is.null(im.Understory$Fs)) {
      im.tmp <- im.Understory$Hed 
      im.tmp$v <- ifelse(is.nan(im.tmp$v), 1,  im.tmp$v)
      im.Understory <- c(list(Fs = im.tmp), im.Understory)
    }
    
    if (names(table(is.nan(im.Understory$Rub$v))) == TRUE)  {
      im.Understory$Rub$v <- ifelse(is.nan(im.Understory$Rub$v), 1,  im.Understory$Rub$v)
    }

    load(paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/Environ_Maps/Canopy_A", i, ".RData"))
    ras.Canopy <- stack(aggregate(ras.Canopy.join[[2]], fact = 5, na.rm = TRUE))
    names(ras.Canopy) <- "Canopy"
    im.Canopy <- normalize_stack(im.Canopy)
    
    #load(paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/Environ_Maps/Lad_A", i, ".RData"))
    
    ras.all <- stack((ras.Light), (ras.Understory), (ras.Canopy)) #ras.MDTS, ras.lad
    matrix.env[[i]] <- layerStats(ras.all,'pearson', na.rm = T)$`pearson correlation coefficient`
    matrix.env[[i]] <- ifelse(is.finite(matrix.env[[i]]), matrix.env[[i]], NA)
    
    corrplot(matrix.env[[i]])
    
    # Setup point data --------------------------------------------------------------
    
    if (Year[j] < 2009) {
      
      data.pinus.year <- data.frame(x = data.pinus.plot$X, y = data.pinus.plot$Y,
                                    SP = "Ps", 
                                    SIZE = predict(m, newdata = data.frame(x = rowMeans(cbind(data.pinus.plot$Dn_NS_07, data.pinus.plot$Dn_EO_07), na.rm = T))),
                                    FATE = ifelse(is.na(data.pinus.plot$Muerto.07) | data.pinus.plot$Muerto.07 == 0, 1, 0))
      
    } else {
      
      data.pinus.plot.tmp <- data.pinus.plot[is.na(data.pinus.plot$Muerto.07) | is.na(data.pinus.plot$Muerto.09),]
      data.pinus.year <- data.frame(x = data.pinus.plot.tmp$X, y = data.pinus.plot.tmp$Y,
                                    SP = "Ps", 
                                    SIZE = predict(m, newdata = data.frame(x = data.pinus.plot.tmp$Dn_09)),
                                    FATE = ifelse(is.na(data.pinus.plot.tmp$Muerto.09) | data.pinus.plot.tmp$Th_09 == 0, 1, 0))
    }
    
    if (Year[j] == 2009) {
      
      data.pinus.year <- data.frame(data.pinus.year, 
                                    GROWTH =  round(data.pinus.plot.tmp$Dn_09 - rowMeans(cbind(data.pinus.plot.tmp$Dn_NS_07, data.pinus.plot.tmp$Dn_EO_07), na.rm = T), 1)) 
      data.pinus.year$GROWTH <- log(ifelse(data.pinus.year$GROWTH < 0, 0, data.pinus.year$GROWTH) + 1)
      
    } else {
      
      data.pinus.year <- data.frame(data.pinus.year, GROWTH =  0)
      
    }
    
    if (j != length(Year)) {
      
      data.plot$Size <- data.plot[paste0("ALT",Year[j])][[1]]
      data.sel <- data.plot[which(grepl("ALT", colnames(data.plot)))][, (j + 1):length(Year)]
      
      if (!is.null(dim(data.sel))) {
        
        data.plot$Size <- ifelse(rowSums(data.sel == 201.0 | data.sel == -9999.0 | is.na(data.plot$Size), na.rm = T) > 0, 
                                 apply(data.sel, 1, FUN = min, na.rm = T), data.plot$Size)
        
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
      data.year <- rbind(data.year, data.pinus.year)
      data.year <- data.year[!is.na(data.year$SIZE),]
      data.year$SP <- as.factor(as.character(data.year$SP))
      
      data.year <- with(data.year, 
                        ppp(x = x, y = y, marks = data.frame(SP = as.factor(SP), SIZE = log(SIZE + 1), FATE = FATE, GROWTH = GROWTH), 
                            owin(xrange = c(p.range[1],p.range[2]), yrange = c(p.range[3],p.range[4]))))
      
      data_tmp <- data.frame(x = data.year$x, y = data.year$y, species = data.year$marks$SP, size = data.year$marks$SIZE, shape = ifelse(data.year$marks$FATE == "0", 19, 1))
      p <- ggplot(data_tmp, aes(x, y, colour = species, size = size)) + xlab("x coords") + ylab("y coords") + geom_point(shape = data_tmp$shape)
      name <- paste("p",i,sep = "_")
      tmp <- list(p)
      map.Plots[[name]] <- tmp
      
      
      layout(matrix(c(1,1,2,3), 2, 2, byrow = FALSE), heights = c(3, 3), widths = c(3, 3))
      par(mar = c(3, 3, 3, 4))
      plot(data.year$x,data.year$y,col = as.integer(data.year$marks$SP), 
           cex = data.year$marks$SIZE, pch = ifelse(data.year$marks$FATE == "0", 19, 1), xlab = "x coords", ylab = "y coords", main = NULL)
      title(paste0("Plot A", i, "; Yr: ", Year[j], "; Th: ",Treat[i], "; n: ", data.year$n, "; Surv : ", round(sum(data.year$marks$FATE == "1") / data.year$n,3)), line = -1, cex.main = 1, outer = TRUE)
      
    } else {
      
      data.plot$Size <- data.plot[paste0("ALT",Year[j])][[1]]
      data.plot$Size <- ifelse(is.infinite(data.plot$Size), NA, data.plot$Size)
      data.plot$Size <- ifelse(data.plot$Size == -9999.0, 201, data.plot$Size)
      
      data.year <- data.plot[is.finite(data.plot$Size), ]
      data.year <- data.frame(x = data.year$MAPX, y = data.year$MAPY, SP = as.factor(data.year$SP2), SIZE = data.year$Size / 100, GROWTH = 0)
      data.year <- rbind(data.year, data.pinus.year[c(1:4,6)])
      data.year$SP <- as.factor(as.character(data.year$SP))
      
      data.year <- with(data.year, 
                        ppp(x = x, y = y, marks = data.frame(SP = as.factor(SP), SIZE = log(SIZE + 1), GROWTH = GROWTH), 
                            owin(xrange = c(p.range[1],p.range[2]), yrange = c(p.range[3],p.range[4]))))
      
      data_tmp <- data.frame(x = data.year$x, y = data.year$y, species = data.year$marks$SP, size = data.year$marks$SIZE, shape = 21)
      p <- ggplot(data_tmp, aes(x, y, colour = species, size = size)) + 
        xlab("x coords") + ylab("y coords") + geom_point(shape = data_tmp$shape) + scale_size(range = c(0,1))
      name <- paste("p",i,sep = "_")
      tmp <- list(p)
      map.Plots[[name]] <- tmp
      
      layout(matrix(c(1,1,2,3), 2, 2, byrow = FALSE), heights = c(2, 2), widths = c(2, 2))
      par(mar = c(3, 3, 3, 4))
      plot(data.year$x,data.year$y,col = as.integer(data.year$marks$SP), 
           cex = data.year$marks$SIZE, xlab = "x coords", ylab = "y coords", main = NULL)
      title(paste0("Plot A", i, "; Yr: ", Year[j], "; Th: ",Treat[i], "; n: ", data.year$n), line = -1, cex.main = 1, outer = TRUE)
      
    }
    
    Q <- quadratcount(data.year, nx = p.range[2]/2, ny = p.range[4]/2)
    quadrat.test(data.year, nx = p.range[2]/2, ny = p.range[4]/2)
    plot(Q, add = TRUE, cex = 0.5, col = "grey")
    
    par(mar = c(3, 0, 3, 3))
    hist(array(Q), xlab = "density plants",  main = "")
    par(mar = c(3, 0, 3, 3))
    hist(data.year$marks$SIZE, xlab = "log(height+1)", main = "")
    abline(v = log(sizejuv + 1), col = "grey"); abline(v = log(sizead + 1), col = "grey")
    layout(matrix(c(1), 1, 1, byrow = FALSE))
    
    
    # Arrange datasets to analyse
    
    dens.all$im[[i]] <- density(data.year, sigma = bw.diggle(data.year))
    
    data.sp$ppp[[i]] <- cut.ppp(data.year, "SP")
    
    dens.sp$im[[i]] <- density(split(data.sp$ppp[[i]]), sigma = bw.diggle(data.sp$ppp[[i]]))
    
    data.size$ppp[[i]] <- cut.ppp(data.year, "SIZE", breaks = c(-0.1, log(sizejuv + 1), log(sizead + 1), Inf), labels = c("Recruit", "Recruit", "Adult"))
    
    data.size$ppp[[i]]$marks <- as.factor(as.character(data.size$ppp[[i]]$marks))
    plot(split(data.size$ppp[[i]]), main = paste0("Plot A", i, "; Year = ", Year[j]))
    
    dens.size$im[[i]] <- density(split(data.size$ppp[[i]]), sigma = bw.diggle(data.size$ppp[[i]]))
    
    data.size.c$ppp[[i]] <- subset.ppp(data.year, select = "SIZE")
    
    dens.size.c$im[[i]] <- density(data.size.c$ppp[[i]], sigma = bw.diggle(data.size.c$ppp[[i]]))
    
    data.sp.rec$ppp[[i]] <- data.sp$ppp[[i]][which(data.size$ppp[[i]]$marks != "Adult"),]
    data.sp.rec.Fs$ppp[[i]] <- data.sp.rec$ppp[[i]][which(data.sp.rec$ppp[[i]]$marks == "Fs"),]
    data.sp.rec.Qh$ppp[[i]] <- data.sp.rec$ppp[[i]][which(data.sp.rec$ppp[[i]]$marks == "Qh"),]
    data.sp.rec.Qi$ppp[[i]] <- data.sp.rec$ppp[[i]][which(data.sp.rec$ppp[[i]]$marks == "Qi"),]
    
    data.size.rec$ppp[[i]] <- data.size$ppp[[i]][which(data.size$ppp[[i]]$marks != "Adult"),]
    data.size.c.rec$ppp[[i]] <- data.size.c$ppp[[i]][which(data.size$ppp[[i]]$marks != "Adult"),]
    
    if (sum(is.na(data.size$ppp[[i]]$marks)) > 0) {
      
      data.size$ppp[[i]] <- data.size$ppp[[i]][!is.na(data.size$ppp[[i]]$marks)]
      
    }
    
    data.size.c$ppp[[i]] <- data.size.c$ppp[[i]][!is.na(as.numeric(data.size.c$ppp[[i]]$marks)),]
    
    data.size.c.ad$ppp[[i]] <- data.size.c$ppp[[i]][which(data.size$ppp[[i]]$marks == "Adult"),]
    dens.size.c.ad$im[[i]] <- density(data.size.c.ad$ppp[[i]], sigma = bw.diggle(data.size.c.ad$ppp[[i]]))
    dens.size.c.rec$im[[i]] <- density(data.size.c.rec$ppp[[i]], sigma = bw.diggle(data.size.c.rec$ppp[[i]]))
    dens.sp.rec$im[[i]] <- density(split(data.sp.rec$ppp[[i]]), sigma = bw.diggle(data.sp.rec$ppp[[i]]))
    
    dens.sp.rec.rich$im[[i]] <- dens.size.c.rec$im[[i]]
    dens.sp.rec.rich$im[[i]]$v <- im.diversity(dens.sp.rec$im[[i]]$Fs$v , dens.sp.rec$im[[i]]$Qh$v , dens.sp.rec$im[[i]]$Qi$v)[[1]]
    
    dens.sp.rec.shan$im[[i]] <- dens.size.c.rec$im[[i]]
    dens.sp.rec.shan$im[[i]]$v <- im.diversity(dens.sp.rec$im[[i]]$Fs$v , dens.sp.rec$im[[i]]$Qh$v , dens.sp.rec$im[[i]]$Qi$v)[[2]]
    
    data.sp.list[[i]] <- 1:data.sp.rec$ppp[[i]]$n
    
    if (sum(table(data.sp.rec$ppp[[i]]$marks) < 6 & table(data.sp.rec$ppp[[i]]$marks) > 0) > 0) {
      
      rm.sp.rec <- names(which(table(data.sp.rec$ppp[[i]]$marks) < 6))
      
      #Add new points (up to 10) randomly and avoid removing species
      for (k in 1:length(rm.sp.rec)) {
        
        data.sp.tmp <- data.sp.rec$ppp[[i]][which(data.sp.rec$ppp[[i]]$marks == rm.sp.rec[k]),]
        data.sp.new <- runifpoint(10 - data.sp.tmp$n, win = dens.all$im[[i]])
        marks(data.sp.new) <- as.factor(rep(rm.sp.rec[k], data.sp.new$n))
        levels(data.sp.new$marks) <- levels(data.sp.tmp$marks)
        data.sp.rec$ppp[[i]] <- superimpose(data.sp.rec$ppp[[i]], data.sp.new)
        
      }
      
    }
    
    if (sum(unique(data.sp.rec$ppp[[i]]$marks) == "Ps")) {
      
      data.sp.list[[i]] <- which(data.sp.rec$ppp[[i]]$marks != "Ps")
      data.sp.rec$ppp[[i]] <- data.sp.rec$ppp[[i]][data.sp.list[[i]],]
      
    }
    
    data.sp.rec$ppp[[i]]$marks <- as.factor(as.character(data.sp.rec$ppp[[i]]$marks))
    data.size.rec$ppp[[i]]$marks <- as.factor(as.character(data.size.rec$ppp[[i]]$marks))
    
    if (j != length(Year)) {
      
      data.fate$ppp[[i]] <- cut.ppp(data.year, "FATE")
      data.fate.rec$ppp[[i]] <- data.fate$ppp[[i]][which(data.size$ppp[[i]]$marks != "Adult"),]
      data.fate.rec.i$ppp[[i]] <- data.fate$ppp[[i]]
      data.fate.rec.i$ppp[[i]]$marks <- ifelse(as.character(data.size$ppp[[i]]$marks) == "Adult", "Tree",
                                               ifelse(as.character(data.size$ppp[[i]]$marks) != "Adult" & data.fate$ppp[[i]]$marks == 0, "RcDead",
                                                      ifelse(as.character(data.size$ppp[[i]]$marks) != "Adult" & data.fate$ppp[[i]]$marks == 1, "RcSurv", NA)))
      
      if (sum(is.na(data.fate.rec.i$ppp[[i]]$marks)) > 0) {
        
        data.fate.rec.i$ppp[[i]] <- data.fate.rec.i$ppp[[i]][-which(is.na(data.fate.rec.i$ppp[[i]]$marks))]
        
      } 
      
      plot(split(data.fate.rec$ppp[[i]]), main = paste0("Plot A", i, "; Year = ", Year[j]))
      
    }
    
    if (Year[j] == 2009) {
      
      data.growth$ppp[[i]] <- data.year
      data.growth$ppp[[i]]$marks <- data.growth$ppp[[i]]$marks[4][[1]]
      data.growth.ad$ppp[[i]] <- data.growth$ppp[[i]][data.size$ppp[[i]]$marks == "Adult" & data.sp$ppp[[i]]$marks == "Ps",]
      
      data.fate.ad$ppp[[i]] <- data.fate$ppp[[i]][data.size$ppp[[i]]$marks == "Adult" & data.sp$ppp[[i]]$marks == "Ps",]
      
    }
    
    
    ras.all <- stack(raster(dens.all$im[[i]]), raster(dens.size$im[[i]]["Adult"][[1]]), raster(dens.size$im[[i]]["Recruit"][[1]]), raster(dens.size.c$im[[i]]), raster(dens.size.c.rec$im[[i]]),
          raster(dens.sp$im[[i]]["Fs"][[1]]), raster(dens.sp$im[[i]]["Qh"][[1]]), raster(dens.sp$im[[i]]["Qi"][[1]]), raster(dens.sp.rec.rich$im[[i]]), raster(dens.sp.rec.shan$im[[i]]))
    
    names(ras.all) <- c("Dens", "Dens.adult", "Dens.rec", "Dens.size", "Dens.size.rec", "Dens.Fs.rec", "Dens.Qh.rec", "Dens.Qi.rec","Dens.rich.rec", "Dens.shan.rec")   

    matrix.dens[[i]] <- layerStats(ras.all,'pearson', na.rm = T)$`pearson correlation coefficient`
    matrix.dens[[i]] <- ifelse(is.finite(matrix.dens[[i]]), matrix.dens[[i]], NA)
    
    
    # Setup covariates
    
    data.sp$covs[[i]] <- listof(#MDT = im.MDTS[[1]], MDS = im.MDTS[[2]],
                                CanOpen = im.Light[[1]], LAI = im.Light[[2]], DirectBelow = im.Light[[3]], DiffBelow = im.Light[[4]], DirectBelow.Yr = im.Light[[5]], DiffBelow.Yr = im.Light[[6]], N.Sunflecks = im.Light[[7]], Mdn.Sunflecks = im.Light[[8]], Max.Sunflecks = im.Light[[9]],
                                Canopy = im.Canopy[[1]], 
                                Fs = im.Understory[[1]], Hed = im.Understory[[2]], Pter = im.Understory[[3]], Rub = im.Understory[[4]], Scl = im.Understory[[5]],
                                Dens = dens.all$im[[i]], Dens.adult = dens.size$im[[i]]["Adult"][[1]], Dens.rec = dens.size$im[[i]]["Recruit"][[1]], Dens.size = dens.size.c$im[[i]], Dens.size.rec = dens.size.c.rec$im[[i]],
                                Dens.Fs.rec = dens.sp$im[[i]]["Fs"][[1]], Dens.Qh.rec = dens.sp$im[[i]]["Qh"][[1]], Dens.Qi.rec = dens.sp$im[[i]]["Qi"][[1]], Dens.rich.rec = dens.sp.rec.rich$im[[i]], Dens.shan.rec = dens.sp.rec.shan$im[[i]])
    
    #Summarize data sets
    
    Q.sp <- quadratcount(split(data.sp.rec$ppp[[i]]), nx = p.range[2]/2, ny = p.range[4]/2)
    
    Q.size.rec <- quadratcount(split(data.size.rec$ppp[[i]]), nx = p.range[2]/2, ny = p.range[4]/2)
    
    Q.size.ad <- quadratcount(data.size.c.ad$ppp[[i]], nx = p.range[2]/2, ny = p.range[4]/2)
    
    if (j != length(Year)) {
      
      table.plot <- data.frame(sp = data.sp$ppp[[i]]$marks, size = data.size$ppp[[i]]$marks, fate = data.fate$ppp[[i]]$marks)
      
      #table.plot <- table.plot[rowSums(is.na(table.plot)) != 1,] #delete row with missing information (sp or size)
      
      summary.plot <- rbind(summary.plot, data.frame(Year = Year[j], Plot = i, Treat = Treat[i], 
                                                     table.plot %>% group_by(sp, size, fate) %>% count()))
      
      if (length(unique(data.fate.rec$ppp[[i]]$marks)) > 1) {
        
        Q.fate <- quadratcount(split(data.fate.rec$ppp[[i]]), nx = p.range[2]/2, ny = p.range[4]/2)
        
        summary.rec.Q <- rbind(summary.rec.Q, data.frame(Year = Year[j], Plot = i, Treat = Treat[i], Sp.Fs = array(Q.sp$Fs), Sp.Qh = array(Q.sp$Qh), Sp.Qi = array(Q.sp$Qi), Size = array(Q.size.rec$Recruit), Fate = array(Q.fate$`1`)))
        
        summary.ad.Q <- rbind(summary.ad.Q, data.frame(Year = Year[j], Plot = i, Treat = Treat[i], Size = array(Q.size.ad)))
        
      } else {
        
        Q.fate <- quadratcount(data.fate.rec$ppp[[i]], nx = p.range[2]/2, ny = p.range[4]/2)
        
        summary.rec.Q <- rbind(summary.rec.Q, data.frame(Year = Year[j], Plot = i, Treat = Treat[i], Sp.Fs = array(Q.sp$Fs), Sp.Qh = array(Q.sp$Qh), Sp.Qi = array(Q.sp$Qi), Size = array(Q.size.rec$Recruit), Fate = array(Q.fate)))
        
        summary.ad.Q <- rbind(summary.ad.Q, data.frame(Year = Year[j], Plot = i, Treat = Treat[i], Size = array(Q.size.ad)))
        
      }
      
      
      
    } else {
      
      table.plot <- data.frame(sp = data.sp$ppp[[i]]$marks, size = ifelse(data.year$marks$SIZE > -0.1 & data.year$marks$SIZE < log(sizejuv + 1), "Recruit", "Adult"))
      
      #table.plot <- table.plot[rowSums(is.na(table.plot)) != 1,] #delete row with missing information (sp or size)
      
      summary.plot <- rbind(summary.plot, data.frame(Year = Year[j], Plot = i, Treat = Treat[i], 
                                                     table.plot %>% group_by(sp, size) %>% count()))
      
      summary.rec.Q <- rbind(summary.rec.Q, data.frame(Year = Year[j], Plot = i, Treat = Treat[i], Sp.Fs = array(Q.sp$Fs), Sp.Qh = array(Q.sp$Qh), Sp.Qi = array(Q.sp$Qi), Size = array(Q.size.rec$Recruit)))
      
      summary.ad.Q <- rbind(summary.ad.Q, data.frame(Year = Year[j], Plot = i, Treat = Treat[i], Size = array(Q.size.ad)))
      
    }
    
    summary.size.tmp <- hist(data.year$marks$SIZE, xlab = "log(height+1)", breaks = seq(0, 2.5, 0.125), plot = F, main = "")
    summary.size.plot <- rbind(summary.size, data.frame(Year = Year[j], Plot = i, Treat = Treat[i], breaks = summary.size.tmp$breaks[-1], Freq = summary.size.tmp$counts))
    
  }
  
  if (save.output == T) {
    
    write.csv(summary.ad.Q, paste0("~/Documentos/Datos NO publicados/BioIntForest/PPA_Results/Reports/summary.ad.Q_", Year[j], ".csv"))
    write.csv(summary.rec.Q, paste0("~/Documentos/Datos NO publicados/BioIntForest/PPA_Results/Reports/summary.recruit.Q_", Year[j], ".csv"))
    write.csv(summary.plot, paste0("~/Documentos/Datos NO publicados/BioIntForest/PPA_Results/Reports/summary.Plot_", Year[j], ".csv"))
    
  }
 
  
  # Generate summary plots for supplementary material
  
  prow <- plot_grid(map.Plots$p_3[[1]] + theme(legend.position = "none"), map.Plots$p_4[[1]] + theme(legend.position = "none"), map.Plots$p_9[[1]] + theme(legend.position = "none"), 
                    map.Plots$p_2[[1]] + theme(legend.position = "none"), map.Plots$p_5[[1]] + theme(legend.position = "none"), map.Plots$p_7[[1]] + theme(legend.position = "none"),
                    align = 'vh',
                    labels = c("a)", "b)", "c)", "d)", "e)", "f)"),
                    hjust = -1)
  
  # extract the legend from one of the plots
  legend <- get_legend(
    # create some space to the left of the legend
    map.Plots$p_2[[1]] + theme(legend.box.margin = margin(0, 0, 0, 12))
  )
  
  # add the legend to the row we made earlier. Give it one-third of 
  # the width of one plot (via rel_widths).
  print(
    plot_grid(prow, legend, rel_widths = c(3, .4)))
  
  var.env <- c("Canopy", "CanOpen", "LAI", "DiffBelow", "N.Sunflecks", "Max.Sunflecks", "Fs", "Hed", "Rub", "Pter", "Scl")
  var.dens <- c("Dens.adult", "Dens.rec", "Dens.rich.rec", "Dens.shan.rec")
  var.all <- c(var.env, var.dens)
  
  for (i in 1:length(var.all)) {
    
    raster.tmp <- stack(raster(data.sp$covs$`3`[[var.all[i]]]), raster(data.sp$covs$`4`[[var.all[i]]]), raster(data.sp$covs$`9`[[var.all[i]]]), raster(data.sp$covs$`2`[[var.all[i]]]), raster(data.sp$covs$`5`[[var.all[i]]]), raster(data.sp$covs$`7`[[var.all[i]]]))
    names(raster.tmp) <- c(paste0("Plot_", c(3,4,9,2,5,7)))
    print( 
      gplot(raster.tmp) + 
        geom_tile(aes(fill = value)) + guides(fill = guide_legend(title = var.all[i])) +
        facet_wrap(~ variable) +
        scale_fill_gradientn(colours = rev(terrain.colors(225))) +
        coord_equal() + theme_classic() 
    )
    
  }
  
  # Test only recruits ----------------------------------------------------
  
  cat("Test only recruits patterns \n")

  L.E.rec <- list(); g.E.rec <- L.E.rec; kNN.E.rec <- L.E.rec; F.E.rec <- L.E.rec
  L.E.rec.cat <- NULL; g.E.rec.cat <- NULL; kNN.E.rec.cat <- NULL; F.E.rec.cat <- NULL
  g.E.rec.gof <- NULL
  
  for (i in Plots) {
    
    g.E.rec[[i]] <- envelope(unmark(data.sp.rec$ppp[[i]]), pcf, r = seq(0,4,0.02), nsim = nsim, fix.n = TRUE, correction = "Ripley", savefuns = TRUE, savepatterns = TRUE)
    kNN.E.rec[[i]] <- envelope(unmark(data.sp.rec$ppp[[i]]), Gest, r = seq(0,1,0.01), nsim = nsim, fix.n = TRUE, correction = "rs", savefuns = TRUE, savepatterns = TRUE)
  
    g.E.rec.cat <- cbind(g.E.rec[[i]]$r, g.E.rec.cat, ppp.cat(g.E.rec[[i]]))
    kNN.E.rec.cat <- cbind(kNN.E.rec[[i]]$r, kNN.E.rec.cat, ppp.cat(kNN.E.rec[[i]]))
      
    g.E.rec.test <- gof.int(g.E.rec[[i]], gof.win)

    g.E.rec.gof <- rbind(g.E.rec.gof, 
                         data.frame(Plot = i, g.E.rec.test))
    
  }
  
  g.E.rec.ctrl <- pool(g.E.rec[[3]], g.E.rec[[4]], g.E.rec[[9]], savefuns = TRUE)
  g.E.rec.thnn <- pool(g.E.rec[[2]], g.E.rec[[5]], g.E.rec[[7]], savefuns = TRUE)
  
  kNN.E.rec.ctrl <- pool(kNN.E.rec[[3]], kNN.E.rec[[4]], kNN.E.rec[[9]])
  kNN.E.rec.thnn <- pool(kNN.E.rec[[2]], kNN.E.rec[[5]], kNN.E.rec[[7]])
  
  g.E.rec.ctrl.test <- gof.int(g.E.rec.ctrl, gof.win)
  g.E.rec.thnn.test <- gof.int(g.E.rec.thnn, gof.win)

  g.E.rec.gof <- rbind(g.E.rec.gof, 
                       data.frame(Plot = "Ctrl", g.E.rec.ctrl.test),
                       data.frame(Plot = "Thnn", g.E.rec.thnn.test))
  
  
  par(mfrow = c(2,2), mar = c(1, 1, 1.25, 1.25), oma = c(4, 4, 2, 2)) 
  #plot(L.E.rec.ctrl, . -r ~ r, shade=c("hi", "lo"), legend = F, main = NULL)
  plot(g.E.rec.ctrl, main = NULL, legend = F)
  plot(kNN.E.rec.ctrl, main = NULL, legend = F)
  plot(g.E.rec.thnn, main = NULL, legend = F)
  plot(kNN.E.rec.thnn, main = NULL, legend = F)
  title("Ctrl (u) vs Thinning (l) plots", line = 0, outer = TRUE)
  old.par <- par(mfrow = c(1,1), mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
  par(old.par)
  
  
  
  # Test Cluster (Thomas) process -------------------------------------------
  
  fit.clust <- list()
  sum.clust <- data.frame(R = double(), mu = double() )
  g.E.t <- list()
  g.E.t.cat <- NULL; g.E.t.gof <- NULL
  
  par(mfrow = c(2,3), mar = c(0.5, 0.5, 0.2, 0.2), oma = c(4, 4, 0.2, 0.2))
  
  for (i in Plots) {
    
    fit.clust[[i]] <- kppm(unmark(data.sp.rec$ppp[[i]]) ~ 1, "Thomas", statistic = "pcf", data = data.sp$covs[[i]], savefuns = TRUE, savepatterns = TRUE)
    sum.clust[i,1] <- fit.clust[[i]]$modelpar[2]
    sum.clust[i,2] <- fit.clust[[i]]$mu
    
    g.E.t[[i]] <- envelope(fit.clust[[i]], pcf, nsim = nsim, fix.n = TRUE, correction = "Ripley", savefuns = TRUE, savepatterns = TRUE)
    
    g.E.t.cat <- cbind(g.E.t[[i]]$r, g.E.t.cat, ppp.cat(g.E.t[[i]]))
    g.E.t.test <- gof.int(g.E.t[[i]], gof.win)

    g.E.t.gof <- rbind(g.E.t.gof, 
                         data.frame(Plot = i, g.E.t.test))
    
  }
  
  par(old.par)
  
  g.E.t.ctrl <- pool(g.E.t[[3]], g.E.t[[4]], g.E.t[[9]], savefuns = TRUE)
  g.E.t.thnn <- pool(g.E.t[[2]], g.E.t[[5]], g.E.t[[7]], savefuns = TRUE)
  
  g.E.t.ctrl.test <- gof.int(g.E.t.ctrl, gof.win)
  g.E.t.thnn.test <- gof.int(g.E.t.thnn, gof.win)

  g.E.t.gof <- rbind(g.E.t.gof, 
                       data.frame(Plot = "Ctrl", g.E.t.ctrl.test),
                       data.frame(Plot = "Thnn", g.E.t.thnn.test))
  
  
  par(mfrow = c(2,2), mar = c(1, 1, 1.25, 1.25), oma = c(4, 4, 2, 2)) 
  plot(g.E.t.ctrl, main = NULL, legend = F)
  plot(g.E.t.thnn, main = NULL, legend = F)
  title("Ctrl (u) vs Thinning (l) plots", line = 0, outer = TRUE)
  old.par <- par(mfrow = c(1,1), mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
  par(old.par)
  
  
  
  # Test patterns of marks --------------------------------------------
  
  # Between species...
  # Al final esta marca se estudiará en detalle en otro paper
  
  # cat("Test differences between species \n")
  # 
  # L.E.sp <- list(); g.E.sp <- L.E.sp; kNN.E.sp <- L.E.sp; Jdif.E.sp <- L.E.sp; markcon.E.sp <- L.E.sp
  # L.E.sp.cat <- NULL; g.E.sp.cat <- NULL; kNN.E.sp.cat <- NULL; Jdif.E.sp.cat <- NULL; markcon.E.sp.cat <- NULL
  # 
  # for (i in Plots) {
  #   
  #   #L.E.sp[[i]] <- alltypes(data.sp$ppp[[i]], Ldot, r = seq(0,10,0.05), nsim = nsim, envelope = TRUE, correction="Ripley", title = NULL, savefuns=TRUE, savepatterns = TRUE)
  #   g.E.sp[[i]] <- alltypes(data.sp$ppp[[i]], pcfdot, r = seq(0,4,0.02), nsim = nsim, envelope = TRUE, correction="Ripley", title = NULL, savefuns=TRUE, savepatterns = TRUE)
  #   kNN.E.sp[[i]] <- alltypes(data.sp$ppp[[i]], Gdot, r = seq(0,1,0.01), nsim = nsim, envelope = TRUE, correction="rs", title = NULL, savefuns=TRUE, savepatterns = TRUE)
  #   markcon.E.sp[[i]] <- alltypes(data.sp$ppp[[i]], markconnect, nsim = nsim, envelope = TRUE, title = NULL, savefuns=TRUE, savepatterns = TRUE, simulate=expression(rlabel(data.sp$ppp[[i]])))
  #   Jdif.E.sp[[i]] <- alltypes(data.sp$ppp[[i]], Jdif, r = seq(0,1,0.01), nsim = nsim, savefuns=TRUE, savepatterns = TRUE, envelope = TRUE, simulate = expression(rlabel(data.sp$ppp[[i]])))
  #   
  #   #L.E.sp.tmp <- do.call(cbind, lapply(L.E.sp[[i]]$fns, ppp.cat))
  #   #colnames(L.E.sp.tmp) <- apply(melt(t(L.E.sp[[i]]$which))[-3], 1, paste, collapse="-")
  #   #L.E.sp.cat <- cbind(L.E.sp[[i]]$fns[[1]]$r, L.E.sp.cat, L.E.sp.tmp)
  #   
  #   g.E.sp.tmp <- do.call(cbind, lapply(g.E.sp[[i]]$fns, ppp.cat))
  #   colnames(g.E.sp.tmp) <- apply(melt(t(g.E.sp[[i]]$which))[-3], 1, paste, collapse="-")
  #   g.E.sp.cat <- cbind(g.E.sp[[i]]$fns[[1]]$r, g.E.sp.cat, g.E.sp.tmp)
  #   
  #   kNN.E.sp.tmp <- do.call(cbind, lapply(kNN.E.sp[[i]]$fns, ppp.cat))
  #   colnames(kNN.E.sp.tmp) <- apply(melt(t(kNN.E.sp[[i]]$which))[-3], 1, paste, collapse="-")
  #   kNN.E.sp.cat <- cbind(kNN.E.sp[[i]]$fns[[1]]$r, kNN.E.sp.cat, kNN.E.sp.tmp)
  #   
  #   markcon.E.sp.tmp <- do.call(cbind, lapply(markcon.E.sp[[i]]$fns, ppp.cat))
  #   colnames(markcon.E.sp.tmp) <- apply(melt(t(markcon.E.sp[[i]]$which))[-3], 1, paste, collapse="-")
  #   markcon.E.sp.cat <- cbind(markcon.E.sp[[i]]$fns[[1]]$r, markcon.E.sp.cat, markcon.E.sp.tmp)
  #   
  #   Jdif.E.sp.tmp <- do.call(cbind, lapply(Jdif.E.sp[[i]]$fns, ppp.cat))
  #   colnames(Jdif.E.sp.tmp) <- apply(melt(t(Jdif.E.sp[[i]]$which))[-3], 1, paste, collapse="-")
  #   Jdif.E.sp.cat <- cbind(Jdif.E.sp[[i]]$fns[[1]]$r, Jdif.E.sp.cat,Jdif.E.sp.tmp)
  #   
  # }
  # 
  # #L.E.sp.ctrl <- pool(L.E.sp[[3]], L.E.sp[[4]], L.E.sp[[9]]); L.E.sp.thnn <- pool(L.E.sp[[2]], L.E.sp[[5]], L.E.sp[[7]])
  # g.E.sp.ctrl <- pool(g.E.sp[[3]], g.E.sp[[4]], g.E.sp[[9]]); g.E.sp.thnn <- pool(g.E.sp[[2]], g.E.sp[[5]], g.E.sp[[7]])
  # kNN.E.sp.ctrl <- pool(kNN.E.sp[[3]], kNN.E.sp[[4]], kNN.E.sp[[9]]); kNN.E.sp.thnn <- pool(kNN.E.sp[[2]], kNN.E.sp[[5]], kNN.E.sp[[7]])
  # Jdif.E.sp.ctrl <- pool(Jdif.E.sp[[3]], Jdif.E.sp[[4]], Jdif.E.sp[[9]]); Jdif.E.sp.thnn <- pool(Jdif.E.sp[[2]], Jdif.E.sp[[5]], Jdif.E.sp[[7]])
  # if (j != length(Year)) { markcon.E.sp.ctrl <- pool(markcon.E.sp[[3]], markcon.E.sp[[4]], markcon.E.sp[[9]]); markcon.E.sp.thnn <- pool(markcon.E.sp[[2]], markcon.E.sp[[5]], markcon.E.sp[[7]])} 
  # 
  # #plot(L.E.sp.ctrl, . -r ~ r, shade=c("hi", "lo"), legend = F, mar.panel=c(1, 1, 1, 1), title = NULL)
  # plot(g.E.sp.ctrl, legend = F, mar.panel=c(1, 1, 1, 1), title = NULL)
  # plot(kNN.E.sp.ctrl, legend = F, mar.panel=c(1, 1, 1, 1), title = NULL)
  # plot(Jdif.E.sp.ctrl, legend = F, mar.panel=c(1, 1, 1, 1), title = NULL)
  # if (j != length(Year)) { plot(markcon.E.sp.ctrl, legend = F, mar.panel=c(1, 1, 1, 1), title = NULL) }
  # 
  # #plot(L.E.sp.thnn, . -r ~ r, shade=c("hi", "lo"), legend = F, mar.panel=c(1, 1, 1, 1), title = NULL) }
  # plot(g.E.sp.thnn, legend = F, mar.panel=c(1, 1, 1, 1), title = NULL)
  # plot(kNN.E.sp.thnn, legend = F, mar.panel=c(1, 1, 1, 1), title = NULL)
  # plot(Jdif.E.sp.thnn, legend = F, mar.panel=c(1, 1, 1, 1), title = NULL)
  # if (j != length(Year)) { plot(markcon.E.sp.thnn, legend = F, mar.panel=c(1, 1, 1, 1), title = NULL) }
  
  
  
  # Betwenn categorical sizes ...
  
  cat("Test differences between categorical sizes \n")
  
  L.E.size <- list(); g.E.size <- L.E.size; kNN.E.size <- L.E.size; g.E.ac <- L.E.size; Jdif.E.size <- L.E.size; markcon.E.size <- L.E.size
  L.E.size.cat <- NULL; g.E.size.cat <- NULL; kNN.E.size.cat <- NULL; g.E.ac.cat <- NULL; Jdif.E.size.cat <- NULL; markcon.E.size.cat <- NULL
  g.E.size.gof <- NULL; g.E.ac.gof <- NULL
  
  for (i in Plots) {
    
    g.E.size[[i]] <- alltypes(data.size$ppp[[i]], pcfdot, r = seq(0,4,0.02), nsim = nsim, envelope = TRUE, correction = "Ripley", title = NULL, savefuns = TRUE, savepatterns = TRUE)
    kNN.E.size[[i]] <- alltypes(data.size$ppp[[i]], Gdot, r = seq(0,1,0.01), nsim = nsim, envelope = TRUE, correction = "rs", title = NULL, savefuns = TRUE, savepatterns = TRUE)
    markcon.E.size[[i]] <- alltypes(data.size$ppp[[i]], markconnect, nsim = nsim, envelope = TRUE, title = NULL, savefuns = TRUE, savepatterns = TRUE, simulate = expression(rlabel(data.size$ppp[[i]])))
    
    #https://r-spatialecology.github.io/onpoint/reference/simulate_antecedent_conditions.html#ref-examples
    null_model_AC <- simulate_antecedent_conditions(data.size$ppp[[i]], i = "Recruit", j = "Adult", nsim = nsim, heterogenous = FALSE, savefuns = TRUE, savepatterns = TRUE)
    g.E.ac[[i]] <- envelope(Y = data.size$ppp[[i]], fun = pcf, nsim = nsim, correction = "Ripley", simulate = null_model_AC, envelope = TRUE, savefuns = TRUE, savepatterns = TRUE)

    ## Para calcular random labelling pero nos sugirieron que no es adecuado como modelo nulo
    ### Jdif.E.size[[i]] <- envelope(data.size$ppp[[i]], Jdif, r = seq(0,1,0.01), nsim = nsim, i = "Recruit", correction="rs", savefuns=TRUE, savepatterns = TRUE, simulate = expression(rlabel(data.size$ppp[[i]])))
    
    g.E.size.tmp <- do.call(cbind, lapply(g.E.size[[i]]$fns, ppp.cat))
    colnames(g.E.size.tmp) <- apply(melt(t(g.E.size[[i]]$which))[-3], 1, paste, collapse = "-")
    g.E.size.cat <- cbind(g.E.size[[i]]$fns[[1]]$r, g.E.size.cat, g.E.size.tmp)
    
    g.E.size.test1 <- gof.int(g.E.size[[i]]$fns[[1]], gof.win)
    g.E.size.test2 <- gof.int(g.E.size[[i]]$fns[[2]], gof.win)

    g.E.size.gof <- rbind(g.E.size.gof, 
                          data.frame(Plot = i, Class = dimnames(g.E.size[[i]][[2]])[[1]][1], g.E.size.test1),
                          data.frame(Plot = i, Class = dimnames(g.E.size[[i]][[2]])[[1]][2], g.E.size.test2))
    
    kNN.E.size.tmp <- do.call(cbind, lapply(kNN.E.size[[i]]$fns, ppp.cat))
    colnames(kNN.E.size.tmp) <- apply(melt(t(kNN.E.size[[i]]$which))[-3], 1, paste, collapse = "-")
    kNN.E.size.cat <- cbind(kNN.E.size[[i]]$fns[[1]]$r, kNN.E.size.cat, kNN.E.size.tmp)
    
    markcon.E.size.tmp <- do.call(cbind, lapply(markcon.E.size[[i]]$fns, ppp.cat))
    colnames(markcon.E.size.tmp) <- apply(melt(t(markcon.E.size[[i]]$which))[-3], 1, paste, collapse = "-")
    markcon.E.size.cat <- cbind(markcon.E.size[[i]]$fns[[1]]$r, markcon.E.size.cat, markcon.E.size.tmp)
    
    g.E.ac.cat <- cbind(g.E.ac[[i]]$r, g.E.ac.cat, ppp.cat(g.E.ac[[i]]))
    g.E.ac.test <- gof.int(g.E.ac[[i]], gof.win)
 
    g.E.ac.gof <- rbind(g.E.ac.gof, 
                        data.frame(Plot = i, g.E.ac.test))

  }
  
  g.E.size.ctrl <- pool(g.E.size[[3]], g.E.size[[4]], g.E.size[[9]])
  g.E.size.ctrl1 <- pool(g.E.size[[3]][[1]][[1]], g.E.size[[4]][[1]][[1]], g.E.size[[9]][[1]][[1]], savefuns = TRUE)
  g.E.size.ctrl2 <- pool(g.E.size[[3]][[1]][[2]], g.E.size[[4]][[1]][[2]], g.E.size[[9]][[1]][[2]], savefuns = TRUE)
  
  g.E.size.thnn <- pool(g.E.size[[2]], g.E.size[[5]], g.E.size[[7]])
  g.E.size.thnn1 <- pool(g.E.size[[2]][[1]][[1]], g.E.size[[5]][[1]][[1]], g.E.size[[7]][[1]][[1]], savefuns = TRUE)
  g.E.size.thnn2 <- pool(g.E.size[[2]][[1]][[2]], g.E.size[[5]][[1]][[2]], g.E.size[[7]][[1]][[2]], savefuns = TRUE)
  
  kNN.E.size.ctrl <- pool(kNN.E.size[[3]], kNN.E.size[[4]], kNN.E.size[[9]])
  kNN.E.size.thnn <- pool(kNN.E.size[[2]], kNN.E.size[[5]], kNN.E.size[[7]])
  
  g.E.ac.ctrl <- pool(g.E.ac[[3]], g.E.ac[[4]], g.E.ac[[9]], savefuns = TRUE)
  g.E.ac.thnn <- pool(g.E.ac[[2]], g.E.ac[[5]], g.E.ac[[7]], savefuns = TRUE)
  
  g.E.size.ctrl.test1 <- gof.int(g.E.size.ctrl1, gof.win)
  g.E.size.ctrl.test2 <- gof.int(g.E.size.ctrl2, gof.win)
  g.E.size.thnn.test1 <- gof.int(g.E.size.thnn1, gof.win)
  g.E.size.thnn.test2 <- gof.int(g.E.size.thnn2, gof.win)
  
  g.E.size.gof <- rbind(g.E.size.gof, 
                        data.frame(Plot = "Ctrl", Class = dimnames(g.E.size[[i]][[2]])[[1]][1], g.E.size.ctrl.test1),
                        data.frame(Plot = "Ctrl", Class = dimnames(g.E.size[[i]][[2]])[[1]][2], g.E.size.ctrl.test2),
                        data.frame(Plot = "Thnn", Class = dimnames(g.E.size[[i]][[2]])[[1]][1], g.E.size.thnn.test1),
                        data.frame(Plot = "Thnn", Class = dimnames(g.E.size[[i]][[2]])[[1]][2], g.E.size.thnn.test2))
  
  g.E.ac.ctrl.test <- gof.int(g.E.ac.ctrl, gof.win)
  g.E.ac.thnn.test <- gof.int(g.E.ac.thnn, gof.win)
  
  g.E.ac.gof <- rbind(g.E.ac.gof, 
                      data.frame(Plot = "Ctrl", g.E.ac.ctrl.test),
                      data.frame(Plot = "Thnn", g.E.ac.thnn.test))

  
  if (Year[j] < 2012) { markcon.E.size.ctrl <- pool(markcon.E.size[[3]], markcon.E.size[[4]], markcon.E.size[[9]]); markcon.E.size.thnn <- pool(markcon.E.size[[2]], markcon.E.size[[5]], markcon.E.size[[7]])}
  
  #plot(L.E.size.ctrl, . -r ~ r, shade=c("hi", "lo"), legend = F, mar.panel=c(1, 1, 1, 1), title = NULL)
  plot(g.E.size.ctrl, legend = F, mar.panel = c(1, 1, 1, 1), title = NULL)
  plot(kNN.E.size.ctrl, legend = F, mar.panel = c(1, 1, 1, 1), title = NULL)
  plot(g.E.ac.ctrl, legend = F, mar.panel = c(1, 1, 1, 1), title = NULL)
  #plot(Jdif.E.size.ctrl, legend = F, mar.panel = c(1, 1, 1, 1), title = NULL)
  if (Year[j] < 2012) { plot(markcon.E.size.ctrl, legend = F, mar.panel = c(1, 1, 1, 1), title = NULL) }
  
  #plot(L.E.size.thnn, . -r ~ r, shade=c("hi", "lo"), legend = F, mar.panel=c(1, 1, 1, 1), title = NULL)
  plot(g.E.size.thnn, legend = F, mar.panel = c(1, 1, 1, 1), title = NULL)
  plot(kNN.E.size.thnn, legend = F, mar.panel = c(1, 1, 1, 1), title = NULL)
  plot(g.E.ac.thnn, legend = F, mar.panel = c(1, 1, 1, 1), title = NULL)
  #plot(Jdif.E.size.thnn, legend = F, mar.panel = c(1, 1, 1, 1), title = NULL)
  if (Year[j] < 2012) { plot(markcon.E.size.thnn, legend = F, mar.panel = c(1, 1, 1, 1), title = NULL) }
  
  
  
  # Betwenn continuous sizes ...
  
  cat("Test differences in continuous sizes \n")
  
  markcor.E.size.c.ad <- list(); markcor.E.size.c.rec <- markcor.E.size.c.ad
  markcor.E.size.c.ad.cat <- NULL; markcor.E.size.c.rec.cat <- NULL
  markcor.E.size.c.rec.gof <- NULL
  
  if (Year[j] < 2012) {
    
    for (i in Plots) {
      
      markcor.E.size.c.ad[[i]] <- envelope(data.size.c.ad$ppp[[i]], markcorr, nsim = nsim, envelope = TRUE, savefuns = TRUE, savepatterns = TRUE)
      markcor.E.size.c.rec[[i]] <- envelope(data.size.c.rec$ppp[[i]], markcorr, nsim = nsim, envelope = TRUE, savefuns = TRUE, savepatterns = TRUE) 
      
      markcor.E.size.c.ad.cat <- cbind(markcor.E.size.c.ad[[i]]$r, markcor.E.size.c.ad.cat, ppp.cat(markcor.E.size.c.ad[[i]]))
      markcor.E.size.c.rec.cat <- cbind(markcor.E.size.c.rec[[i]]$r, markcor.E.size.c.rec.cat, ppp.cat(markcor.E.size.c.rec[[i]]))
      
      markcor.E.size.c.rec.test <- gof.int(markcor.E.size.c.rec[[i]], gof.win)
      
      markcor.E.size.c.rec.gof <- rbind(markcor.E.size.c.rec.gof, 
                          data.frame(Plot = i, markcor.E.size.c.rec.test))
      
    }
    
    markcor.E.size.c.ad.ctrl <- pool(markcor.E.size.c.ad[[3]], markcor.E.size.c.ad[[4]], markcor.E.size.c.ad[[9]])
    markcor.E.size.c.ad.thnn <- pool(markcor.E.size.c.ad[[2]], markcor.E.size.c.ad[[5]], markcor.E.size.c.ad[[7]])
    markcor.E.size.c.rec.ctrl <- pool(markcor.E.size.c.rec[[3]], markcor.E.size.c.rec[[4]], markcor.E.size.c.rec[[9]], savefuns = TRUE)
    markcor.E.size.c.rec.thnn <- pool(markcor.E.size.c.rec[[2]], markcor.E.size.c.rec[[5]], markcor.E.size.c.rec[[7]], savefuns = TRUE)
    
    markcor.E.size.c.rec.ctrl.test <- gof.int(markcor.E.size.c.rec.ctrl, gof.win)
    markcor.E.size.c.rec.thnn.test <- gof.int(markcor.E.size.c.rec.thnn, gof.win)
  
    markcor.E.size.c.rec.gof <- rbind(markcor.E.size.c.rec.gof, 
                       data.frame(Plot = "Ctrl", markcor.E.size.c.rec.ctrl.test),
                       data.frame(Plot = "Thnn", markcor.E.size.c.rec.thnn.test))
    
    par(mfrow = c(2,2), mar = c(1, 1, 1.25, 1.25), oma = c(4, 4, 2, 2)) 
    plot(markcor.E.size.c.ad.ctrl, legend = F)
    plot(markcor.E.size.c.ad.thnn, legend = F)
    plot(markcor.E.size.c.rec.ctrl, legend = F)
    plot(markcor.E.size.c.rec.thnn, legend = F)
    
  }
  
  Jdif.E.rec <- list(); Jdif.E.sp.rec <- Jdif.E.rec; Jdif.E.size.c.rec <- Jdif.E.rec; markcon.E.sp <- Jdif.E.rec 
  Jdif.E.rec.cat <- NULL; markcon.E.sp.cat <- NULL; Jdif.E.sp.rec.cat <- NULL; Jdif.E.size.c.rec.cat <- NULL
  Jdif.E.rec.gof <- NULL
  
  for (i in Plots) {
    
    Jdif.E.rec[[i]] <- envelope(data.sp.rec$ppp[[i]], Jdif, r = seq(0,1,0.01), nsim = nsim, savefuns = TRUE, savepatterns = TRUE, simulate = expression(rlabel(data.sp.rec$ppp[[i]])))
    Jdif.E.size.c.rec[[i]] <- envelope(data.size.c.rec$ppp[[i]], markcorr, nsim = nsim, envelope = TRUE, savefuns = TRUE, savepatterns = TRUE) 
    
    Jdif.E.rec.cat <- cbind(Jdif.E.rec[[i]]$r, Jdif.E.rec.cat, ppp.cat(Jdif.E.rec[[i]]))
    Jdif.E.size.c.rec.cat <- cbind(Jdif.E.size.c.rec[[i]]$r, Jdif.E.size.c.rec.cat, ppp.cat(Jdif.E.size.c.rec[[i]]))
    Jdif.E.rec.test <- gof.int(Jdif.E.rec[[i]], gof.win)
     
    Jdif.E.rec.gof <- rbind(Jdif.E.rec.gof, 
                                      data.frame(Plot = i, Jdif.E.rec.test))
    
  }
  
  Jdif.E.rec.ctrl <- pool(Jdif.E.rec[[3]], Jdif.E.rec[[4]], Jdif.E.rec[[9]], savefuns = TRUE)
  Jdif.E.rec.thnn <- pool(Jdif.E.rec[[2]], Jdif.E.rec[[5]], Jdif.E.rec[[7]], savefuns = TRUE)
  
  Jdif.E.rec.ctrl.test <- gof.int(markcor.E.size.c.rec.ctrl, gof.win)
  Jdif.E.rec.thnn.test <- gof.int(markcor.E.size.c.rec.thnn, gof.win)
  
  Jdif.E.rec.gof <- rbind(Jdif.E.rec.gof, 
                                    data.frame(Plot = "Ctrl", Jdif.E.rec.ctrl.test),
                                    data.frame(Plot = "Thnn", Jdif.E.rec.thnn.test))
  
  if (j != length(Year)) { Jdif.E.size.c.rec.ctrl <- pool(Jdif.E.size.c.rec[[3]], Jdif.E.size.c.rec[[4]], Jdif.E.size.c.rec[[9]]); Jdif.E.size.c.rec.thnn <- pool(Jdif.E.size.c.rec[[2]], Jdif.E.size.c.rec[[5]], Jdif.E.size.c.rec[[7]]) }
  
  
  par(mfrow = c(2,1), mar = c(1, 1, 1.25, 1.25), oma = c(4, 4, 2, 2)) 
  plot(Jdif.E.rec.ctrl, legend = F)
  plot(Jdif.E.rec.thnn, legend = F)
  
  
  
  # Between fates ...
  
  cat("Test differences in fate \n")
  
  Jdif.E.fate.rec <- list(); K012.E.fate.rec.i <- Jdif.E.fate.rec
  markcor.E.size.c.fate.rec <- list(); null_model_Kmulti <- list()
  Jdif.E.fate.rec.cat <- NULL; K012.E.fate.rec.i.cat <- NULL
  Jdif.E.fate.rec.ctrl <- NULL; Jdif.E.fate.rec.thnn <- NULL
  K012.E.fate.rec.i.ctrl <- NULL; K012.E.fate.rec.i.thnn <- NULL
  Jdif.E.fate.rec.gof <- NULL; K012.E.fate.rec.i.gof <- NULL
  markcor.E.size.c.alive.rec.cat <- NULL; markcor.E.size.c.dead.rec.cat <- NULL
  
  if (Year[j] < 2011) {
    
    for (i in Plots) {
      
      Jdif.E.fate.rec[[i]] <- envelope(data.fate.rec$ppp[[i]], Jdif, i = "0", r = seq(0,1,0.01), nsim = nsim, savefuns = TRUE, savepatterns = TRUE, simulate = expression(rlabel(data.fate.rec$ppp[[i]]))) 
      Jdif.E.fate.rec.cat <- cbind(Jdif.E.fate.rec[[i]]$r, Jdif.E.fate.rec.cat, ppp.cat(Jdif.E.fate.rec[[i]]))
      
      Jdif.E.fate.rec.test <- gof.int(Jdif.E.fate.rec[[i]], gof.win)
      
      Jdif.E.fate.rec.gof <- rbind(Jdif.E.fate.rec.gof, 
                              data.frame(Plot = i, Jdif.E.fate.rec.test) )
      
      K012.E.fate.rec.i[[i]] <- K012(data.fate.rec.i$ppp[[i]], fijo = "Tree", i = "RcSurv", j = "RcDead", r = seq(0, 8, le = length(markcor.E.size.c.ad.ctrl$r)), nsim = nsim, nrank = 5, correction = "Ripley")
      #Kmulti.ls(data.fate.rec.i$ppp[[i]], I = "Tree", J = "RcSurv", r = seq(0, 8, le = 51), corre = "Ripley")
      K012.E.fate.rec.i.cat <- cbind(K012.E.fate.rec.i[[i]]$k01$r, K012.E.fate.rec.i.cat, ppp.cat(K012.E.fate.rec.i[[i]]$k01))
      
      #repito el codigo de la funcion "gof.int" porque la forma de calcularlo es algo diferente
      cons.values <- rollsum(abs(ppp.cat(K012.E.fate.rec.i[[i]]$k01)), gof.win, fill = NA, align = "center") == gof.win
      K012.E.fate.rec.i.interval <- c(min(K012.E.fate.rec.i[[i]]$k01$r[cons.values], na.rm = T), max(K012.E.fate.rec.i[[i]]$k01$r[cons.values], na.rm = T))
      
      K012.E.fate.rec.i.interval0 <- K012.E.fate.rec.i.interval
      
      if (is.infinite(K012.E.fate.rec.i.interval[1]) | (K012.E.fate.rec.i.interval[1]) == (K012.E.fate.rec.i.interval[2])) {
        
        K012.E.fate.rec.i.interval <- c(min(K012.E.fate.rec.i[[i]]$k01$r), max(K012.E.fate.rec.i[[i]]$k01$r))
        
      }

      # lo he hecho de manera mas o menos apañada, metiendo la función base que usa K012; no sé si estará bien hecho, pero es lo mejor que puedo hacer
      marx <- marks(data.fate.rec.i$ppp[[i]])
      fijo <- (marx == "Tree")
      I = (marx == "RcSurv")
      J <- (marx == "RcDead")

      kk <- Kmulti.ls(data.fate.rec.i$ppp[[i]], fijo, I, corre = "isotropic")
      plot(kk)
      
      #k <- envelope(data.fate.rec.i$ppp[[i]], fun = K012, nsim = nsim, i = "RcSurv", j = "RcDead", fijo = "Tree", nrank = 5)
      
      null_model_Kmulti[[i]] <- envelope(data.fate.rec.i$ppp[[i]], fun = Kmulti, nsim = nsim, I = fijo, J = I, nrank = 5, corre = "isotropic", savefuns = TRUE, savepatterns = TRUE)
      
      cons.values <- rollsum(abs(ppp.cat(null_model_Kmulti[[i]])), gof.win, fill = NA, align = "center") == gof.win
      null_model_Kmulti.interval <- c(min(null_model_Kmulti[[i]]$r[cons.values], na.rm = T), max(null_model_Kmulti[[i]]$r[cons.values], na.rm = T))
      
      null_model_Kmulti.interval0 <- null_model_Kmulti.interval
      
      if (is.infinite(null_model_Kmulti.interval[1]) | (null_model_Kmulti.interval[1]) == (null_model_Kmulti.interval[2])) {
        
        null_model_Kmulti.interval <- c(min(null_model_Kmulti[[i]]$r), max(null_model_Kmulti[[i]]$r))
        
      }

      null_model_Kmulti.test <- dclf.test(null_model_Kmulti[[i]], rinterval = null_model_Kmulti.interval, fun = Kmulti, nsim = nsim, I = fijo, J = I, nrank = 5, corre = "isotropic")
      
      if (null_model_Kmulti.interval0[1] == null_model_Kmulti.interval0[2]) {
        
        null_model_Kmulti.test$statistic[1] <- 0
        null_model_Kmulti.test$p.value <- 1
        attributes(null_model_Kmulti.test)$rinterval[1] <- null_model_Kmulti.interval0[1]
        attributes(null_model_Kmulti.test)$rinterval[2] <- null_model_Kmulti.interval0[2]
        
      } 
      
      K012.E.fate.rec.i.gof <- rbind(K012.E.fate.rec.i.gof, 
                          data.frame(Plot = i, r.min = attributes(null_model_Kmulti.test)$rinterval[1], r.max = attributes(null_model_Kmulti.test)$rinterval[2], null_model_Kmulti.test$statistic[1], p.value = null_model_Kmulti.test$p.value))
      
      mark.corr.c.fate <- data.size.c.rec$ppp[[i]] 
      mark.corr.c.fate$marks <- data.frame(SIZE = data.size.c.rec$ppp[[i]]$marks, FATE = data.fate.rec$ppp[[i]]$marks) 
      markcor.E.size.c.fate.rec[[i]] <- markcrosscorr(mark.corr.c.fate, nsim = nsim, correction = "Ripley", savefuns = TRUE, savepatterns = TRUE)
      
      markcor.E.size.c.alive.rec.cat <- cbind(markcor.E.size.c.alive.rec.cat, markcor.E.size.c.fate.rec[[i]]$fns[[3]]$iso)
      markcor.E.size.c.dead.rec.cat <- cbind(markcor.E.size.c.dead.rec.cat, markcor.E.size.c.dead.rec.cat[[i]]$fns[[2]]$iso)
      #remotes::install_github("petrkeil/spasm")
      #PCFr(mark.corr.c.fate, 1, 0, 100, 10)
      #envelope(mark.corr.c.fate, markcrosscorr, nsim = nsim)
      #bivariate mark-correlation function, calcular envelope
      
    }
    
    Jdif.E.fate.rec.ctrl <- pool(Jdif.E.fate.rec[[3]], Jdif.E.fate.rec[[4]], Jdif.E.fate.rec[[9]], savefuns = TRUE)
    Jdif.E.fate.rec.thnn <- pool(Jdif.E.fate.rec[[2]], Jdif.E.fate.rec[[5]], Jdif.E.fate.rec[[7]], savefuns = TRUE)
    
    Jdif.E.fate.rec.ctrl.test <- gof.int(markcor.E.size.c.rec.ctrl, gof.win)
    Jdif.E.fate.rec.thnn.test <- gof.int(markcor.E.size.c.rec.thnn, gof.win)
    
    Jdif.E.fate.rec.gof <- rbind(Jdif.E.fate.rec.gof, 
                            data.frame(Plot = "Ctrl", Jdif.E.fate.rec.ctrl.test),
                            data.frame(Plot = "Thnn", Jdif.E.fate.rec.thnn.test))
    
    K012.E.fate.rec.i.ctrl <- pool(K012.E.fate.rec.i[[3]]$k01, K012.E.fate.rec.i[[4]]$k01, K012.E.fate.rec.i[[9]]$k01)
    K012.E.fate.rec.i.thnn <- pool(K012.E.fate.rec.i[[2]]$k01, K012.E.fate.rec.i[[5]]$k01, K012.E.fate.rec.i[[7]]$k01)
    
    null_model_Kmulti.ctrl <- pool(null_model_Kmulti[[3]], null_model_Kmulti[[4]], null_model_Kmulti[[9]], savefuns = TRUE)
    null_model_Kmulti.thnn <- pool(null_model_Kmulti[[2]], null_model_Kmulti[[5]], null_model_Kmulti[[7]], savefuns = TRUE)
    
    null_model_Kmulti.ctrl.test <- gof.int(null_model_Kmulti.ctrl, gof.win)
    null_model_Kmulti.thnn.test <- gof.int(null_model_Kmulti.thnn, gof.win)

    # K012.E.fate.rec.i.ac.test <- dclf.test(K012.E.fate.rec.i.ctrl, rinterval = K012.E.fate.rec.i.ac.interval, fun = Kmulti.ls, nsim = nsim, I = fijo, J = I, nrank = 5)
    
    K012.E.fate.rec.i.gof <- rbind(K012.E.fate.rec.i.gof,
                          data.frame(Plot = "Ctrl", null_model_Kmulti.ctrl.test),
                          data.frame(Plot = "Thnn", null_model_Kmulti.thnn.test))
    
    markcor.E.size.c.alive.rec.ctrl <- pool(markcor.E.size.c.fate.rec[[3]]$fns[[3]], markcor.E.size.c.fate.rec[[4]]$fns[[3]], markcor.E.size.c.fate.rec[[9]]$fns[[3]])
    markcor.E.size.c.alive.rec.thnn <- pool(markcor.E.size.c.fate.rec[[2]]$fns[[3]], markcor.E.size.c.fate.rec[[5]]$fns[[3]], markcor.E.size.c.fate.rec[[7]]$fns[[3]])
    
    markcor.E.size.c.dead.rec.ctrl <- pool(markcor.E.size.c.fate.rec[[3]]$fns[[2]], markcor.E.size.c.fate.rec[[4]]$fns[[2]], markcor.E.size.c.fate.rec[[9]]$fns[[2]])
    markcor.E.size.c.dead.rec.thnn <- pool(markcor.E.size.c.fate.rec[[2]]$fns[[2]], markcor.E.size.c.fate.rec[[5]]$fns[[2]], markcor.E.size.c.fate.rec[[7]]$fns[[2]])
    
    par(mfrow = c(2,2), mar = c(1, 1, 1.25, 1.25), oma = c(4, 4, 2, 2)) 
    plot(Jdif.E.fate.rec.ctrl, main = NULL, legend = F)
    plot(Jdif.E.fate.rec.thnn, main = NULL, legend = F)
    
    col.trans <- hsv(1, 1, 1, 0)
    plot(K012.E.fate.rec.i.ctrl, sqrt(./pi)-r~r, ylab = expression(L[12]), main = NULL, legend = F, 
         col = c(col.trans, col.trans, col.trans, col.trans, col.trans, col.trans, "grey","black","grey"),
         lty = c(0, 0, 0, 0, 0, 0, 3, 1, 3))
    plot(K012.E.fate.rec.i.thnn, sqrt(./pi)-r~r, ylab = expression(L[12]), main = NULL, legend = F, 
         col = c(col.trans, col.trans, col.trans, col.trans, col.trans, col.trans, "grey","black","grey"),
         lty = c(0, 0, 0, 0, 0 , 0, 3, 1, 3))
    title(paste0("Year ", Year[j], " // Distance-dependent of surviving recruits"), line = -1, cex.main = 1, outer = TRUE)
    par(old.par)
    
    par(mfrow = c(2,2), mar = c(1, 1, 1.25, 1.25), oma = c(4, 4, 2, 2)) 
    plot(markcor.E.size.c.alive.rec.ctrl, main = NULL, legend = F)
    plot(markcor.E.size.c.alive.rec.thnn, main = NULL, legend = F)
    plot(markcor.E.size.c.dead.rec.ctrl, main = NULL, legend = F)
    plot(markcor.E.size.c.dead.rec.thnn, main = NULL, legend = F)
    title(paste0("Year ", Year[j], " // Correlation between size and surviving (above) 
                 and dead (below) recruits"), line = -1, cex.main = 1, outer = TRUE)
    
  }
  
  markcor.E.growth.ad <- list(); Jdif.E.fate.ad <- markcor.E.growth.ad
  markcor.E.growth.ad.cat <- NULL
  markcor.E.growth.ad.ctrl <- NULL; markcor.E.growth.ad.thnn <- NULL
  Jdif.E.fate.ad.cat <- NULL
  
  if (Year[j] == 2009) {
    
    for (i in Plots) {
      
      markcor.E.growth.ad[[i]] <- envelope(data.growth.ad$ppp[[i]], markcorr, nsim = nsim, envelope = TRUE, savefuns = TRUE, savepatterns = TRUE) 
      
      markcor.E.growth.ad.cat <- cbind(markcor.E.growth.ad[[i]]$r, markcor.E.growth.ad.cat, ppp.cat(markcor.E.growth.ad[[i]]))
      
    }
    
    markcor.E.growth.ad.ctrl <- pool(markcor.E.growth.ad[[3]], markcor.E.growth.ad[[4]], markcor.E.growth.ad[[9]]); markcor.E.growth.ad.thnn <- pool(markcor.E.growth.ad[[2]], markcor.E.growth.ad[[5]], markcor.E.growth.ad[[7]])
    
    # Para calcular la mortalidad de los árboles (si es que hay)
    #
    # for (i in c(2,5,7)) {
    #   
    #   Jdif.E.fate.ad[[i]] <- envelope(data.fate.ad$ppp[[i]], Jdif, r = seq(0,1,0.01), nsim = nsim, savefuns=TRUE, savepatterns = TRUE, simulate=expression(rlabel(data.fate.ad$ppp[[i]])))  
    #   
    #   Jdif.E.fate.ad.cat <- cbind(Jdif.E.fate.ad[[i]]$r, Jdif.E.fate.ad.cat, ppp.cat(Jdif.E.fate.ad[[i]]))
    #
    # }
    # 
    # Jdif.E.fate.ad.thnn <- pool(Jdif.E.fate.ad[[2]], Jdif.E.fate.ad[[5]], Jdif.E.fate.ad[[7]])
    # 
    # par(mfrow=c(2,2), mar=c(1, 1, 1.25, 1.25), oma = c(4, 4, 2, 2)) 
    # plot(Jdif.E.fate.ad.thnn, main= NULL, legend = F)
    # 
    # plot(Jdif.E.fate.ad.thnn, main= NULL, legend = F)
    # 
    # title(paste0("Year ", Year[j], " // Growth and fate of trees"), line = -1, cex.main = 1, outer = TRUE)
    # par(old.par)
    
  }
  
  
  
  # Test covariates ----------------------------------------------------
  
  cat("Test environmental covariates \n")
  
  fit.env <- vector("list", length = 9)
  title.env <- character(9)
  L.E.env <- vector("list", length = 9); g.E.env <- L.E.env; kNN.E.env <- L.E.env
  L.E.env.cat <- NULL; g.E.env.cat <- NULL; kNN.E.env.cat <- NULL
  pred.env.sp <- data.frame()
  pred.env.fate <- data.frame()
  g.E.env.gof <- NULL
  
  data.env <- hyperframe(Canopy = listof(rep(NA, 9)), CanOpen = listof(rep(NA, 9)), LAI = listof(rep(NA, 9)), DiffBelow = listof(rep(NA, 9)), DiffBelow.Yr = listof(rep(NA, 9)), N.Sunflecks = listof(rep(NA, 9)), Mdn.Sunflecks = listof(rep(NA, 9)), Max.Sunflecks = listof(rep(NA, 9)), Fs = listof(rep(NA, 9)), Hed = listof(rep(NA, 9)), Rub = listof(rep(NA, 9)), Pter = listof(rep(NA, 9)), Scl = listof(rep(NA, 9)))
  # MDT = listof(rep(NA, 9)), MDS = listof(rep(NA, 9)), Lad.var = listof(rep(NA, 9)), Lad.cv = listof(rep(NA, 9)), Lad.shan = listof(rep(NA, 9)), 
  
  data.sp.rec <- cbind.hyperframe(data.sp.rec, data.env)
  
  matrix.env.sum <- array(unlist(matrix.env), dim = c(15,15))
  colnames(matrix.env.sum) <- gsub('.pred', '', colnames(matrix.env[[2]]))
  rownames(matrix.env.sum) <- gsub('.pred', '', rownames(matrix.env[[2]]))
  sort(colSums(matrix.env.sum > 0.7) - 1)
  sort(rowSums(matrix.env.sum > 0.7) - 1)
  
  par(mfrow = c(1,1))
  corrplot(matrix.env.sum, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)
  
  berman.env <- data.frame()
  var.env.Plot <- list()
  
  for (i in Plots) {
    
    # if (!is.null(data.sp$covs[[i]]$Fs)) {
    
    if (fit.gam == F & i != 9) fit.env.full <- ppm(data.sp.rec$ppp[[i]] ~ 1 + Canopy + CanOpen + LAI + DiffBelow + N.Sunflecks + Max.Sunflecks + Fs + Hed + Pter + Rub + Scl, data = data.sp$covs[[i]])
    if (fit.gam == F & i == 9) fit.env.full <- ppm(data.sp.rec$ppp[[i]] ~ 1 + Canopy + CanOpen + LAI + DiffBelow + Max.Sunflecks + Fs + Hed + Pter + Rub + Scl, data = data.sp$covs[[i]]) # la variable N.Sunflecks en la parcela 9 es homogenea asi que estimaria NA y daría error mas abajo
    if (fit.gam == T) fit.env.full <- ppm(data.sp.rec$ppp[[i]] ~ 1 + bs(Canopy,3) + bs(CanOpen,3) + bs(LAI,3) + bs(DiffBelow,3) + bs(N.Sunflecks,3) + bs(Max.Sunflecks,3) + bs(Fs,3) + bs(Hed,3) + bs(Pter,3) + bs(Rub,3) + bs(Scl,3), use.gam = TRUE, data = data.sp$covs[[i]])
    
    fit.env[[i]] <- step(fit.env.full, trace = 0) 
    title.env[[i]] <- as.character(fit.env[[i]]$trend[[2]])[[2]]
    
    #Para grabar, simplificar objeto (una tabla), ya que ocupa muchisimo espacio (gigas)
    var.env.Plot[[i]] <- summary(fit.env[[i]])$coefs.SE.CI
    
    # for (j in 1:length(var.env2)) {
    #   
    #   berman_tmp <- berman.test(fit.env[[i]],  data.sp$covs[[i]][[c(var.env2[j])]], "Z1")
    #   berman.env <- rbind(berman.env, data.frame(Plot = i, Variable = var.env2[j], Z1 = berman_tmp$statistic, p.value = berman_tmp$p.value, hypothesis = berman_tmp$alternative))
    #   
    # }
    
    data.sp.rec$pred.env[[i]] <- predict(fit.env[[i]], type = "intensity")
    
    # add variables to hyperframe
    data.sp.rec$Canopy[[i]] <- data.sp$covs[[i]]$Canopy
    data.sp.rec$CanOpen[[i]] <- data.sp$covs[[i]]$CanOpen
    data.sp.rec$LAI[[i]] <- data.sp$covs[[i]]$LAI
    data.sp.rec$DiffBelow[[i]] <- data.sp$covs[[i]]$DiffBelow
    data.sp.rec$DiffBelow.Yr[[i]] <- data.sp$covs[[i]]$DiffBelow.Yr
    data.sp.rec$N.Sunflecks[[i]] <- data.sp$covs[[i]]$N.Sunflecks
    data.sp.rec$Mdn.Sunflecks[[i]] <- data.sp$covs[[i]]$Mdn.Sunflecks
    data.sp.rec$Max.Sunflecks[[i]] <- data.sp$covs[[i]]$Max.Sunflecks
    data.sp.rec$Fs[[i]] <- data.sp$covs[[i]]$Fs
    data.sp.rec$Hed[[i]] <- data.sp$covs[[i]]$Hed
    data.sp.rec$Pter[[i]] <- data.sp$covs[[i]]$Pter
    data.sp.rec$Rub[[i]] <- data.sp$covs[[i]]$Rub
    data.sp.rec$Scl[[i]] <- data.sp$covs[[i]]$Scl
    
    g.E.env[[i]] <- envelope(fit.env[[i]], pcf, r = seq(0,4,0.02), nsim = nsim, envelope = TRUE, correction = "Ripley", title = NULL, savefuns = TRUE, savepatterns = TRUE)
    kNN.E.env[[i]] <- envelope(fit.env[[i]], Gest, r = seq(0,1,0.01), nsim = nsim, envelope = TRUE, correction = "rs", title = NULL, savefuns = TRUE, savepatterns = TRUE)
    
    g.E.env.cat <- cbind(g.E.env[[i]]$r, g.E.env.cat, ppp.cat(g.E.env[[i]]))
    kNN.E.env.cat <- cbind(kNN.E.env[[i]]$r, kNN.E.env.cat, ppp.cat(kNN.E.env[[i]]))
    
    g.E.env.test <- gof.int(g.E.env[[i]], gof.win)

    g.E.env.gof <- rbind(g.E.env.gof, 
                         data.frame(Plot = i, g.E.env.test))
    
    pred.env.sp <- rbind(pred.env.sp, data.frame(Plot = i, Treat = Treat[i], Sp = data.sp.rec$ppp[[i]]$marks[data.sp.list[[i]]], 
                                                 pred = fitted(fit.env[[i]], dataonly = T)[data.sp.list[[i]]] )) 
    if (Year[j] < 2011) {
      
      pred.env.fate <- rbind(pred.env.fate, data.frame(Plot = i, Treat = Treat[i], fate = data.fate.rec$ppp[[i]]$marks[data.sp.list[[i]]], 
                                                       pred = fitted(fit.env[[i]], dataonly = T)[data.sp.list[[i]]] ))
      
    }
    
  }
  
  #Plot predicted for each plot
  raster.tmp <- stack(
    raster(data.sp.rec$pred.env$`3`$Qh), raster(data.sp.rec$pred.env$`4`$Qh), raster(data.sp.rec$pred.env$`9`$Qh),
    raster(data.sp.rec$pred.env$`2`[[2]]$Qh), raster(data.sp.rec$pred.env$`5`$Qh), raster(data.sp.rec$pred.env$`7`$Qh))
    
  names(raster.tmp) <- c(paste0("Plot_", c(3,4,9,2,5,7)))
  print( 
    gplot(raster.tmp) + 
      geom_tile(aes(fill = value)) + guides(fill = guide_legend(title = "Prob.")) +
      facet_wrap(~ variable) +
      scale_fill_gradientn(colours = rev(terrain.colors(225))) +
      coord_equal() + theme_classic() 
  )
  
  
  #plot(data.sp.rec$pred.env[Plots], main = "Predicted Environment") # si se quita del modelo especie
  data.sp.rec.Fs <- cbind.hyperframe(data.sp.rec.Fs, data.sp.rec[,7:24])
  data.sp.rec.Qh <- cbind.hyperframe(data.sp.rec.Qh, data.sp.rec[,7:24])
  data.sp.rec.Qi <- cbind.hyperframe(data.sp.rec.Qi, data.sp.rec[,7:24])
  
  # https://www.r-graph-gallery.com/4-barplot-with-error-bar.html
  # Calculates mean, sd, se and IC
  
  sum.env.sp <- pred.env.sp[!is.na(pred.env.sp$Sp),] %>%
    group_by(Sp, Treat) %>%
    summarise(
      n = n(),
      mean = mean(pred, na.rm = T),
      sd = sd(pred, na.rm = T)
    ) %>%
    mutate( se = sd/sqrt(n))  %>%
    mutate( ic = se * qt((1 - 0.05)/2 + .5, n - 1))
  
  if (Year[j] < 2011) {
    
    sum.env.fate <- pred.env.fate[!is.na(pred.env.fate$fate),] %>%
      group_by(fate, Treat) %>%
      summarise(
        n = n(),
        mean = mean(pred, na.rm = T),
        sd = sd(pred, na.rm = T)
      ) %>%
      mutate( se = sd/sqrt(n))  %>%
      mutate( ic = se * qt((1 - 0.05)/2 + .5, n - 1))
    
  }
  
  g.E.env.ctrl <- pool(g.E.env[[3]], g.E.env[[4]], g.E.env[[9]], savefuns = TRUE)
  g.E.env.thnn <- pool(g.E.env[[2]], g.E.env[[5]], g.E.env[[7]], savefuns = TRUE)
  kNN.E.env.ctrl <- pool(kNN.E.env[[3]], kNN.E.env[[4]], kNN.E.env[[9]])
  kNN.E.env.thnn <- pool(kNN.E.env[[2]], kNN.E.env[[5]], kNN.E.env[[7]])
  
  g.E.env.ctrl.test <- gof.int(g.E.env.ctrl, gof.win)
  g.E.env.thnn.test <- gof.int(g.E.env.thnn, gof.win)
  
  g.E.env.gof <- rbind(g.E.env.gof, 
                       data.frame(Plot = "Ctrl", g.E.env.ctrl.test),
                       data.frame(Plot = "Thnn", g.E.env.thnn.test))
  
  par(mfrow = c(2,2), mar = c(1, 1, 1.25, 1.25), oma = c(4, 4, 2, 2)) 
  plot(g.E.env.ctrl, legend = F, mar.panel = c(1, 1, 1, 1), title = NULL)
  plot(kNN.E.env.ctrl, legend = F, mar.panel = c(1, 1, 1, 1), title = NULL)
  
  plot(g.E.env.thnn, legend = F, mar.panel = c(1, 1, 1, 1), title = NULL)
  plot(kNN.E.env.thnn, legend = F, mar.panel = c(1, 1, 1, 1), title = NULL)
  par(old.par)
  
  if (fit.gam == F) fit.env.full <- mppm(ppp ~ 1 + Treat + Canopy + CanOpen + LAI + DiffBelow + N.Sunflecks + Max.Sunflecks + Fs + Hed + Pter + Rub + Scl, data.sp.rec.Fs[Plots,], iformula = ~Interaction:id)
  if (fit.gam == T) fit.env.full <- mppm(ppp ~ 1 + Treat + bs(Canopy,3) + bs(CanOpen,3) + bs(LAI,3) + bs(DiffBelow,3) + bs(N.Sunflecks,3) + bs(Max.Sunflecks,3) + bs(Fs,3) + bs(Hed,3) + bs(Pter,3) + bs(Rub,3) + bs(Scl,3), use.gam = TRUE, data.sp.rec[Plots,], iformula = ~Interaction*id)
  
  fit.env.red <- stepAIC(fit.env.full)
  dev.env <- anova(fit.env.red)
  
  grid.newpage()
  grid.table(dev.env)
  
  if (!is.null(summary(fit.env.red)$coefs.SE.CI)) d <- summary(fit.env.red)$coefs.SE.CI
  if (!is.null(summary(fit.env.red)$Fit)) d <- summary(fit.env.red)$Fit$FIT$coefficients
  d <- data.frame(Variable = rownames(d), d)
  
  grid.newpage()
  estimates.env.red = as_tibble(d) %>% mutate_if(is.numeric, ~sprintf("%.4f",.))
  grid.table(estimates.env.red)
  
  res <- residuals(fit.env.red, type = "pearson")
  smor <- with(hyperframe(res = res), Smooth(res, sigma = 4))
  plot(smor, equal.ribbon = TRUE, main = "Residuals environment")
  
  
  
  # Density covariates
  
  cat("Test density covariates \n")
  
  fit.dens <- vector("list", length = 9)
  title.dens <- character(9)
  L.E.dens <- vector("list", length = 9); g.E.dens <- L.E.dens; kNN.E.dens <- L.E.dens
  L.E.dens.cat <- NULL; g.E.dens.cat <- NULL; kNN.E.dens.cat <- NULL
  pred.dens.sp <- data.frame()
  pred.dens.fate <- data.frame()
  g.E.dens.gof <- NULL
  
  data.dens <- hyperframe(Dens = listof(rep(NA, 9)), Dens.adult = listof(rep(NA, 9)), Dens.rec = listof(rep(NA, 9)), Dens.size = listof(rep(NA, 9)), Dens.size.rec = listof(rep(NA, 9)), Dens.Fs.rec = listof(rep(NA, 9)), Dens.Qh.rec = listof(rep(NA, 9)), Dens.Qi.rec = listof(rep(NA, 9)), Dens.rich.rec = listof(rep(NA, 9)), Dens.shan.rec = listof(rep(NA, 9)))
  
  data.sp.rec <- cbind.hyperframe(data.sp.rec, data.dens)
  
  matrix.dens.sum <- array(unlist(matrix.dens), dim = c(10,10))
  colnames(matrix.dens.sum) <- colnames(matrix.dens[[2]])
  rownames(matrix.dens.sum) <- rownames(matrix.dens[[2]])
  sort(colSums(matrix.dens.sum > 0.7) - 1)
  sort(rowSums(matrix.dens.sum > 0.7) - 1)
  
  par(mfrow = c(1,1))
  corrplot(matrix.dens.sum, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)
  
  berman.dens <- data.frame()
  
  var.dens.Plot <- list()
  
  for (i in Plots) {
    
    if (fit.gam == F) fit.dens.full <- ppm(data.sp.rec$ppp[[i]] ~ 1 + Dens.adult + Dens.size.rec + Dens.rich.rec + Dens.shan.rec, data = data.sp$covs[[i]])
    if (fit.gam == T) fit.dens.full <- ppm(data.sp.rec$ppp[[i]] ~ 1 + bs(Dens.adult,3) + bs(Dens.size.rec,3) + bs(Dens.rich.rec,3) + bs(Dens.shan.rec,3), use.gam = TRUE, data = data.sp$covs[[i]])
    
    fit.dens[[i]] <- step(fit.dens.full, trace = 0)
    title.dens[[i]] <- as.character(fit.dens[[i]]$trend)[[2]]
    
    var.dens.Plot[[i]] <- summary(fit.dens[[i]])$coefs.SE.CI
    
    # for (j in 1:length(var.dens)) {
    #   
    #   berman_tmp <- berman.test(fit.dens[[i]],  data.sp$covs[[i]][[c(var.dens2[j])]], "Z2")
    #   berman.dens <- rbind(berman.dens, data.frame(Plot = i, Variable = var.dens2[j], Z2 = berman_tmp$statistic, p.value = berman_tmp$p.value, hypothesis = berman_tmp$alternative))
    #   
    # }
    
    data.sp.rec$pred.dens[[i]] <- predict(fit.dens[[i]], type = "intensity")
    
    # add variables to hyperframe
    data.sp.rec$Dens[[i]] <- data.sp$covs[[i]]$Dens
    data.sp.rec$Dens.adult[[i]] <- data.sp$covs[[i]]$Dens.adult
    data.sp.rec$Dens.rec[[i]] <- data.sp$covs[[i]]$Dens.rec
    data.sp.rec$Dens.size[[i]] <- data.sp$covs[[i]]$Dens.size
    data.sp.rec$Dens.size.rec[[i]] <- data.sp$covs[[i]]$Dens.size.rec
    data.sp.rec$Dens.Fs.rec[[i]] <- data.sp$covs[[i]]$Dens.Fs.rec
    data.sp.rec$Dens.Qh.rec[[i]] <- data.sp$covs[[i]]$Dens.Qh.rec
    data.sp.rec$Dens.Qi.rec[[i]] <- data.sp$covs[[i]]$Dens.Qi.rec
    data.sp.rec$Dens.rich.rec[[i]] <- data.sp$covs[[i]]$Dens.rich.rec
    data.sp.rec$Dens.shan.rec[[i]] <- data.sp$covs[[i]]$Dens.shan.rec
    
    g.E.dens[[i]] <- envelope(fit.dens[[i]], pcf, r = seq(0,4,0.02), nsim = nsim, envelope = TRUE, correction = "Ripley", title = NULL, savefuns = TRUE, savepatterns = TRUE)
    kNN.E.dens[[i]] <- envelope(fit.dens[[i]], Gest, r = seq(0,1,0.01), nsim = nsim, envelope = TRUE, correction = "rs", title = NULL, savefuns = TRUE, savepatterns = TRUE)
    
    g.E.dens.cat <- cbind(g.E.dens[[i]]$r, g.E.dens.cat, ppp.cat(g.E.dens[[i]]))
    kNN.E.dens.cat <- cbind(kNN.E.dens[[i]]$r, kNN.E.dens.cat, ppp.cat(kNN.E.dens[[i]]))
    
    g.E.dens.test <- gof.int(g.E.dens[[i]], gof.win)

    g.E.dens.gof <- rbind(g.E.dens.gof, 
                         data.frame(Plot = i, g.E.dens.test))
    
    pred.dens.sp <- rbind(pred.dens.sp, data.frame(Plot = i, Treat = Treat[i], Sp = data.sp.rec$ppp[[i]]$marks[data.sp.list[[i]]], 
                                                   pred = fitted(fit.dens[[i]], dataonly = T)[data.sp.list[[i]]] )) 
    if (j != length(Year)) {
      
      pred.dens.fate <- rbind(pred.dens.fate, data.frame(Plot = i, Treat = Treat[i], fate = data.fate.rec$ppp[[i]]$marks[data.sp.list[[i]]], 
                                                         pred = fitted(fit.dens[[i]], dataonly = T)[data.sp.list[[i]]] ))
      
    }
    
  }
  
  #Plot predicted for each plot
  raster.tmp <- stack(
    raster(data.sp.rec$pred.dens$`3`$Qh), raster(data.sp.rec$pred.dens$`4`$Qh), raster(data.sp.rec$pred.dens$`9`$Qh),
    raster(data.sp.rec$pred.dens$`2`[[2]]$Qh), raster(data.sp.rec$pred.dens$`5`$Qh), raster(data.sp.rec$pred.dens$`7`$Qh))
  
  names(raster.tmp) <- c(paste0("Plot_", c(3,4,9,2,5,7)))
  print( 
    gplot(raster.tmp) + 
      geom_tile(aes(fill = value)) + guides(fill = guide_legend(title = "Prob.")) +
      facet_wrap(~ variable) +
      scale_fill_gradientn(colours = rev(terrain.colors(225))) +
      coord_equal() + theme_classic() 
  )
  

  #plot(data.sp.rec$pred.dens[Plots], main = "Predicted Aggregation")
  data.sp.rec.Fs <- cbind.hyperframe(data.sp.rec.Fs, data.sp.rec[,25:35])
  data.sp.rec.Qh <- cbind.hyperframe(data.sp.rec.Qh, data.sp.rec[,25:35])
  data.sp.rec.Qi <- cbind.hyperframe(data.sp.rec.Qi, data.sp.rec[,25:35])
  
  
  # https://www.r-graph-gallery.com/4-barplot-with-error-bar.html
  # Calculates mean, sd, se and IC
  sum.dens.sp <- pred.dens.sp[!is.na(pred.dens.sp$Sp),] %>%
    group_by(Sp, Treat) %>%
    summarise(
      n = n(),
      mean = mean(pred, na.rm = T),
      sd = sd(pred, na.rm = T)
    ) %>%
    mutate( se = sd/sqrt(n))  %>%
    mutate( ic = se * qt((1 - 0.05)/2 + .5, n - 1))
  
  if (j != length(Year)) {
    
    sum.dens.fate <- pred.dens.fate[!is.na(pred.dens.fate$fate),] %>%
      group_by(fate, Treat) %>%
      summarise(
        n = n(),
        mean = mean(pred, na.rm = T),
        sd = sd(pred, na.rm = T)
      ) %>%
      mutate( se = sd/sqrt(n))  %>%
      mutate( ic = se * qt((1 - 0.05)/2 + .5, n - 1))
    
  }
  
  g.E.dens.ctrl <- pool(g.E.dens[[3]], g.E.dens[[4]], g.E.dens[[9]], savefuns = TRUE)
  g.E.dens.thnn <- pool(g.E.dens[[2]], g.E.dens[[5]], g.E.dens[[7]], savefuns = TRUE)
  kNN.E.dens.ctrl <- pool(kNN.E.dens[[3]], kNN.E.dens[[4]], kNN.E.dens[[9]])
  kNN.E.dens.thnn <- pool(kNN.E.dens[[2]], kNN.E.dens[[5]], kNN.E.dens[[7]])
  
  g.E.dens.ctrl.test <- gof.int(g.E.dens.ctrl, gof.win)
  g.E.dens.thnn.test <- gof.int(g.E.dens.thnn, gof.win)
  
  g.E.dens.gof <- rbind(g.E.dens.gof, 
                       data.frame(Plot = "Ctrl", g.E.dens.ctrl.test),
                       data.frame(Plot = "Thnn", g.E.dens.thnn.test))
  
  par(mfrow = c(2,2), mar = c(1, 1, 1.25, 1.25), oma = c(4, 4, 2, 2)) 
  plot(g.E.dens.ctrl, legend = F, mar.panel = c(1, 1, 1, 1), title = NULL)
  plot(kNN.E.dens.ctrl, legend = F, mar.panel = c(1, 1, 1, 1), title = NULL)
  plot(g.E.dens.thnn, legend = F, mar.panel = c(1, 1, 1, 1), title = NULL)
  plot(kNN.E.dens.thnn, legend = F, mar.panel = c(1, 1, 1, 1), title = NULL)
  par(old.par)
  
  if (fit.gam == F) fit.dens.full <- mppm(ppp ~ 1 + Treat +  Dens.adult + Dens.size.rec + Dens.rich.rec + Dens.shan.rec, data.sp.rec[Plots,], iformula = ~Interaction*id)
  if (fit.gam == T) fit.dens.full <- mppm(ppp ~ 1 + Treat + bs(Dens.adult,3) + bs(Dens.size.rec,3) + bs(Dens.rich.rec,3) + bs(Dens.shan.rec,3),  use.gam = TRUE, data.sp.rec[Plots,], iformula = ~Interaction*id)
  
  fit.dens.red <- stepAIC(fit.dens.full)
  dev.dens <- anova(fit.dens.red)
  grid.newpage()
  grid.table(dev.dens)
  
  if (!is.null(summary(fit.dens.red)$coefs.SE.CI)) d <- summary(fit.dens.red)$coefs.SE.CI
  if (!is.null(summary(fit.dens.red)$Fit)) d <- summary(fit.dens.red)$Fit$FIT$coefficients
  d <- data.frame(Variable = rownames(d), d)
  
  grid.newpage()
  estimates.dens.red = as_tibble(d) %>% mutate_if(is.numeric, ~sprintf("%.4f",.))
  grid.table(estimates.dens.red)
  
  res <- residuals(fit.dens.red, type = "pearson")
  smor <- with(hyperframe(res = res), Smooth(res, sigma = 4))
  plot(smor, equal.ribbon = TRUE, main = "Residuals density")
  
  if (save.output == T) dev.off()
  
  #object.size 
  if (save.output == T) save(
    
    fit.clust, sum.clust,
    g.E.t.ctrl.test, g.E.t.thnn.test,
    g.E.t.cat, g.E.t.gof,
    
    #L.E.rec.cat, 
    g.E.rec.cat, kNN.E.rec.cat, #F.E.rec.cat,
    #L.E.rec.ctrl, 
    g.E.rec.ctrl, kNN.E.rec.ctrl, #F.E.rec.ctrl, L.E.rec.thnn, 
    g.E.rec.thnn, kNN.E.rec.thnn, 
    g.E.rec.gof,
    #F.E.rec.thnn,
    
    #L.E.sp.cat, g.E.sp.cat, kNN.E.sp.cat, Jdif.E.sp.cat, markcon.E.sp.cat,
    #L.E.sp.ctrl, g.E.sp.ctrl, kNN.E.sp.ctrl, Jdif.E.sp.ctrl, L.E.sp.thnn, g.E.sp.thnn, kNN.E.sp.thnn, Jdif.E.sp.thnn, 
    
    #L.E.size.cat, 
    g.E.size.cat, kNN.E.size.cat, Jdif.E.size.cat, markcon.E.size.cat, g.E.ac.cat,
    #L.E.size.ctrl, 
    g.E.size.ctrl, kNN.E.size.ctrl, g.E.ac.ctrl, #Jdif.E.size.ctrl, 
    #L.E.size.thnn, 
    g.E.size.thnn, kNN.E.size.thnn, g.E.ac.thnn, #Jdif.E.size.thnn,
    g.E.size.ctrl1, g.E.size.ctrl2, g.E.size.thnn1, g.E.size.thnn2,
    g.E.size.gof,
    g.E.ac.gof,
    markcon.E.size.ctrl, #markcon.E.sp.thnn, 
    
    markcor.E.size.c.ad, markcor.E.size.c.rec, 
    markcor.E.size.c.rec.thnn, markcor.E.size.c.rec.ctrl,
    markcor.E.size.c.rec.gof, markcor.E.size.c.rec.cat,
    
    #####
    Jdif.E.rec, Jdif.E.sp.rec, Jdif.E.size.c.rec, markcon.E.sp,
    Jdif.E.rec.cat, Jdif.E.rec.gof, 
    markcon.E.sp.cat, Jdif.E.sp.rec.cat, Jdif.E.size.c.rec.cat,
    Jdif.E.rec.ctrl, Jdif.E.rec.thnn, #Jdif.E.sp.rec.ctrl, Jdif.E.sp.rec.thnn,
    
    #markcon.E.sp.ctrl, markcon.E.sp.thnn,
    Jdif.E.size.c.rec.ctrl, Jdif.E.size.c.rec.thnn,
    Jdif.E.fate.rec, K012.E.fate.rec.i,
    Jdif.E.fate.rec.cat, Jdif.E.fate.rec.gof,
    K012.E.fate.rec.i.cat, K012.E.fate.rec.i.gof,
    Jdif.E.fate.rec.ctrl, Jdif.E.fate.rec.thnn, K012.E.fate.rec.i.ctrl, K012.E.fate.rec.i.thnn,
    null_model_Kmulti.ctrl, null_model_Kmulti.thnn, 
    Jdif.E.fate.ad.cat, 
    markcor.E.size.c.alive.rec.cat, markcor.E.size.c.dead.rec.cat,
    markcor.E.size.c.alive.rec.thnn, markcor.E.size.c.alive.rec.ctrl,
    markcor.E.size.c.dead.rec.thnn, markcor.E.size.c.dead.rec.ctrl,
    #####

    markcor.E.growth.ad, Jdif.E.fate.ad,
    markcor.E.growth.ad.cat,
    markcor.E.growth.ad.ctrl, markcor.E.growth.ad.thnn,
    
    data.sp.rec,
    
    matrix.env.sum,
    #fit.env, title.env, L.E.env, g.E.env, kNN.E.env, 
    pred.env.sp, sum.env.sp,
    pred.env.fate, sum.env.fate, #####
    
    #L.E.env.cat, 
    g.E.env.cat, kNN.E.env.cat,
    #L.E.env.ctrl, L.E.env.thnn, 
    g.E.env.ctrl, g.E.env.thnn, kNN.E.env.ctrl, kNN.E.env.thnn, 
    g.E.env.gof,
    estimates.env.red, dev.env, #fit.env.red, 
    berman.env, var.env.Plot,
    
    matrix.dens.sum,
    #fit.dens, title.dens, L.E.dens, g.E.dens, kNN.E.dens, 
    pred.dens.sp, sum.dens.sp, 
    pred.dens.fate, sum.dens.fate, #####
    
    #L.E.dens.cat, 
    g.E.dens.cat, kNN.E.dens.cat,
    #L.E.dens.ctrl, L.E.dens.thnn, 
    g.E.dens.ctrl, g.E.dens.thnn, kNN.E.dens.ctrl, kNN.E.dens.thnn, 
    g.E.dens.gof,
    estimates.dens.red, dev.dens, #fit.dens.red,
    berman.dens, var.dens.Plot,

    compress = TRUE,
    file = paste0("~/Documentos/Datos NO publicados/BioIntForest/PPA_Results/Reports/PPA_", Year[j],".RData"))
  
}