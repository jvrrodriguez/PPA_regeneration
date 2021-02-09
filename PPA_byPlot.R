##########################################
####   ANALIZE POINT DATA           ######
##########################################

library(shapefiles)
library(spatstat)
library(spatstat.local)
library(spatstat.utils)
library(raster)
#library(spatialEco) #https://github.com/jeffreyevans/spatialEco
library(shar) #https://r-spatialecology.github.io/shar/
library(ecespa)

#datos con la informacion de posicion de las plantas
setwd("~/Documentos/Datos NO publicados/BioIntForest/Data/Inventarios_floristicos")
data.frond <- data.frame(read.table(file="cruce_sombras_plantulas_2008_2011_JVR.csv", header=T, sep=";",dec=","))
data.pinus <- data.frame(read.table(file="PS09Medidas pino 2009 con coordenadas.csv", header=T, sep=";",dec=","))
data.pinus99 <- data.frame(read.table(file="Medidas pinos 1999.csv", header=T, sep=";",dec=","))

#Treatment plots
Treat <-  c("40%", "30%", "0%", "0%", "30%", "40%", "30%", "40%", "0%") #ONLY 0% and 30%

data.frond$SP2 <- ifelse(data.frond$SP=="Qi" | data.frond$SP=="Qir" | data.frond$SP=="Qic" | data.frond$SP=="Qi?", "Qi",
                         ifelse(data.frond$SP=="Qh" | data.frond$SP=="Qhr", "Qh",
                                ifelse(data.frond$SP=="Fs" | data.frond$SP=="Fsr", "Fs",NA)))
data.frond$SP2<-as.factor(data.frond$SP2)
data.frond$RES <- ifelse(data.frond$SP=="Qir" | data.frond$SP=="Qhr" | data.frond$SP=="Fsr", 1, 0)

range01 <- function(x){
  x <- ifelse(x == -9999, NA, x) 
  x0 <- (x-min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T))
  return(x0)}

#convert to an hyperframe plots
sizejuv <- 0.2 # umbral (m) para discriminar juvenil de plantula
sizead <- 2 # umbral (m) para discriminar adulto de juvenil
Year <- c(2008, 2009, 2010, 2011)
nsim <- 199


for (i in unique(data.frond$PARCELA)) {
  
  cat("Simulating Plot...", i)
  
  pdf(paste0("~/Documentos/Datos NO publicados/BioIntForest/Results/Results.PPA_A", i,".pdf"), width=7, height=5)
  
  data.plot <- subset(data.frond, PARCELA==i & SP2!="<NA>") 
  
  #stantardize relative positions only for decidious trees
  data.plot$MAPX <- data.plot$MAPX - min(data.plot$MAPX)
  data.plot$MAPY <- data.plot$MAPY - min(data.plot$MAPY)
  p.range <- c(round(min(data.plot$MAPX)), round(max(data.plot$MAPX)), round(min(data.plot$MAPY)), round(max(data.plot$MAPY)))
  
  data.pinus99.plot <- subset(data.pinus99, Parcela==paste0("A",i) & is.finite(DM))
  data.pinus99.plot$Size <- rowMeans(cbind(data.pinus99.plot$P.MeanHt, data.pinus99.plot$A.MeanHt), na.rm =T)
  x <- data.pinus99.plot$DM[!is.nan(data.pinus99.plot$Size)]
  y <- data.pinus99.plot$Size[!is.nan(data.pinus99.plot$Size)]
  
  # Modelo para calcular las relacion alometrica entre dbh y altura
  m <- nls(y~a*x/(b+x)) 
  cor(y, predict(m))
  
  data.pinus.plot <- subset(data.pinus, PARCELA==paste0("A",i))
  data.pinus.plot <- data.pinus.plot[!is.na(data.pinus.plot$Dn_NS_07) | !is.na(data.pinus.plot$Dn_EO_07),]
  
  #Load environmental maps
  
  range01.raster <- function(x){(x-min(x))/(max(x)-min(x))}
  
  load(paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/Environ_Maps/MDS-MDT_A", i, ".RData"))
  
  load(paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/Environ_Maps/Light_A", i, ".RData"))
  
  # DELETE AFTER NEW CALCULATIONS
  if (i < 7) {
    raster.Light <- raster.krig
    im.Light <- im.krig 
  } else {
    raster.Light <- ras.Light
  }
  # DELETE AFTER NEW CALCULATIONS
  
  load(paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/Environ_Maps/Understory_A", i, ".RData"))
  
  load(paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/Environ_Maps/Canopy_A", i, ".RData"))
  
  load(paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/Environ_Maps/Lad_A", i, ".RData"))
  
  
  
  # Analysis for each year --------------------------------------------------
  
  for (j in 1:length(Year)){
    
    cat("Simulating...", Year[j])
    
    # Setup data --------------------------------------------------------------
    
    
    if (Year[j] < 2009){
      
      data.pinus.year <- data.frame(x = data.pinus.plot$X, y = data.pinus.plot$Y,
                                    SP = "Ps", 
                                    SIZE = predict(m, newdata = data.frame(x = rowMeans(cbind(data.pinus.plot$Dn_NS_07, data.pinus.plot$Dn_EO_07), na.rm=T))),
                                    FATE = ifelse(is.na(data.pinus.plot$Muerto.07) | data.pinus.plot$Muerto.07 == 0, 1, 0))
      
    } else {
      
      data.pinus.year <- data.pinus.plot[is.na(data.pinus.plot$Muerto.07) | is.na(data.pinus.plot$Muerto.09),]
      data.pinus.year <- data.frame(x = data.pinus.year$X, y = data.pinus.year$Y,
                                    SP = "Ps", 
                                    SIZE = predict(m, newdata = data.frame(x = data.pinus.year$Dn_09)),
                                    FATE = ifelse(is.na(data.pinus.year$Muerto.09) | data.pinus.year$Th_09 == 0, 1, 0))
    }
    
    if (Year[j] == 2009){
    
      data.pinus.year <- data.frame(data.pinus.year, 
                                    GROWTH =  round(data.pinus.plot$Dn_09 - rowMeans(cbind(data.pinus.plot$Dn_NS_07, data.pinus.plot$Dn_EO_07), na.rm=T), 1)) 
      data.pinus.year$GROWTH <- log(ifelse(data.pinus.year$GROWTH < 0, 0, data.pinus.year$GROWTH) + 1)
    
    } else {
      
      data.pinus.year <- data.frame(data.pinus.year, GROWTH =  0)
      
    }
    
    
    if (j != length(Year)){
      
      data.plot$Size <- data.plot[paste0("ALT",Year[j])][[1]]
      data.sel <- data.plot[which(grepl("ALT", colnames(data.plot)))][(j+1):length(Year)]
      data.plot$Size <- ifelse(rowSums(data.sel== 201.0 | data.sel == -9999.0 | is.na(data.plot$Size), na.rm=T) > 0, 
                               apply(data.sel, 1, FUN=min, na.rm=T), data.plot$Size)
      data.plot$Size <- ifelse(is.infinite(data.plot$Size), NA, data.plot$Size)
      data.plot$Size <- ifelse(data.plot$Size == -9999.0, 201, data.plot$Size)
      data.plot$Fate <- as.factor(ifelse(is.na(rowSums(data.sel)), 0, 1))
      
      data.year <- data.plot[is.finite(data.plot$Size), ]
      data.year <- data.frame(x = data.year$MAPX, y = data.year$MAPY, SP = as.factor(data.year$SP2), SIZE = data.year$Size / 100, FATE = data.year$Fate, GROWTH = 0)
      data.year <- rbind(data.year, data.pinus.year)
      data.year <- data.year[!is.na(data.year$SIZE),]
      data.year$SP <- as.factor(as.character(data.year$SP))
      
      data.year <- with(data.year, 
                        ppp(x = x, y = y, marks = data.frame(SP = as.factor(SP), SIZE = log(SIZE + 1), FATE = FATE, GROWTH = GROWTH), 
                            owin(xrange=c(p.range[1],p.range[2]), yrange=c(p.range[3],p.range[4]))))
      
      layout(matrix(c(1,1,2,3), 2, 2, byrow = FALSE), heights = c(3, 3), widths = c(3, 3))
      par(mar = c(3, 3, 3, 4))
      plot(data.year$x,data.year$y,col=as.integer(data.year$marks$SP), 
           cex = data.year$marks$SIZE, pch = ifelse(data.year$marks$FATE=="0", 19, 1), xlab = "x coords", ylab = "y coords", main = NULL)
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
                            owin(xrange=c(p.range[1],p.range[2]), yrange=c(p.range[3],p.range[4]))))
      
      layout(matrix(c(1,1,2,3), 2, 2, byrow = FALSE), heights = c(2, 2), widths = c(2, 2))
      par(mar = c(3, 3, 3, 4))
      plot(data.year$x,data.year$y,col=as.integer(data.year$marks$SP), 
           cex = data.year$marks$SIZE, xlab = "x coords", ylab = "y coords", main = NULL)
      title(paste0("Plot A", i, "; Yr: ", Year[j], "; Th: ",Treat[i], "; n: ", data.year$n), line = -1, cex.main = 1, outer = TRUE)
      
    }
    
    Q <- quadratcount(data.year, nx = p.range[2]/2, ny = p.range[4]/2)
    quadrat.test.ppp(data.year, nx = p.range[2]/2, ny = p.range[4]/2)
    plot(Q, add = TRUE, cex = 0.5, col="grey")
    
    par(mar = c(3, 0, 3, 3))
    hist(array(Q), xlab="density plants",  main="")
    par(mar = c(3, 0, 3, 3))
    hist(data.year$marks$SIZE, xlab="log(height+1)", main="")
    abline (v=log(sizejuv + 1), col="grey"); abline(v=log(sizead + 1), col="grey")
    layout(matrix(c(1), 1, 1, byrow = FALSE))
    
    dens.all <- density(data.year, sigma=bw.diggle(data.year))
    
    data.sp <- cut.ppp(data.year, "SP")
    # if (sum(data.sp$marks == "Fs") < 6) { #Este es el grupo menos frecuente
    #   data.sp <- data.sp[as.character(data.sp$marks) != "Fs",]
    # }
    data.sp$marks <- as.factor(as.character(data.sp$marks))
    plot(split(data.sp), main = paste0("Plot A", i, "; Year = ", Year[j]))
    
    dens.sp <- density(split(data.sp), sigma=bw.diggle(data.sp))
    
    data.size <- cut.ppp(data.year, "SIZE", breaks=c(-0.1, log(sizejuv + 1), log(sizead + 1), Inf), labels=c("Recruit", "Recruit", "Adult"))
    
    if (sum(is.na(data.size$marks)) > 0) {
      data.size <- data.size[!is.na(data.size$marks)]}
    
    data.size$marks <- as.factor(as.character(data.size$marks))
    plot(split(data.size), main = paste0("Plot A", i, "; Year = ", Year[j]))
    
    dens.size <- density(split(data.size), sigma=bw.diggle(data.size))
    
    data.size.c <- subset.ppp(data.year, select = "SIZE")
    data.size.c <- data.size.c[!is.na(as.numeric(data.size.c$marks)),]
    dens.size.c <- density(data.size.c, sigma=bw.diggle(data.size.c))
    
    data.sp.reg <- data.sp[which(data.size$marks != "Adult"),]
    data.size.reg <- data.size[which(data.size$marks!= "Adult"),]
    data.size.c.reg <- data.size.c[which(data.size$marks != "Adult"),]
    
    data.size.c.ad <- data.size.c[which(data.size$marks == "Adult"),]
    
    dens.size.c.ad <- density(data.size.c.ad, sigma=bw.diggle(data.size.c.ad))
    dens.size.c.reg <- density(data.size.c.reg, sigma=bw.diggle(data.size.c.reg))
    
    if (sum(table(data.sp.reg$marks) < 6) > 0) { #Este es el grupo menos frecuente
      
      rm.sp.reg <- names(which(table(data.sp.reg$marks) < 6)) 
      
      for (k in 1:length(rm.sp.reg)) {
        data.sp.reg <- data.sp.reg[which(data.sp.reg$marks != rm.sp.reg[k]),] 
      }
      
    }
    
    data.sp.reg$marks <- as.factor(as.character(data.sp.reg$marks))
    data.size.reg$marks <- as.factor(as.character(data.size.reg$marks))
    
    if (j != length(Year)) {
      
      data.fate <- cut.ppp(data.year, "FATE")
      data.fate.reg <- data.fate[which(data.size$marks != "Adult"),]
      data.fate.reg.i <- data.fate
      data.fate.reg.i$marks <- ifelse(as.character(data.size$marks) == "Adult", "Tree",
                                      ifelse(as.character(data.size$marks) != "Adult" & data.fate$marks == 0, "RcDead",
                                             ifelse(as.character(data.size$marks) != "Adult" & data.fate$marks == 1, "RcSurv", NA)))
      
      if (sum(is.na(data.fate.reg.i$marks)) > 0) {
        
        data.fate.reg.i <- data.fate.reg.i[-which(is.na(data.fate.reg.i$marks))] } 
        plot(split(data.fate), main = paste0("Plot A", i, "; Year = ", Year[j]))
      
    }
    
    if (Year[j] == 2009){
      
      data.growth <- data.year
      data.growth$marks <- data.growth$marks[4][[1]]
      data.growth.ad <- data.growth[data.size$marks == "Adult" & data.sp$marks == "Ps",]
      
      data.fate.ad <- data.fate[data.size$marks == "Adult" & data.sp$marks == "Ps",]
      
    }
    
    
    plot(dens.sp)
    plot(dens.size)
    plot(dens.size.c)
    plot(dens.size.c.ad)
    plot(dens.size.c.reg)
    
    
    # Test general patterns ----------------------------------------------------
    
    L.E <- envelope(data.size, Lest, nsim=nsim, fix.n=TRUE, correction="Ripley")
    g.E <- envelope(data.size, pcf, nsim=nsim, fix.n=TRUE, correction="Ripley")
    kNN.E <- envelope(data.size, Gest, nsim=nsim, fix.n=TRUE, correction="rs")
    F.E <- envelope(data.size, Fest, nsim=nsim, fix.n=TRUE, correction="rs")
    
    par(mfrow=c(2,2), mar=c(1, 1, 1.25, 1.25), oma = c(4, 4, 2, 2)) 
    plot(L.E, . -r ~ r, shade=c("hi", "lo"), legend = F, main = NULL)
    plot(g.E, main = NULL, legend = F)
    plot(kNN.E, main = NULL, legend = F)
    plot(F.E, main = NULL, legend = F)
    title(paste("Plot", i), line = 0, outer = TRUE)
    old.par <- par(mfrow=c(1,1), mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
    
    
    
    # Test spatial patterns of marks --------------------------------------------
    
    # Compare species
    
    #data.size.sp <- split(data.size)$Fs
    L.E.sp <- alltypes(data.sp, Ldot, nsim = nsim, envelope = TRUE, correction="Ripley", title = NULL) 
    g.E.sp <- alltypes(data.sp, pcfdot, nsim = nsim, envelope = TRUE, correction="Ripley", title = NULL) 
    kNN.E.sp <- alltypes(data.sp, Gdot, nsim = nsim, envelope = TRUE, correction="rs", title = NULL)
    #F.E.sp <- alltypes(data.sp, Fdot, nsim = nsim, envelope = TRUE, correction="rs")
    
    plot(L.E.sp, . -r ~ r, shade=c("hi", "lo"), legend = F, mar.panel=c(1, 1, 1, 1), title = NULL)
    plot(g.E.sp, legend = F, mar.panel=c(1, 1, 1, 1), title = NULL)
    plot(kNN.E.sp, legend = F, mar.panel=c(1, 1, 1, 1), title = NULL)
    #plot(F.E.sp, main = paste("Plot", i), legend = F, mar.panel=c(0.5, 0.5, 0.2, 0.2))
    
    # Compare sizes
    
    L.E.size <- alltypes(data.size, Ldot, nsim = nsim, envelope = TRUE, correction="Ripley") 
    g.E.size <- alltypes(data.size, pcfdot, nsim = nsim, envelope = TRUE, correction="Ripley") 
    kNN.E.size <- alltypes(data.size, Gdot, nsim = nsim, envelope = TRUE, correction="rs")
    #F.E.size <- alltypes(data.size, Fdot, nsim = nsim, envelope = TRUE, correction="rs")
    
    plot(L.E.size, . -r ~ r, shade=c("hi", "lo"), legend = F, mar.panel=c(1, 1, 1, 1), title = NULL)
    plot(g.E.size, legend = F, mar.panel=c(1, 1, 1, 1), title = NULL)
    plot(kNN.E.size, legend = F, mar.panel=c(1, 1, 1, 1), title = NULL)
    #plot(F.E.size, main = paste("Plot", i), legend = F), mar.panel=c(0.5, 0.5, 0.2, 0.2)
    old.par
    
    
    # Test covariates ----------------------------------------------------
    
    # Environmental covariates
    
    if (length(im.Understory) == 5) {
      
      im.covs <- listof(MDT = im.MDTS[[1]], MDS = im.MDTS[[2]],
                        Lad.var = im.lad[[1]], Lad.cv= im.lad[[2]], Lad.shan = im.lad[[3]],
                        CanOpen = im.Light[[1]], LAI = im.Light[[2]], DirectBelow = im.Light[[3]], DiffBelow = im.Light[[4]], DirectBelow.Year = im.Light[[5]], DiffBelow.Year = im.Light[[6]],
                        Canopy = im.Canopy[[1]], 
                        Fs = im.Understory[[1]], Hed = im.Understory[[2]], Pter = im.Understory[[3]], Rub = im.Understory[[4]], Scl = im.Understory[[5]],
                        Dens = dens.all, Dens.Adult = dens.size["Adult"][[1]], Dens.Rec = dens.size["Recruit"][[1]], Dens.size = dens.size.c, Dens.size.reg = dens.size.c.reg)
    } else {
      
      im.covs <- listof(MDT = im.MDTS[[1]], MDS = im.MDTS[[2]],
                        Lad.var = im.lad[[1]], Lad.cv= im.lad[[2]], Lad.shan = im.lad[[3]],
                        CanOpen = im.Light[[1]], LAI = im.Light[[2]], DirectBelow = im.Light[[3]], DiffBelow = im.Light[[4]], DirectBelow.Year = im.Light[[5]], DiffBelow.Year = im.Light[[6]],
                        Canopy = im.Canopy[[1]], 
                        Hed = im.Understory[[1]], Pter = im.Understory[[2]], Rub = im.Understory[[3]], Scl = im.Understory[[4]],
                        Dens = dens.all, Dens.Adult = dens.size["Adult"][[1]], Dens.Rec = dens.size["Recruit"][[1]], Dens.size = dens.size.c, Dens.size.reg = dens.size.c.reg)
      
    }
    
    if (length(im.Understory) == 5) {
      
      fit.env8 <- ppm(unmark(data.sp.reg) ~ 1 + x + y + MDT + MDS + Lad.var + Lad.cv + Lad.shan + Canopy + CanOpen + LAI + DiffBelow + DiffBelow + DirectBelow.Year + Fs + Hed + Pter + Rub + Scl, data = im.covs)
      
    } else {
      
      fit.env8 <- ppm(unmark(data.sp.reg) ~ 1 + x + y + MDT + MDS + Lad.var + Lad.cv + Lad.shan + Canopy + CanOpen + LAI + DiffBelow + DiffBelow + DirectBelow.Year + Hed + Pter + Rub + Scl, data = im.covs)
      #You can use Gam using the following example bs(MDT, 5) and use.gam=TRUE
    }
    
    fit.env.mix <- ppm(unmark(data.sp.reg) ~ 1 + Canopy, data = im.covs)
    fit.env <- step(fit.env8, trace = 0)
    title.env <- fit.env$trend[[2]]
    
    layout(matrix(c(1,1,2,3), 2, 2, byrow = FALSE), heights = c(2, 2), widths = c(2, 2))
    par(mar = c(0, 0, 3, 4))
    plot(predict(fit.env, type = "intensity"),  main = NULL)
    plot(data.sp.reg, add = TRUE, cex = 1, col="grey")
    fit.env.sp <- data.frame(Sp = data.sp.reg$marks, pred = fitted(fit.env, dataonly=T))
    fit.mix.sp <- data.frame(Sp = data.sp.reg$marks, pred = fitted(fit.env.mix, dataonly=T))
    par(mar = c(0, 0, 3, 4))
    boxplot(pred ~ Sp, data = fit.env.sp, outline=F)
    par(mar = c(0, 0, 3, 4))
    boxplot(pred ~ Sp, data = fit.mix.sp, outline=F)
    title(paste0("Plot A", i, " // Best model: ", unlist(title.env)), line = -1, cex.main = 1, outer = TRUE)
    old.par
    
    L.E.env <- envelope(fit.env, Lest, nsim=nsim, fix.n=TRUE, correction="Ripley")
    g.E.env <- envelope(fit.env, pcf, nsim=nsim, fix.n=TRUE, correction = "Ripley")
    kNN.E.env <- envelope(fit.env, Gest, nsim=nsim, fix.n=TRUE, correction="rs")
    F.E.env <- envelope(fit.env, Fest, nsim=nsim, fix.n=TRUE, correction="rs")
    
    par(mfrow=c(2,2), mar=c(1, 1, 1.25, 1.25), oma = c(4, 4, 2, 2)) 
    plot(L.E.env, . -r ~ r, shade=c("hi", "lo"), main = NULL, legend = F)
    plot(g.E.env, main = NULL, legend = F)
    plot(kNN.E.env, main = NULL, legend = F)
    plot(F.E.env, main = NULL, legend = F)
    title(paste("Plot", i), line = 0, outer = TRUE)
    old.par
    
    
    # Density covariates
    
    if (!is.null(dens.size["Adult"][[1]])) {
      
      fit.dens2 <- ppm(as.formula(paste0("unmark(data.sp.reg) ~ ", formula(fit.env), "+ Dens.Adult + Dens.Rec + Dens.size.reg")[[2]]), data = im.covs)
    
    } else {
      
      fit.dens2 <- ppm(as.formula(paste0("unmark(data.sp.reg) ~ ", formula(fit.env), "+ Dens.Rec + Dens.size.reg")[[2]]), data = im.covs)
      
    }
    
    fit.dens <- step(fit.dens2, trace = 0)
    title.dens <- fit.dens$trend[[2]]
    
    layout(matrix(c(1,1,2), 1, 3, byrow = FALSE), heights = c(2, 2), widths = c(2, 2))
    par(mar = c(0, 0, 3, 4))
    plot(predict(fit.dens, type = "intensity"), main = NULL)
    plot(data.sp.reg, add = TRUE, cex =  1, col = "grey")
    fit.dens.sp <- data.frame(Sp = data.sp.reg$marks, pred = fitted(fit.dens, dataonly=T))
    par(mar = c(0, 0, 2.5, 0))
    boxplot(pred ~ Sp, data = fit.dens.sp, outline=F)
    title(paste0("Plot A", i, " // Best model: ", unlist(title.dens)), line = -1, cex.main = 1, outer = TRUE)
    old.par
    
    L.E.dens <- envelope(fit.dens, Lest, nsim=nsim, fix.n=TRUE, correction="Ripley")
    g.E.dens <- envelope(fit.dens, pcf, nsim=nsim, fix.n=TRUE, correction = "Ripley")
    kNN.E.dens <- envelope(fit.dens, Gest, nsim=nsim, fix.n=TRUE, correction="rs")
    F.E.dens <- envelope(fit.dens, Fest, nsim=nsim, fix.n=TRUE, correction="rs")
    
    par(mfrow=c(2,2), mar=c(1, 1, 1.25, 1.25), oma = c(4, 4, 2, 2)) 
    plot(L.E.dens, . -r ~ r, shade=c("hi", "lo"), legend = F, main=NULL)
    plot(g.E.dens, legend = F, main=NULL)
    plot(kNN.E.dens, legend = F, main=NULL)
    plot(F.E.dens, legend = F, main=NULL)
    old.par
    
    
    
    # Test Cluster (Thomas) process -------------------------------------------
    
    fit.t0 <- kppm(unmark(data.sp.reg) ~ 1, "Thomas", statistic="pcf") #Homogeneous
    plot(listof(Observed = fit.t0$X, simulate(fit.t0, 7)), main = paste0("Thomas / R (scale): ", round(fit.t0$modelpar[2], 3), " / cl. size (mu): ",round(fit.t0$mu, 3)))
    
    L.E.t <- envelope(fit.t0, Lest, nsim = nsim, correction = "Ripley")
    g.E.t <- envelope(fit.t0, pcf, nsim=nsim, fix.n=TRUE, correction = "Ripley")
    kNN.E.t <- envelope(fit.t0, Gest, nsim=nsim, fix.n=TRUE, correction="rs")
    F.E.t <- envelope(fit.t0, Fest, nsim=nsim, fix.n=TRUE, correction="rs")
    
    par(mfrow=c(2,2), mar=c(0.5, 0.5, 0.2, 0.2), oma = c(4, 4, 0.2, 0.2))
    plot(L.E.t, . -r ~ r, shade=c("hi", "lo"), main=paste("Plot", i, "// Thomas process"), legend = F)
    plot(g.E.t, legend = F)
    plot(kNN.E.t, legend = F)
    plot(F.E.t, legend = F)
    old.par
    
    
    
    # Test effect of marks ----------------------------------------------------
    
    # http://rosanaferrero.blogspot.com/2010/12/procesos-espaciales-puntuales-con-r.html
    # https://www.rdocumentation.org/packages/ads/versions/1.5-3/topics/k12fun
    # Change markconnect by Jdif
    
    Jdif <- function(X, ..., i) {
      Jidot <- Jdot(X, ..., i = i)
      J <- Jest(X, ...)
      dif <- eval.fv(Jidot - J)
      return(dif)
    }
    
    rlab.E.size.c.ad <- envelope(data.size.c.ad, markcorr, nsim = nsim, envelope = TRUE) 
    
    fit.size.ad <- ppm(unmark(data.size.c.ad) ~ 1, DiggleGratton(r=2))
    cif.size.ad <- predict(fit.size.ad, type = "cif", ngrid=512)
    cif.size.ad.R <- residuals(fit.size.ad)
    
    layout(matrix(c(1,1,2,3), 2, 2, byrow = FALSE), heights = c(2, 2), widths = c(2, 2))
    plot(cif.size.ad.R, legend =F)
    
    par(mar = c(0, 7, 2.5, 0))
    plot(rlab.E.size.c.ad, legend = F)
    
    fit.size.ad.sum <- data.frame(pred = fitted(fit.size.ad, dataonly=T), Size = data.size.c.ad$marks)
    reg1 <- lm(Size~pred,data=fit.size.ad.sum) 
    
    par(mar = c(0, 7, 3, 0))
    plot(fit.size.ad.sum)
    abline(reg1, col="red")
    old.par
    
    # diagnose.ppm(fit.size.ad, type="Pearson", envelope=TRUE, nsim=10)
    
    rlab.E.size.c.reg <- envelope(data.size.c.reg, markcorr, nsim = nsim, envelope = TRUE) 
    
    layout(matrix(c(1,2,2), 1, 3, byrow = FALSE), heights = c(2, 2), widths = c(2, 2))
    plot(rlab.E.size.c.reg, legend = F)
    
    fit.size.reg <- ppm(unmark(data.size.c.reg) ~ 1, DiggleGratton(r=0.2))
    cif.size.reg <- predict(fit.size.reg, type = "cif", ngrid=512)
    cif.size.reg.R <- residuals(fit.size.reg)
    plot(cif.size.reg.R)
    old.par
    
    rlab.E.reg <- envelope(data.sp.reg, markconnect, nsim=nsim)
    rlab.E.sp.reg <- alltypes(data.sp.reg, markconnect, nsim = nsim, envelope = TRUE, simulate=expression(rlabel(data.sp.reg)))
    plot(rlab.E.sp.reg, legend = F, transpose_= T, mar.panel=c(1, 1, 1, 1))
    
    rlab.E.size.c.reg <- envelope(data.size.c.reg, markcorr, nsim = nsim, envelope = TRUE) 
    plot(rlab.E.size.c.reg, legend = F)
    
    #rlab.E.size.reg <- alltypes(data.size.reg, markconnect, nsim = nsim, envelope = TRUE, simulate=expression(rlabel(data.size.reg)))
    #plot(rlab.E.size.reg, legend = F, transpose_= T, mar.panel=c(1, 1, 1, 1))
    
    if (j != length(Year)){
      
      layout(matrix(c(1,1,2,3), 2, 2, byrow = FALSE), heights = c(2, 2), widths = c(2, 2))
      par(mar = c(3, 3, 3, 4))
      plot(data.fate.reg.i,  leg.side = "bottom", main = NULL, bg = 19, cols = c("gray", "red", "black"))
      rlab.E.fate.reg <- envelope(data.fate.reg, markconnect, nsim = nsim, simulate=expression(rlabel(data.fate.reg)))  
      par(mar = c(3, 0, 3, 3))
      plot(rlab.E.fate.reg, main= NULL, legend = F)
      
      rlab.E.fate.reg.i <- K012(data.fate.reg.i, fijo="Tree", i="RcDead", j="RcSurv",
                                r=seq(0,8,le=51), nsim=nsim, nrank=5, correction="isotropic")
      par(mar = c(3, 0, 3, 3))
      plot(rlab.E.fate.reg.i$k01, sqrt(./pi)-r~r,  col=c(1, 2, 1), lty=c(3, 1, 3), las=1,
           ylab=expression(L[12]), xlim=c(0, 8), legend = F,
           main=NULL)
      title(paste0("Plot A", i, " // Distance-dependent fate of recruits"), line = -1, cex.main = 1, outer = TRUE)
      old.par
      
    }
    
    if (Year[j] == 2009){
      
      layout(matrix(c(1,1,1, 2, 3, 4), 3, 2, byrow = FALSE), heights = c(1, 1), widths = c(1, 1))
      par(mar = c(3, 3, 3, 4))
      
      plot(dens.size[["Adult"]])
      plot(data.growth.ad,  leg.side = "bottom", pch=21, cols = ifelse(as.character(data.fate.ad$marks)=="0", "black", "grey"), main = NULL, add=T)
      rlab.E.growth.ad <- envelope(data.growth.ad, markcorr, nsim=nsim)
      
      par(mar = c(3, 3, 3, 4))
      plot(rlab.E.growth.ad, legend = F)
      
      if (sum((data.fate.ad$marks) == 0) > 0){ #if there is any dead tree
        rlab.E.fate.ad <- envelope(data.fate.ad, markconnect, nsim = nsim, simulate=expression(rlabel(data.fate.ad)))  
        plot(rlab.E.fate.ad, legend = F)
      }
      title(paste0("Plot A", i, " // Growth and fate of trees"), line = -1, cex.main = 1, outer = TRUE)
      
      fit.growth.c.sum <-  data.frame(Density = dens.size[["Adult"]][data.growth.ad], Growth = data.growth.ad$marks)
      reg1 <- lm(Density~Growth,data=fit.growth.c.sum)
          
      par(mar = c(0, 3, 3, 0))
      plot(fit.growth.c.sum)
      abline(reg1, col="red")
      
      old.par
    }
    
  }
  
  dev.off()
  
}

