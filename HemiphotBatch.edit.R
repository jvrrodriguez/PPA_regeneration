##########          Hemiphot.R         ##########
# a script to calculate light indeces and more from hemispherical images
# read the manual and helfile to use this script.
# HemiphotTest.R provides examples of all single calculations
# HemiphotBatch.R provides examples of batch calculations for a 
# number of images in one directory
#
# to be able to use Hemophot.R yo will need to install the functions
# with the command >source("Hemiphot.R")    
# provided the script is foudn in our working directory 
#######END          Hemiphot.R         ##########

### batch processing in terminal
# cd /
# cd /home/javier/Documentos/Datos\ NO\ publicados/BioIntForest/Analysis/
# R CMD BATCH HemiphotBatch.edit.R



##########          How to cite the use of Hemiphot         ##########
#If you use Hemiphot.R for your research or work please cite as:
#Hans ter Steege (2018) Hemiphot.R: Free R scripts to analyse hemispherical 
#photographs for canopy openness, leaf area index and photosynthetic 
#active radiation under forest canopies.  
#Unpublished report. Naturalis Biodiversity Center, Leiden, The Netherlands
#https://github.com/Hans-ter-Steege/Hemiphot
#######END          How to cite the use of Hemiphot         ##########





### This is the batch script for Hemiphot.R
### Here you can run all functions and store  
### all data by file


### clear memory and set working directory
rm(list = ls())




##########          load libraries and source files          ##########

library(jpeg)      			     # read and write jpg's

library(ggplot2)
library(dplyr)

# https://cran.r-project.org/web/packages/imager/vignettes/gettingstarted.html
# https://dahtah.github.io/imager/imager.html
# library(imager)

# analysing canopy geometry and solar radiation regimes using hemispherical photographs
# Por ahora no calcula sunflecks (el sofwtare original, CIMES-FISHEYE, si)
# https://github.com/ggranath/cimesr #http://jmnw.free.fr/
library("devtools") #devtools::install_github("cimesr", "ggranath", dependencies=TRUE) 

#install.packages("BiocManager")
#BiocManager::install("EBImage")
library(EBImage)
#https://www.bioconductor.org/packages/release/bioc/vignettes/EBImage/inst/doc/EBImage-introduction.html
library(Sky) #Canopy Openness Analyzer Package (needs EBImage)

# CAnopy IMage ANalysis
# To correct images affected by sunlight canopy, multiple scattering, vignetting, blooming, and chromatic aberration
# https://github.com/GastonMauroDiaz/caiman
# devtools::install_github("GastonMauroDiaz/caiman")
library(caiman)

# BIO7 Former ImageJ software. 
# Allows to calclate many other binarization (threshold) indexes to separate canopy-sky
# https://bio7.org/screenshots/

# Algorithms for automatically finding appropriatethresholds for numerical data
# Based on ImageJ sowtware
library(autothresholdr)

# library(shadow) # SR Package for GeometricShadow Calculations in an UrbanEnvironment 
# https://cran.r-project.org/web/packages/shadow/vignettes/introduction.pdf

# extract metadata from files, The more sophisticated is exiftool used by exiftoolr package
#library(exifr) 
#library(EXIFr)
#library(exiftoolr)

source("~/Documentos/Datos NO publicados/BioIntForest/Analysis/Hemiphot.R")  # functions to calculate all results

Detect.Sunspots <- FALSE
Improve.images <- FALSE
save.output <- TRUE


#######END          load libraries and source files          ##########



##########          initialize site and image data          ##########

### Location parameters
location.latitude   = 42.70861111111111
location.altitude   = 642
location.longitude  = 1.1444444444444444
location.day        = 275 # 1st of october; no se exactamente el dia pero por lo que conton Bosco, deberia ser septiembre o octubre.
location.days       = seq(15,360,30)   # roughly each mid of the 12 months

#location <- SpatialPointsDataFrame(coords = data.frame(x = 42.70861111111111, y = 1.1444444444444444), data = location, proj4string = CRS("+proj=longlat +datum=WGS84"))

### Image parameters
## determine in Hemiphot.R and fill in here for batch processing
location.cx         = 0             # x coordinate of center
location.cy         = 0             # y coordinate of center
location.cr         = 0             # radius of circle
#location.threshold  = 0.65
location.magnetic <- c(90-55, 90-49, 90-40, 90-50, 90-51, 90-60, 90-47, 90-49, 90-40)   #reorientation for each image; now varies for each site but it can slightly varies for each picture

### atmospheric parameters (Calculados el 22/jun/2013 un dia despejado en Aspurz); un dia nuboso podria ser 0.6
location.tau        = 0.9           # Clear-sky transmission coefficient. Varies tipycally from 0.6 to 0.7 in dust-free regions
location.uoc        = 0.1           # Sky-region brightness

#######END          initialize site and image data          ##########





##########          load image namess          ##########

### We assume colour images
### in a subdirectory images
### load all JPG file names in a list

Plot <- paste0("A",1:9)
Plot.data <- NULL
# setwd("~/Documentos/Datos NO publicados/BioIntForest/Data")
# Plot.data <- read.csv("Hemiphot_Plot.csv")

for (j in 1:length(Plot)) {
  
  all.images = list.files(paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/Fotos_hemisfericas/Fotos hemisféricas espaciales Aspurz 2004/Fotos hemisféricas Aspurz 2004. Red fina previa a corta/",Plot[j]," 2004/"),pattern = ".JPG")
  nr.images = length(all.images); nr.images
  
  ## Create data frame to hold all results
  all.data = data.frame(matrix(0, nr.images, 16))
  names(all.data) = c("File", "Plot", "SunSpot",
                      "CanOpen", "LAI",
                      "DirectAbove", "DiffAbove", "DirectBelow", "DiffBelow", 
                      "N.Sunflecks", "Mdn.Sunflecks", "Max.Sunflecks",
                      "DirectAbove.Yr", "DiffAbove.Yr", "DirectBelow.Yr", "DiffBelow.Yr")
  all.data[,1] = all.images
  
  if (save.output == T) pdf(paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/Fotos_hemisfericas/Hemiphot_A", j,"_2004.pdf"), width=7, height=5, compress = TRUE)
  
  ## now the batch can start
  t1 = Sys.time()
  
  for(i in 1:nr.images){    
    
    ## read file
    
    all.data[,2] = Plot[j]
    
    #image.data <- exif_read(paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/Fotos hemisféricas espaciales Aspurz 2004/",Plot[j]," 2004/",all.images[i],sep = ""))
    #location.day <- image.data$ModifyDate

    
    
    # Detect sun spot in pictures  --------------------------------------------
    
    if (Detect.Sunspots == TRUE) {
      
      img = load.image(paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/Fotos_hemisfericas/Fotos hemisféricas espaciales Aspurz 2004/",Plot[j]," 2004/",all.images[i],sep = ""))
      #plot(img)
      
      # bdf <- as.data.frame(img)
      # bdf <- mutate(bdf,channel=factor(cc,labels=c('R','G','B')))
      # p_white <- c(length(bdf$value[bdf$channel == "R" & bdf$value == 1]) / length(bdf$value[bdf$channel == "R"]),
      #   length(bdf$value[bdf$channel == "G" & bdf$value == 1]) / length(bdf$value[bdf$channel == "G"]),
      #   length(bdf$value[bdf$channel == "B" & bdf$value == 1]) / length(bdf$value[bdf$channel == "B"]))
      # p_white
      # 
      # ggplot(bdf,aes(value,col=channel))+geom_histogram(bins=30)+facet_wrap(~ channel)
  
      imgChT <- imsplit(img, "c")
      imgChT[[1]] <- threshold(imgChT[[1]],"99%")
      imgChT[[2]] <- threshold(imgChT[[2]],"99%")
      imgChT[[3]] <- threshold(imgChT[[3]],"99%")
      #plot(imgChT[[1]])
      
      # Select the first channel. Red??
      imgB <- isoblur(imgChT[[1]], 100) #very high blurrness
      
      
      #AUTOMATIC THRESHOLDING ENTRE VEGETACION Y CIELO PARA CADA FOTOGRAFIA, SERIA LO MAS CORRECTO
      #ALGUNO DE LOS PAQUETES LO HACE? sky?)
      
      Hdet <- with(imhessian(imgB),(xx*yy - xy^2))
      #plot(Hdet,main="Determinant of Hessian")
      
      threshold(Hdet,"99.9%") %>% plot(main="Determinant: 1% highest values")
      
      lab <- threshold(Hdet,"99.9%") %>% label #get a strongest threshold
      plot(lab,main="Labelled regions")
      
      df <- as.data.frame(lab) %>% subset(value>0)
      
      all.data[,3] <- ifelse(identical(unique(df$value), numeric(0)) == TRUE, 0, 1)
      
      df <- as.data.frame(lab) %>% subset(value>0)
      centers <- dplyr::group_by(df,value) %>% dplyr::summarise(mx=mean(x),my=mean(y))
      plot(img)
      with(centers,points(mx,my,col="red"))
    
    } else {
      
      all.data[,3] <- NA
    
    }
    
    
    
    # Improve images with caiman -----------------------------
    
    if (Improve.images == T) {
      
      # see the help of loadPhoto with ?loadPhoto. 
      # Maybe you need to use the arguments upperLeft, width and height.
      x <- loadPhoto(paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/Fotos_hemisfericas/Fotos hemisféricas espaciales Aspurz 2004/",Plot[j]," 2004/",all.images[i],sep = ""))
      x <- crop(x, extent(200, 1400, 0, 1200))
      x <- as(x, "CanopyPhoto")
      fisheye(x) <- newFishEye(TRUE, TRUE, FALSE)
      x <- normalize(x, 0, 255)
      # Here, I assume that your lens has perfect polor projection.
      z <- makeZimage(nrow(x), lensPolyCoef())
      m <- doMask(z)
      bin <- autoThr(enhanceHP(x, sharpen = TRUE))
      # This takes a while but you can see the progres
      seg <- doPolarQtree(x, z, scaleParameter = 0.2)
      name <- strsplit(i, "\\.")[[1]][1]
      out <- doOBIA(x, bin, z, seg, zlim = asAngle(c(20, 80)))
      writeRaster(out * 255, paste0(path_out, "/", name, ".TIF"), datatype = "INT1U", overwrite = TRUE)
      
    } 
  
    

    # Automatic detection of canopy-sky threshold -----------------------------

    # Read image & calculate sky-canopy threshold
    
    image <- readImage(paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/Fotos_hemisfericas/Fotos hemisféricas espaciales Aspurz 2004/",Plot[j]," 2004/",all.images[i],sep = ""))
    image.threshold <- Ridler(image,p = TRUE) #(ISodata algorithm in ImageJ)
    
    # image <- readJPEG(paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/Fotos_hemisfericas/Fotos hemisféricas espaciales Aspurz 2004/",Plot[j]," 2004/",all.images[i],sep = ""))
    # writeImage(image, "~/Descargas/image_tmp.tiff", type = "tiff")
    # image <- ijtiff::read_tif("image_tmp.tiff")
    # mask <- auto_thresh_mask(image, "Isodata")
    # image.threshold <- auto_thresh(image, "Otsu")

    
    
    # Compute Luminic indexes -------------------------------------------------
    
    ## convert to Hemi image
    
    image = readJPEG(paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/Fotos_hemisfericas/Fotos hemisféricas espaciales Aspurz 2004/",Plot[j]," 2004/",all.images[i],sep = ""), native = F)     #if native = T creates a raster, else an array
    
    par(mfrow = c(1, 2))
    
    image = Image2Hemiphot(image)
    PlotHemiImage(image = image, draw.circle = T)
    
    arrows(x0 = dim(image[[1]])[[2]]/2, y0 = dim(image[[1]])[[1]]/2, x1 = dim(image[[1]])[[2]], y1 = dim(image[[1]])[[2]]/2 + dim(image[[1]])[[2]] * sin(-location.magnetic[j]), col="red")
    
    ## set cirlce parameters
    image = SetCircle(image, cx = location.cx, cy = location.cy)
    #PlotHemiImage(image)
    
    imageB <- SelectRGB(image, "B")
    # imageR <- SelectRGB(image, "R")
    # imageG <- SelectRGB(image, "G")
    # 
    # imageR[[1]] <- 20*exp(imageR[[1]])/max(20*exp(imageR[[1]]))
    # imageG[[1]] <- 20*exp(imageG[[1]])
    # PlotHemiImage(imageR)
    
    #select each channel and apply a threshold
    imageB = ThresholdImage(im = imageB, th = image.threshold, draw.image = T)
    #imageR = ThresholdImage(im = imageR, th = image.threshold, draw.image = F)
    #imageG = ThresholdImage(im = imageG, th = image.threshold, draw.image = F)
    
    # canopy openness
    gap.fractions = CalcGapFractions(imageB)
    all.data[i,4] = CalcOpenness(fractions = gap.fractions)
    
    ## calculate LAI according to Licor's LAI Analyzer 
    all.data[i,5] = CalcLAI(fractions = gap.fractions)
    
    ## Photosynthetic Photon Flux Density (PPDF, umol m-1 s-1) P
    rad.avg = CalcPAR.Day(im = imageB,
                      lat = location.latitude, 
                      lon = location.longitude, 
                      d = location.days,
                      tau = location.tau, uoc = location.uoc, 
                      draw.tracks = F, full.day = F)
    
    PlotHemiImage(imageB, draw.circle = T)
    
    ### plot location.angle
    mtext(all.data[i,1], side = 3, line = -3, outer = TRUE)
    
    for(k in 3:21){
      DrawSolarTracks(imageB, lat = location.latitude, lon = location.longitude, time.zone = +1,
                      d = location.day, magn.corr = location.magnetic[j], sun.location = T, h = k)
    }
    
    par(mfrow = c(1,2))
    
    #PlotPAR.Day(radiation = rad, real.time = T)
    
    Rad = CalcPAR.Year(imageB, lat = location.latitude, tau = location.tau, uoc = location.uoc, magn.corr = location.magnetic[j])
    
    PlotPAR.Year(radiation = Rad)
    abline(v=c(172, 264), col="red", lty= 2)
    
    mtext(all.data[i,1], side = 3, line = -3, outer = TRUE)
    
    rad.max = CalcPAR.Day(im = imageB,
                          lat = location.latitude, 
                          lon = location.longitude, 
                          d = which(Rad[,2] == max(Rad[,2])),
                          tau = location.tau, uoc = location.uoc, 
                          draw.tracks = T, full.day = T)
    
    PlotPAR.Day(radiation = rad.max, real.time = F)
    
    
    ## Compute sunflecks and some summary statistics
    
    rad.max.red <- rad.max[rad.max[,"DirectUnder"] > 0, ]
    
    event <- NULL
    cnt <- 0
    
    for (k in 2:nrow(rad.max.red)) {
      
      time.tmp <- rad.max.red[,"Solar Time"][k] - rad.max.red[,"Solar Time"][k-1] 
      if (time.tmp > 0.04) cnt <- cnt + 1
      event <- c(event, cnt)
      
    }
    
    sun.time <- array()
    
    for (k in unique(event)[-1]) {
      
      tmp <- rad.max.red[which(event == unique(event)[[k]]),]
      
      if (length(tmp) > 6){
        
        sun.time <- c(sun.time, 0.033333 * nrow(tmp))
        
      } else {
      
        sun.time <- c(sun.time, 0.033333)
        
      }
      
    }
     
    hist(sun.time)
    
    par(mfrow = c(1,1))

    rad.max = CalcPAR.Day(im = imageB,
                          lat = location.latitude, 
                          lon = location.longitude, 
                          d = which(Rad[,2] == max(Rad[,2])),
                          tau = location.tau, uoc = location.uoc, 
                          draw.tracks = F, full.day = F)
    
    all.data[i,6] = rad.max[1]
    all.data[i,7] = rad.max[2]
    all.data[i,8] = rad.max[3]
    all.data[i,9] = rad.max[4]
    all.data[i,10] = max(event)
    all.data[i,11] = median(sun.time, na.rm=T) #better summarizing not normal distribucions
    all.data[i,12] = max(sun.time, na.rm=T)
    
    all.data[i,13] = mean(Rad[172:264,2])
    all.data[i,14] = mean(Rad[172:264,3])
    all.data[i,15] = mean(Rad[172:264,4])
    all.data[i,16] = mean(Rad[172:264,5])
  }
  
  t2 = Sys.time()
  
  ##time per image
  (t2 - t1)/nr.images
  
  if (save.output == T) dev.off()
  
  
  ## save data
  Plot.data <- rbind(Plot.data, all.data)
  
  if (file.exists("~/Documentos/Datos NO publicados/BioIntForest/Data/Fotos_hemisfericas/Hemiphot_Plots_2004.csv")){
    read.table("~/Documentos/Datos NO publicados/BioIntForest/Data/Fotos_hemisfericas/Hemiphot_Plots_2004.csv", sep = ",")
  }
  
  if (save.output == T) write.table(Plot.data, paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/Fotos_hemisfericas/Hemiphot_Plots_2004.csv"), sep = ",")
  
}
