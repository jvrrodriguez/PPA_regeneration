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

# Algorithms for automatically finding appropriate thresholds for numerical data
# Based on ImageJ sowtware
library(autothresholdr)

# library(shadow) # SR Package for GeometricShadow Calculations in an UrbanEnvironment 
# https://cran.r-project.org/web/packages/shadow/vignettes/introduction.pdf

# extract metadata from files, The more sophisticated is exiftool used by exiftoolr package
#library(exifr) 
#library(EXIFr)
#library(exiftoolr)

#Image segmentation
#https://imagej.net/_images/8/87/Arganda-Carreras-Segmentation-Bioimage-course-MDC-Berlin-2016.pdf
# https://scholar.google.es/citations?hl=es&user=3XVchKQAAAAJ&view_op=list_works&sortby=pubdate


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
location.threshold  =  0.65
location.magnetic <- c(90-55, 90-49, 90-40, 90-50, 90-51, 90-60, 90-47, 90-49, 90-40)   #reorientation for each image; now varies for each site but it can slightly varies for each picture

### atmospheric parameters (Calculados el 22/jun/2013 un dia despejado en Aspurz); un dia nuboso podria ser 0.6
location.tau        = 0.9           # Clear-sky transmission coefficient. Varies tipycally from 0.6 to 0.7 in dust-free regions
location.uoc        = 0.1           # Sky-region brightness

#######END          initialize site and image data          ##########





##########          load image namess          ##########

### We assume colour images
### in a subdirectory images
### load all JPG file names in a list

folder <- paste("~/Documentos/Datos NO publicados/BioIntForest/Data/Fotos_hemisfericas/Aspurz_2009/")
Year <- "2009" #para nombrar el archivo de salida

all.images = list.files(folder, pattern = ".JPG")
nr.images = length(all.images); nr.images
  
## Create data frame to hold all results
all.data = data.frame(matrix(0, nr.images, 16))
names(all.data) = c("File", "Plot", "threshold", 
                      "CanOpen", "LAI",
                      "DirectAbove", "DiffAbove", "DirectBelow", "DiffBelow", 
                      "N.Sunflecks", "Mdn.Sunflecks", "Max.Sunflecks",
                      "DirectAbove.Yr", "DiffAbove.Yr", "DirectBelow.Yr", "DiffBelow.Yr")
all.data[,1] = all.images
  
#if (save.output == T) pdf(paste0("folder, "Hemiphot_2008.pdf"), width=7, height=5, compress = TRUE)
  
## now the batch can start
t1 = Sys.time()
  
for(i in 1:nr.images) {    
    
    ## read file
    
    name.tmp <- strsplit(all.images[i], "_")[[1]][[1]]  
    all.data[i,2] = strsplit(name.tmp, "A ")[[1]][[2]]
    
    #image.data <- exif_read(paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/Fotos hemisféricas espaciales Aspurz 2004/",Plot[i]," 2004/",all.images[i],sep = ""))
    #location.day <- image.data$ModifyDate

    
    # Automatic detection of canopy-sky threshold -----------------------------

    # Read image & calculate sky-canopy threshold
    
    image.threshold <- location.threshold
    
    if (is.null(location.threshold)) {
      
      image <- readImage(paste0(folder, all.images[i],sep = ""))
      image.threshold <- Ridler(image,p = TRUE) #(ISodata algorithm in ImageJ) Este algoritmo generalmente da valores bastante bajos, por lo que se considera buena parte del dosel como vegetación
      all.data[i,3] = image.threshold
      
    } else {
      
      all.data[i,3] = image.threshold
      
    }
    
    
    
    # Compute Luminic indexes -------------------------------------------------
    
    ## convert to Hemi image
    
    image = readJPEG(paste0(folder, all.images[i],sep = ""), native = F)     #if native = T creates a raster, else an array
    
    par(mfrow = c(1, 2))
    
    image = Image2Hemiphot(image)
    PlotHemiImage(image = image, draw.circle = T)
    
    arrows(x0 = dim(image[[1]])[[2]]/2, y0 = dim(image[[1]])[[1]]/2, x1 = dim(image[[1]])[[2]], y1 = dim(image[[1]])[[2]]/2 + dim(image[[1]])[[2]] * sin(-location.magnetic[i]), col="red")
    
    ## set cirlce parameters
    image = SetCircle(image, cx = location.cx, cy = location.cy)
    #PlotHemiImage(image)
    
    imageB <- SelectRGB(image, "B")
    # imageR <- SelectRGB(image, "R")
    # imageG <- SelectRGB(image, "G")
    # 
    # imageR[[1]] <- 20*exp(imageR[[1]])/max(20*exp(imageR[[1]]))
    # imageG[[1]] <- 20*exp(imageG[[1]])
    PlotHemiImage(imageB)
    
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
                      d = location.day, magn.corr = location.magnetic[i], sun.location = T, h = k)
    }
    
    par(mfrow = c(1,2))
    
    #PlotPAR.Day(radiation = rad, real.time = T)
    
    Rad = CalcPAR.Year(imageB, lat = location.latitude, tau = location.tau, uoc = location.uoc, magn.corr = location.magnetic[i])
    
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
    if (is.null(dim(rad.max.red))) rad.max.red <- matrix(rad.max.red, nrow = 1)
    
    event <- NULL
    cnt <- 0
    
    if (dim(rad.max.red)[1] > 2) {
    
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
      
      if (!is.na(sun.time)) hist(sun.time)
        
    } else {
      
      sun.time <- 0
      
    }
    
      
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
    
    
    ## save data

    if (file.exists(paste0(folder, "Hemiphot_Plots_", Year, ".csv"))) {
      
      read.table(paste0(folder, "Hemiphot_Plots_", Year, ".csv"), sep = ",")
      
    }
    
    if (save.output == T) write.table(all.data, paste0(folder, "Hemiphot_Plots_", Year, ".csv"), sep = ",")
    
  }
  
t2 = Sys.time()
  
##time per image
(t2 - t1)/nr.images
  
#if (save.output == T) dev.off()
