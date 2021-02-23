library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)

Year <- c(2008, 2009, 2010, 2011)
data <- c("g.E", "kNN.E", "Jdif.E.rec","Jdif.E.size", "Jdif.E.fate.rec", "K012.E.fate.rec.i", "g.E.env", "g.E.dens")
Treat <-  c("20%", "30%", "0%", "0%", "30%", "20%", "30%", "20%", "0%")

Plot.list <- list()
meanData <- NULL
save.output <- FALSE

for (i in 1:length(data)){

  data.Plot <- NULL
  
  if (save.output == T) pdf(paste0("~/Documentos/Datos NO publicados/BioIntForest/PPA_Results/Figures/PPA_", data[i],"_Plot.pdf"), width=9, height=5)  
  
  if (grepl("fate", data[i], fixed = TRUE)) {
    
    year.tmp <- length(Year) - 1
    
  } else {
    
    year.tmp <- length(Year)
    
  }
  
  for (j in 1:year.tmp){
    
    load(file=paste0("~/Documentos/Datos NO publicados/BioIntForest/PPA_Results/Reports/PPA_", Year[j],".RData"))
    
    data.cat <- get(paste0(data[i],".cat"))

    if (j == 1) {
      
      data.Plot <- NULL

    }
    
    data.Plot <- rbind(data.Plot, data.frame(Year = Year[j], Plot = 3, Treat = "Ctrl", r = data.cat[,2], value = data.cat[,6+2]),
                      data.frame(Year = Year[j], Plot = 4, Treat = "Ctrl", r = data.cat[,3], value = data.cat[,6+3]),
                      data.frame(Year = Year[j], Plot = 9, Treat = "Ctrl", r = data.cat[,6], value = data.cat[,6+6]),
                      data.frame(Year = Year[j], Plot = 2, Treat = "Thnn", r = data.cat[,1], value = data.cat[,6+1]),
                      data.frame(Year = Year[j], Plot = 5, Treat = "Thnn", r = data.cat[,4], value = data.cat[,6+4]),
                      data.frame(Year = Year[j], Plot = 7, Treat = "Thnn", r = data.cat[,5], value = data.cat[,6+5]))
    
  }

  #longData<-melt(data[,c(8,9,12, 7,10,11)]) #primero controles (abajo) y luego thinning
  #longData<-longData[longData$value!=0,]
  #longData$Var1 <- data[,1]
  
  plots <- unique(data.Plot$Plot)
  
  for (j in plots){
    
    Plot.list[[j]] <- ggplot(data.Plot[data.Plot$Plot == j, ], aes(x = r, y = Year)) + 
      geom_raster(aes(fill=value)) + 
      scale_fill_gradient2(low = "grey", mid = "white", high = "indianred2", midpoint = .0) +
      labs(x="Distance (m)", y="Year", title=paste("Plot", j, "//", Treat[j])) +
      theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                         axis.text.y=element_text(size=9),
                         plot.title=element_text(size=11))
    
  } 
  
  # meanData <- rbind(meanData, data.frame(Var1 = data[,1], value=rowMeans(data[,c(8,9,12, 7,10,11)]), Year = Year[i]))

  #https://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots
  
  grid_arrange_shared_legend <- function(..., nrow = 1, ncol = length(list(...)), position = c("bottom", "right")) {
    
    plots <- list(...)
    position <- match.arg(position)
    g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
    gl <- c(gl, nrow = nrow, ncol = ncol)
    
    combined <- switch(position,
                       "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                              legend,
                                              ncol = 1,
                                              heights = unit.c(unit(1, "npc") - lheight, lheight)),
                       "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                             legend,
                                             ncol = 2,
                                             widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
    grid.newpage()
    grid.draw(combined)
    
  }
  
  grid_arrange_shared_legend(Plot.list[[3]], Plot.list[[4]], Plot.list[[9]], Plot.list[[2]], Plot.list[[5]], Plot.list[[7]], nrow = 2, ncol = 3, position = "right")
  
  if (save.output == T) dev.off()
  
}

# ggplot(meanData, aes(x = Var1, y = Year)) + 
#   geom_raster(aes(fill=value)) + 
#   scale_fill_gradient2(low = "grey", mid = "white", high = "indianred2", midpoint = .0)+
#   labs(x="Distance (m)", y="Plot", title= "Yearly averaged") +
#   theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
#                      axis.text.y=element_text(size=9),
#                      plot.title=element_text(size=11))
