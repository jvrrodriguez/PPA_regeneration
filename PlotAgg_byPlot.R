
library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)

Year <- c(2008, 2009, 2010, 2011)
Treat <-  c("20%", "30%", "0%", "0%", "30%", "20%", "30%", "20%", "0%")

Plot.list <- list()
meanData <- NULL

for (j in 1:length(Year)){

 load(file=paste0("~/Documentos/Datos NO publicados/BioIntForest/Results/Results.PPA_", Year[j],".RData"))
  
  #L.E.cat, g.E.cat, kNN.E.cat, F.E.cat,
  data <- Jdif.E.size.cat
  longData<-melt(data[,c(8,9,12, 7,10,11)]) #primero controles (abajo) y luego thinning
  #longData<-longData[longData$value!=0,]
  longData$Var1 <- data[,1]
  
  Plot.list[[j]] <- ggplot(longData, aes(x = Var1, y = Var2)) + 
    geom_raster(aes(fill=value)) + 
    scale_fill_gradient2(low = "grey", mid = "white", high = "indianred2", midpoint = .0)+
    labs(x="Distance (m)", y="Plot", title=Year[j]) +
    theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                       axis.text.y=element_text(size=9),
                       plot.title=element_text(size=11))
  
  meanData <- rbind(meanData, data.frame(Var1 = data[,1], value=rowMeans(data[,c(8,9,12, 7,10,11)]), Year = Year[j]))

}


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

grid_arrange_shared_legend(Plot.list[[1]], Plot.list[[2]], Plot.list[[3]], Plot.list[[4]], nrow = 1, ncol = 4, position = "right")

ggplot(meanData, aes(x = Var1, y = Year)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient2(low = "grey", mid = "white", high = "indianred2", midpoint = .0)+
  labs(x="Distance (m)", y="Plot", title= "Yearly averaged") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))
