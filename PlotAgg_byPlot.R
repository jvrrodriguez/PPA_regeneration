library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)

Year <- 2007 + 1:4
data <- c("g.E.rec", "Jdif.E.rec", "g.E.ac", "markcor.E.size.c.rec", "Jdif.E.fate.rec", "K012.E.fate.rec.i", "g.E.env", "g.E.dens")
Treat <-  c("20%", "30%", "0%", "0%", "30%", "20%", "30%", "20%", "0%")
Treat2 <-  c("20%", "Thinned", "Control", "Control", "Thinned", "20%", "Thinned", "20%", "Control")

Plot.list <- list()
meanData <- NULL
save.output <- TRUE

signif.num <- function(x) {
  symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
         cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
         symbols = c("***", "**", "*", ".", "n.s"))
}

for (i in 1:length(data)) {

  data.cat <- NULL
  data.gof <- NULL
  
  if (save.output == T) pdf(paste0("~/Documentos/Datos NO publicados/BioIntForest/PPA_Results/Figures/PPA_", data[i],"_Plot.pdf"), 
                            width = 9, height = 8)  
  
  if (grepl("fate", data[i], fixed = TRUE)) {
    
    year.tmp <- length(Year) - 1
    
  } else {
    
    year.tmp <- length(Year)
    
  }
  
  for (j in 1:year.tmp) {
    
    load(file = paste0("~/Documentos/Datos NO publicados/BioIntForest/PPA_Results/Reports/PPA_", Year[j],".RData"))
    
    data.cat.tmp <- get(paste0(data[i],".cat"))
    data.gof.tmp <- get(paste0(data[i],".gof"))
    
    if (j == 1) {
      
      data.Plot <- NULL

    }
    
    data.cat <- rbind(data.cat, 
                      data.frame(Year = Year[j], Plot = 3, Treat = "Ctrl", r = data.cat.tmp[,2], value = data.cat.tmp[,6 + 2]),
                      data.frame(Year = Year[j], Plot = 4, Treat = "Ctrl", r = data.cat.tmp[,3], value = data.cat.tmp[,6 + 3]),
                      data.frame(Year = Year[j], Plot = 9, Treat = "Ctrl", r = data.cat.tmp[,6], value = data.cat.tmp[,6 + 6]),
                      data.frame(Year = Year[j], Plot = 2, Treat = "Thnn", r = data.cat.tmp[,1], value = data.cat.tmp[,6 + 1]),
                      data.frame(Year = Year[j], Plot = 5, Treat = "Thnn", r = data.cat.tmp[,4], value = data.cat.tmp[,6 + 4]),
                      data.frame(Year = Year[j], Plot = 7, Treat = "Thnn", r = data.cat.tmp[,5], value = data.cat.tmp[,6 + 5]))
    
    data.gof <- rbind(data.gof, 
                      data.frame(Year = Year[j], Plot = 3, Treat = "Ctrl", data.gof.tmp[data.gof.tmp$Plot == 3,]),
                      data.frame(Year = Year[j], Plot = 4, Treat = "Ctrl", data.gof.tmp[data.gof.tmp$Plot == 4,]),
                      data.frame(Year = Year[j], Plot = 9, Treat = "Ctrl", data.gof.tmp[data.gof.tmp$Plot == 9,]),
                      data.frame(Year = Year[j], Plot = 2, Treat = "Thnn", data.gof.tmp[data.gof.tmp$Plot == 2,]),
                      data.frame(Year = Year[j], Plot = 5, Treat = "Thnn", data.gof.tmp[data.gof.tmp$Plot == 5,]),
                      data.frame(Year = Year[j], Plot = 7, Treat = "Thnn", data.gof.tmp[data.gof.tmp$Plot == 7,]))
    
  }

  #longData<-melt(data[,c(8,9,12, 7,10,11)]) #primero controles (abajo) y luego thinning
  #longData<-longData[longData$value!=0,]
  #longData$Var1 <- data[,1]
  
  plots <- unique(data.cat$Plot)
  
  for (j in plots) {
    
    Plot.list[[j]] <- ggplot(data.cat[data.cat$Plot == j, ], aes(x = r, y = Year)) + 
      geom_raster(aes(fill = value)) + 
      
      scale_fill_gradient2(low = "indianred3", mid = "white", high = "darkolivegreen3", midpoint = .0) +
      
      geom_segment(aes(x = r.min, y = Year, xend = r.max, yend = Year, size = 0.005, alpha = 1/10), 
                   data = data.gof[data.gof$Plot == j,]) +
      
      geom_label(aes(x = r.min, y = Year, label = ""), data = data.gof[data.gof$Plot == j,], alpha = 1/10) +
      geom_label(aes(x = r.max, y = Year, label = ""), data = data.gof[data.gof$Plot == j,], alpha = 1/10) +
      
      geom_text(aes(label = paste("GoF p =", signif.num(p.value)), r.min + 0.0, size = 0.005, vjust = -0.8, hjust = "left"), data = data.gof[data.gof$Plot == j,]) +
      
      labs(x = "Distance (m)", y = "Year", title = paste("Plot", j, "//", Treat2[j])) +
      theme_bw() + theme(axis.text.x = element_text(size = 9, angle = 0, vjust = 0.3),
                         axis.text.y = element_text(size = 9),
                         plot.title = element_text(size = 11))
    
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
