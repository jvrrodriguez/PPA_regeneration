library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)

source("~/Documentos/Datos NO publicados/BioIntForest/PPA_regeneration/function quantumPlot.R")

Year <- c(2008, 2009, 2010, 2011)
data <- c("g.E.rec", "kNN.E.rec", "Jdif.E.rec","g.E.ac", "Jdif.E.fate.rec", "K012.E.fate.rec.i", "g.E.env", "g.E.dens")
y_stats <- list("g(r)", "D(r)", "Ji· (r) - J(r)", "g(r))", "Ji· (r) - J(r)", expression(K[012] (r)), "g(r)", "g(r)")
log.data <- c(T, F, F, T, F, F, T, T)
y_inter <- c(1, 0, 0, 1, 0, 0, 1, 1)
save.output <- T


for (i in 1:length(data)) {

  data.cat <- NULL

  if (save.output == T) pdf(paste0("~/Documentos/Datos NO publicados/BioIntForest/PPA_Results/Figures/PPA_", data[i],".pdf"), width = 7, height = 5)  
  
  if (grepl("fate", data[i], fixed = TRUE)) {
    
    year.tmp <- length(Year) - 1
    
  } else {
    
    year.tmp <- length(Year)
    
  }
    
  for (j in 1:year.tmp) {
  
    load(file = paste0("~/Documentos/Datos NO publicados/BioIntForest/PPA_Results/Reports/PPA_", Year[j], ".RData"))
    
    data.ctrl <- get(paste0(data[i],".ctrl"))
    data.thnn <- get(paste0(data[i],".thnn"))
    
    if (j == 1) {
      
      plot.ctrl <- data.ctrl
      plot.thnn <- data.thnn
      
    }
      
    data.cat <- rbind(data.cat, data.frame(Year = Year[j], Treat = "Ctrl", r = data.ctrl$r, value = ppp.cat(data.ctrl)),
                      data.frame(Year = Year[j], Treat = "Thnn", r = data.thnn$r, value = ppp.cat(data.thnn)))
  
  }
  
  Treat <- unique(data.cat$Treat)
  Plot.list <- list()
  
  for (j in Treat) { 
    
    Plot.list[[j]] <- ggplot(data.cat[data.cat$Treat == j,], aes(x = r, y = Year)) + 
      geom_raster(aes(fill = value)) + 
      scale_fill_gradient2(low = "indianred3", mid = "white", high = "darkolivegreen3", midpoint = .0) +
      labs(x = "Distance r (m)", y = "Year") +
      scale_y_continuous(trans = "reverse") + 
      theme_bw() + theme(axis.text.x = element_text(size = 9, angle = 0, vjust = 0.3),
                         axis.text.y = element_text(size = 7, angle = 90),
                         plot.title = element_text(size = 11)) + theme(legend.position = "none")
    
    
  }
  
  #https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html
  grid.arrange(quantumPlot(plot.ctrl, log.y = log.data[i], y_intercept = y_inter[i], y_name = y_stats[[i]], quantum = F), 
               quantumPlot(plot.thnn, log.y = log.data[i], y_intercept = y_inter[i], y_name = y_stats[[i]], quantum = F), 
               Plot.list$Ctrl, Plot.list$Thnn, 
               top = paste(data[i], "//", Treat[1], "vs", Treat[2]), widths = c(2,2), heights = c(1.5, 1), nrow = 2) #, colour=c("indianred2", "grey", "white")
  
  if (save.output == T) dev.off()

}
