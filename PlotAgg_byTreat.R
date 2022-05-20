library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)

source("~/Documentos/Datos NO publicados/BioIntForest/PPA_regeneration/FunsAgg.R")

Year <- 2007 + 1:4
data <- c("g.E.rec", "Jdif.E.rec", "g.E.ac", "markcor.E.size.c.rec", "Jdif.E.fate.rec", "Kmulti.E.fate.rec.i", "g.E.env", "g.E.dens")
y_stats <- list("g(r)", "Ji· (r) - J(r)", "g(r)", expression(K[mm] (r)), "Ji· (r) - J(r)", expression(K[012] (r)), "g(r)", "g(r)")
log.data <- c(T, F, T, F, F, F, T, T)
y_inter <- c(1, 0, 1, 1, 0, 0, 1, 1)
save.output <- T

for (i in 1:length(data)) {
  
  data.cat <- NULL
  data.gof <- NULL
  
  if (save.output == T) pdf(paste0("~/Documentos/Datos NO publicados/BioIntForest/PPA_Results/Figures/PPA_", data[i],".pdf"), 
                            width = 8, height =  ifelse(i != 5 & i != 6, 7, 6))  
  
  if (grepl("fate", data[i], fixed = TRUE)) {
    
    year.tmp <- length(Year) - 1
    
  } else {
    
    year.tmp <- length(Year)
    
  }
  
  for (j in 1:year.tmp) {
    
    load(file = paste0("~/Documentos/Datos NO publicados/BioIntForest/PPA_Results/Reports/PPA_", Year[j], ".RData"))
    
    data.ctrl <- get(paste0(data[i],".ctrl"))
    data.thnn <- get(paste0(data[i],".thnn"))
    data.gof.tmp <- get(paste0(data[i],".gof"))

    if (j == 1) {
      
      plot.ctrl <- data.ctrl
      plot.thnn <- data.thnn
      
    }
    
    data.cat <- rbind(data.cat, 
                      data.frame(Year = Year[j], Treat = "Ctrl", r = data.ctrl$r, value = ppp.cat(data.ctrl)),
                      data.frame(Year = Year[j], Treat = "Thnn", r = data.thnn$r, value = ppp.cat(data.thnn)) )
    
    data.gof <- rbind(data.gof,
                      data.frame(Year = Year[j], Treat = "Ctrl", data.gof.tmp[data.gof.tmp$Plot == "Ctrl",]),
                      data.frame(Year = Year[j], Treat = "Thnn", data.gof.tmp[data.gof.tmp$Plot == "Thnn",]) )
    
  }
  
  Treat <- unique(data.cat$Treat)
  Plot.list <- list()
  #https://plotnine.readthedocs.io/en/stable/generated/plotnine.geoms.geom_segment.html
  
  for (j in Treat) { 
    
    Plot.list[[j]] <- ggplot(aes(x = r, y = Year), data = data.cat[data.cat$Treat == j,]) + 
      geom_raster(aes(fill = value)) + 
      # geom_vline(xintercept = data.gof$r.max[data.gof$Plot == j]) +
      scale_fill_gradient2(low = "indianred3", mid = "white", high = "darkolivegreen3", midpoint = .0) +
      
      #para que salga proporcional, el atributo "size" hay que ir cambiandolo en las dos lineas, sino se queda uno mas pequeño que el otro.
      geom_segment(aes(x = r.min, y = Year, xend = r.max, yend = Year, size = 0.005, alpha = 1/10), 
                   data = data.gof[data.gof$Plot == j,]) +
      
      # para incluir los intervalos de confianza
      geom_label(aes(x = r.min, y = Year, label = ""), data = data.gof[data.gof$Plot == j,], alpha = 1/10) +
      geom_label(aes(x = r.max, y = Year, label = ""), data = data.gof[data.gof$Plot == j,], alpha = 1/10) +
      
      #para incluir el valor del test GoF
      geom_text(aes(label = paste("GoF", signif.num(p.value)), r.min + 0.0, size = 0.005, vjust = -0.8, hjust = "left"), data = data.gof[data.gof$Plot == j,]) +

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