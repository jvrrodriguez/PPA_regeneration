### ADDITIONAL FUNCTIONS FOR MAPPING, ANALYSIS AND PLOTTING


# functions for mapping ------------------------------------------------------

range01 <- function(x){
  
  x <- ifelse(x == -9999, NA, x) 
  x0 <- (x - min(x, na.rm = T))/(max(x, na.rm = T) - min(x, na.rm = T))
  return(x0)
  
}


im.diversity <- function(x, y, z){
  
  require(vegan)
  
  results <- matrix(NA,nrow(x),nrow(y)); results2 <- results
  thr <- c(median(x), median(y), median(z))
  
  for (i in 1:nrow(x)) {
    
    for (j in 1:nrow(y)) {
      
      r1 <- x[i,j]
      r2 <- y[i,j]
      r3 <- z[i,j]
      results[i,j] <- sum(c(ifelse(r1 >= thr[1], 1, 0), ifelse(r2 >= thr[2], 1, 0), ifelse(r3 >= thr[3], 1, 0)))  ## Example function
      results2[i,j] <- diversity(c(abs(r1) , abs(r2) , abs(r3)), index = "shan")  ## Example function
      
    }
    
  }
  
  return(list(results, results2))
  
}



# Functions for Point-pattern analyses ------------------------------------

# http://rosanaferrero.blogspot.com/2010/12/procesos-espaciales-puntuales-con-r.html
# https://www.rdocumentation.org/packages/ads/versions/1.5-3/topics/k12fun
# Change markconnect by Jdif. Why?

Jdif <- function(X, ..., i) {
  
  Jidot <- Jdot(X, ..., i = i) #se selecciona generalmente la correccion de Ripley
  J <- Jest(X, ...)
  dif <- eval.fv(Jidot - J)
  return(dif)
  
}



# Functions for plotting --------------------------------------------------

#https://www.r-bloggers.com/quantumplots-with-ggplot2-and-spatstat/

quantumPlot <- function(x, quantum = T, log.y = T, y_intercept = 1, y_name = "summary statistic", colour=c("#d73027", "#ffffbf", "#91bfdb")){
  
  # load Packages
  require(ggplot2)
  require(ggthemes)
  
  # convert fv to dataframe
  data <- as.data.frame(x)
  data <- data[-1,]

  if (is.null(data$pooliso) == F) {
    
    data$obs <- sqrt(data$pooliso/pi) - data$r
    data$hi <- sqrt(data$poolhi/pi) - data$r
    data$lo <- sqrt(data$poollo/pi) - data$r

  }
  
  if (attributes(x)$fname[1] == "K") {
    
    data$obs <- sqrt(data$obs/pi) - data$r
    data$hi <- sqrt(data$hi/pi) - data$r
    data$lo <- sqrt(data$lo/pi) - data$r
    
  }
  
  
  # plot it
  
  if (log.y == F) { 
    
  gg_quantomPlot <- ggplot(data, aes(r, obs)) +
    # plot observed value
    geom_line(colour = c("#4d4d4d")) +
    # plot simulation envelopes
    geom_ribbon(aes(ymin = lo,ymax = hi),alpha = 0.1, colour = c("#e0e0e0")) +
    # axes names and limits
    #ylim(min(data$obs)-0.5, max(data$obs)+2) +
    labs(x = "Distance r (m)", y = y_name) +
    # plot expected value, according to null model
    geom_hline(yintercept = y_intercept, linetype = "dashed", colour = c("#999999")) + 
    # make it look beautiful
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  }
    
   if (log.y == T) { 
     
     gg_quantomPlot <- ggplot(data, aes(r, obs)) +
       # plot observed value
       geom_line(colour = c("#4d4d4d")) + scale_y_log10() +
       #ylim(min(data$obs, na.rm=T) - 0.5, max(data$hi, na.rm=T) - max(data$hi, na.rm=T) * 0.1) +
       # plot simulation envelopes
       geom_ribbon(aes(ymin = lo,ymax = hi), alpha = 0.1, colour = c("#e0e0e0")) +
       # axes names and limits
       #ylim(min(data$obs)-0.5, max(data$obs)+2) +
       labs(x = "Distance r (m)", y = y_name) +
       # plot expected value, according to null model
       geom_hline(yintercept = y_intercept, linetype = "dashed", colour = c("#999999")) +
       # make it look beautiful
       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black"))

     }

   if (quantum == T) {
     
     # plot 'Quantums'
      gg_quantomPlot + 
       geom_rug(data = data[data$obs > data$hi,], sides = "b", colour = colour[1]) +
       geom_rug(data = data[data$obs < data$lo,], sides = "b", colour = colour[2]) +
       geom_rug(data = data[data$obs >= data$lo & data$obs <= data$hi,], sides = "b", color = colour[3])
     
   }
    
  return(gg_quantomPlot)
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


# para categorizar cuando el valor observado sale fuera del intervalo de confianza

ppp.cat <- function(data) {
  
  if (is.null(data$obs) == F) output <- ifelse(data$obs < data$lo, -1, ifelse(data$obs > data$hi, 1, 0))
  if (is.null(data$pooliso) == F) output <- ifelse( (sqrt(data$pooliso/pi) - data$r) < (sqrt(data$poollo/pi) - data$r), -1, 
                                                    ifelse( (sqrt(data$pooliso/pi) - data$r) > (sqrt(data$poolhi/pi) - data$r), 1, 0))
  if (is.null(data$iso) == F) output <- ifelse( (sqrt(data$iso/pi) - data$r) < (sqrt(data$lo/pi) - data$r), -1, 
                                                    ifelse((sqrt(data$iso/pi) - data$r) > (sqrt(data$hi/pi) - data$r), 1, 0))
  
  return(output)
  
}

 # ppp.cat <- function(data) {
 #   
 #   if (is.null(data$obs) == F) output <- ifelse(data$obs < data$lo, -1, ifelse(data$obs > data$hi, 1, 0))
 #   if (is.null(data$iso) == F) output <- ifelse(data$iso < data$lo, -1, ifelse(data$iso > data$hi, 1, 0))
 #   
 #   return(output)
 #   
 # }


# para obtener un rango de valores donde sale fuera el intervalo de confianza y calcular el GoF

gof.int <- function(data, r.int) {
  
  require(spatstat.core)
  
  data.cat <- abs(ppp.cat(data))
  
  cons.values <- rollsum(data.cat, r.int, fill = NA, align = "center") == r.int
  data.interval <- c(min(data$r[cons.values], na.rm = T), max(data$r[cons.values], na.rm = T))
  
  if (is.infinite(data.interval[1]) | (data.interval[1]) == (data.interval[2])) {
    
    #si es infinito, o son iguales los dos valores, se le asigna valores de cero a todo
    
    data.sum <- data.frame(r.min = 0, r.max = 0, 
                           e.out = sum(data.cat) / length(data.cat), 
                           u = 0, p.value = 1)
    
  } else {
    
    data.test <- dclf.test(data, rinterval = data.interval)
    
    data.sum <- data.frame(r.min = attributes(data.test)$rinterval[1], r.max = attributes(data.test)$rinterval[2], 
                           e.out = sum(data.cat) / length(data.cat), 
                           u = data.test$statistic[1], p.value = data.test$p.value)
    
  }

  return(data.sum)
  
}


# categorizar en simbolos los valores de significación

signif.num <- function(x) {
  symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
         cutpoints = c(0, 0.001, 0.01, 0.05, 1), #0.1, para rehacerlo, ponerlo el segundo por la cola 
         symbols = c("***", "**", "*", "n.s")) #".", 
}