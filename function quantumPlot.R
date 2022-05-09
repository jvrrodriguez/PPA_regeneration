#https://www.r-bloggers.com/quantumplots-with-ggplot2-and-spatstat/

quantumPlot <- function(x, quantum = T, log.y = T, y_intercept = 1, y_name = "summary statistic", colour=c("#d73027", "#ffffbf", "#91bfdb")){
  
  # load Packages
  require(ggplot2)
  require(ggthemes)
  
  # convert fv to dataframe
  env.data <- as.data.frame(x)
  env.data <- env.data[-1,]

  if (is.null(env.data$pooliso) == F) {
    
    env.data$obs <- sqrt(env.data$pooliso/pi) - env.data$r
    env.data$hi <- sqrt(env.data$poolhi/pi) - env.data$r
    env.data$lo <- sqrt(env.data$poollo/pi) - env.data$r

  }
  
  # plot it
  
  if (log.y == F) { 
    
  gg_quantomPlot <- ggplot(env.data, aes(r, obs)) +
    # plot observed value
    geom_line(colour=c("#4d4d4d")) +
    # plot simulation envelopes
    geom_ribbon(aes(ymin=lo,ymax=hi),alpha=0.1, colour=c("#e0e0e0")) +
    # axes names and limits
    #ylim(min(env.data$obs)-0.5, max(env.data$obs)+2) +
    labs(x="Distance r (m)", y=y_name) +
    # plot expected value, according to null model
    geom_hline(yintercept = y_intercept, linetype = "dashed", colour = c("#999999")) + 
    # make it look beautiful
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  }
    
   if (log.y == T) { 
     
     gg_quantomPlot <- ggplot(env.data, aes(r, obs)) +
       # plot observed value
       geom_line(colour=c("#4d4d4d")) + scale_y_log10() +
       #ylim(min(env.data$obs, na.rm=T) - 0.5, max(env.data$hi, na.rm=T) - max(env.data$hi, na.rm=T) * 0.1) +
       # plot simulation envelopes
       geom_ribbon(aes(ymin=lo,ymax=hi),alpha=0.1, colour=c("#e0e0e0")) +
       # axes names and limits
       #ylim(min(env.data$obs)-0.5, max(env.data$obs)+2) +
       labs(x="Distance r (m)", y=y_name) +
       # plot expected value, according to null model
       geom_hline(yintercept = y_intercept, linetype = "dashed", colour = c("#999999")) +
       # make it look beautiful
       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black"))

     }

   if (quantum == T) {
     
     # plot 'Quantums'
      gg_quantomPlot + 
       geom_rug(data=env.data[env.data$obs > env.data$hi,], sides="b", colour=colour[1]) +
       geom_rug(data=env.data[env.data$obs < env.data$lo,], sides="b", colour=colour[2]) +
       geom_rug(data=env.data[env.data$obs >= env.data$lo & env.data$obs <= env.data$hi,], sides="b", color=colour[3])
     
   }
    
  return(gg_quantomPlot)
}


ppp.cat <- function(data) {
  
  if (is.null(data$obs) == F) output <- ifelse(data$obs < data$lo, -1, ifelse(data$obs > data$hi, 1, 0))
  if (is.null(data$pooliso) == F) output <- ifelse( (sqrt(data$pooliso/pi) - data$r) < (sqrt(data$poollo/pi) - data$r), -1, 
                                                    ifelse( (sqrt(data$pooliso/pi) - data$r) > (sqrt(data$poolhi/pi) - data$r), 1, 0))
  if (is.null(data$iso) == F) output <- ifelse( (sqrt(data$iso/pi) - data$r) < (sqrt(data$lo/pi) - data$r), -1, 
                                                    ifelse((sqrt(data$iso/pi) - data$r) > (sqrt(data$hi/pi) - data$r), 1, 0))
  
  return(output)
  
}