library(sp)
library(raster)
library(spatstat.geom)

# Function to convert raster to im (sp) file ------------------------------

raster.as.im <- function(im) {
  
  r <- raster::res(im)
  orig <- sp::bbox(im)[, 1] + 0.5 * r
  dm <- dim(im)[2:1]
  xx <- unname(orig[1] + cumsum(c(0, rep(r[1], dm[1] - 1))))
  yy <- unname(orig[2] + cumsum(c(0, rep(r[2], dm[2] - 1))))
  
  return(spatstat.geom::im(matrix(raster::values(im), ncol = dm[1], 
                             nrow = dm[2], byrow = TRUE)[dm[2]:1, ], 
                      xcol = xx, yrow = yy))
}

#https://gis.stackexchange.com/questions/115159/converting-raster-to-im-object-for-point-process-model-covariate-in-r



#unction to normalize stacks

normalize_stack <- function(x) {
  
  range01.raster <- function(x) { (x - min(x, na.rm = T))/(max(x, na.rm = T) - min(x, na.rm = T)) }
  
  if (is.im(x[[1]]) == TRUE) {
    
    n.ls <- length(x)  
    
  } else {
    
    n.ls <- dim(x)[3]  
    
  }
  
  ls <- x
  
  for (i in 1:n.ls) {
    
    if (is.im(ls[[i]]) == TRUE) {
      
      ls[[i]]$v <- range01.raster(x[[i]]$v)
      
    } else {
      
      ls[[i]]$v <- range01.raster(ls[[i]]$v)
      
    }
    
    if (unique(is.na(ls[[i]]$v)) == TRUE) ls[[i]]$v <- 0
    
  }
  
  return(ls)
  
}
