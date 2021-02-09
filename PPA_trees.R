# Add required packages 
require(sp)
require(rgdal)
require(spatstat)
library(maptools)

setwd("~/Documentos/Datos NO publicados/BioIntForest/Data/Plots")

Studyarea <- readOGR("parcela 7 poli.shp")

data_1999 <- readOGR("Trees/Pinos_vivos_P7_1999.shp")
data_2009 <- readOGR("Trees/Pinos_vivos_P7.shp")
data_2020a <- readOGR("Trees/Arboles_P7_2020.shp")
data_2020 <- subset(data_2020a, Trat.Espe=="P")

Studyarea <- spTransform(Studyarea, "+proj=utm +zone=30 +ellps=GRS80 +units=m +no_defs")

# Coerce study area to owin object
w <- as.owin(Studyarea)

ppp_1999 <- with(data_1999, ppp(x = data_1999$X, y = data_1999$Y, owin(w)))
ppp_2009 <- with(data_2009, ppp(x = data_2009$X, y = data_2009$Y, owin(w)))
ppp_2020a <- with(data_2020a, ppp(x = data_2020a$X, y = data_2020a$Y, owin(w)))
ppp_2020 <- with(data_2020, ppp(x = data_2020$X, y = data_2020$Y, owin(w)))

nsim = 199

L.E99 <- envelope(ppp_1999, Lest, r = seq(0,10,0.05), nsim=nsim, fix.n=TRUE, correction="Ripley")
g.E99 <- envelope(ppp_1999, pcf, r = seq(0,4,0.02), nsim=nsim, fix.n=TRUE, correction="Ripley")
kNN.E99 <- envelope(ppp_1999, Gest, r = seq(0,1,0.01), nsim=nsim, fix.n=TRUE, correction="rs")
F.E99 <- envelope(ppp_1999, Fest, r = seq(0,1,0.01), nsim=nsim, fix.n=TRUE, correction="rs")

L.E09 <- envelope(ppp_2009, Lest, r = seq(0,10,0.05), nsim=nsim, fix.n=TRUE, correction="Ripley")
g.E09 <- envelope(ppp_2009, pcf, r = seq(0,4,0.02), nsim=nsim, fix.n=TRUE, correction="Ripley")
kNN.E09 <- envelope(ppp_2009, Gest, r = seq(0,1,0.01), nsim=nsim, fix.n=TRUE, correction="rs")
F.E09 <- envelope(ppp_2009, Fest, r = seq(0,1,0.01), nsim=nsim, fix.n=TRUE, correction="rs")

L.E20 <- envelope(ppp_2020, Lest, r = seq(0,10,0.05), nsim=nsim, fix.n=TRUE, correction="Ripley")
g.E20 <- envelope(ppp_2020, pcf, r = seq(0,4,0.02), nsim=nsim, fix.n=TRUE, correction="Ripley")
kNN.E20 <- envelope(ppp_2020, Gest, r = seq(0,1,0.01), nsim=nsim, fix.n=TRUE, correction="rs")
F.E20 <- envelope(ppp_2020, Fest, r = seq(0,1,0.01), nsim=nsim, fix.n=TRUE, correction="rs")

L.E20a <- envelope(ppp_2020a, Lest, r = seq(0,10,0.05), nsim=nsim, fix.n=TRUE, correction="Ripley")
g.E20a <- envelope(ppp_2020a, pcf, r = seq(0,4,0.02), nsim=nsim, fix.n=TRUE, correction="Ripley")
kNN.E20a <- envelope(ppp_2020a, Gest, r = seq(0,1,0.01), nsim=nsim, fix.n=TRUE, correction="rs")
F.E20a <- envelope(ppp_2020a, Fest, r = seq(0,1,0.01), nsim=nsim, fix.n=TRUE, correction="rs")


par(mfrow=c(4,3), mar=c(1, 1, 1.25, 1.25), oma = c(4, 4, 2, 2)) 
plot(L.E99, . -r ~ r, shade=c("hi", "lo"), legend = F, main = NULL)
plot(g.E99, main = NULL, legend = F)
plot(kNN.E99, main = NULL, legend = F)
#plot(F.E99, main = NULL, legend = F)

plot(L.E09, . -r ~ r, shade=c("hi", "lo"), legend = F, main = NULL)
plot(g.E09, main = NULL, legend = F)
plot(kNN.E09, main = NULL, legend = F)
#plot(F.E09, main = NULL, legend = F)

plot(L.E20, . -r ~ r, shade=c("hi", "lo"), legend = F, main = NULL)
plot(g.E20, main = NULL, legend = F)
plot(kNN.E20, main = NULL, legend = F)
#plot(F.E20, main = NULL, legend = F)

plot(L.E20a, . -r ~ r, shade=c("hi", "lo"), legend = F, main = NULL)
plot(g.E20a, main = NULL, legend = F)
plot(kNN.E20a, main = NULL, legend = F)
#plot(F.E20a, main = NULL, legend = F)