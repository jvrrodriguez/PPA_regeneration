##########################################
####   SUMMARIZE COUNTS             ######
##########################################

library(ggplot2)
library(dplyr)
library(ggpubr)
library(vegan)

setwd("~/Documentos/Datos NO publicados/BioIntForest/Results")

Year <- c(2008, 2009, 2010, 2011)

Plots <- c(2,3,4,5,7,9)

Treat <-  c("30%", "0%", "0%", "30%", "30%", "0%") #ONLY 0% and 30%

data.sp <- NULL; data.size <- NULL
data.rec.sp <- NULL; data.rec.size <- NULL; data.rec.fate <- NULL
data.rec.rich <- NULL; data.rec.shan <- NULL

save.output <- TRUE

if (save.output == T) pdf("~/Documentos/Datos NO publicados/BioIntForest/Results/Count_data.pdf", width=10, height=5)  


for (j in 1:length(Year)) {
  
  plot.year <- read.csv(paste0("summary.Plot_", Year[j],".csv"))
  qu.year <- read.csv(paste0("summary.Quadrat_", Year[j],".csv"))
  qu.year <- cbind(qu.year, sp.shan = apply(qu.year[,c(5,6,7)], 1, diversity),
                   sp.rich = rowSums(ifelse(as.matrix(qu.year[,c(5,6,7)]), 1, 0)),
                   Dens = rowSums(as.matrix(qu.year[,c(5,6,7)])))
  
  plot.year.sp <- plot.year %>% group_by(Year, Plot, Treat, sp, size) %>% summarise(n = sum(n))
  
  plot.year.size <- plot.year %>% group_by(Year, Plot, Treat, size) %>% summarise(n = sum(n))
  
  prop.rec.sp <- NULL; prop.rec.size <- NULL; prop.rec.fate <- NULL
  prop.rec.rich <- NULL; prop.rec.shan <- NULL
  
  for (i in Plots) {
    
    prop.rec.rich <-  c(prop.rec.rich, mean(qu.year[qu.year$Plot == i,]$sp.rich[qu.year[qu.year$Plot == i,]$sp.rich > 0], na.rm=T))
    prop.rec.shan <- c(prop.rec.shan, mean(qu.year[qu.year$Plot == i,]$sp.shan[qu.year[qu.year$Plot == i,]$sp.shan > 0], na.rm=T))
    
    plot.rec.sp <- plot.year.sp[plot.year.sp$Plot == i & plot.year.sp$size == "Recruit",]
    prop.rec.sp <- c(prop.rec.sp, plot.rec.sp$n / sum(plot.rec.sp$n))
    
    prop.rec.size <- c(prop.rec.size, mean(qu.year[qu.year$Plot == i,]$Dens, na.rm=T))
    #plot.size <- plot.year.size[plot.year.size$Plot == i,]
    #prop.size <- c(prop.size, plot.size$Freq / sum(plot.size$Freq))
    
    if (j != length(Year)) {
      
      plot.year.fate <- plot.year %>% group_by(Year, Plot, Treat, fate, size) %>% summarise(n = sum(n))
      plot.rec.fate <- plot.year.fate[plot.year.fate$Plot == i & plot.year.fate$size == "Recruit",]
      prop.rec.fate<- c(prop.rec.fate, plot.rec.fate$n / sum(plot.rec.fate$n))  
      
    }
    
  } 
  
  plot.year.rec.sp <- data.frame(plot.year.sp[plot.year.sp$size == "Recruit",], Prop = prop.rec.sp)
  plot.year.rec.rich <- data.frame(Year = Year[j], Plot = Plots, Treat = Treat, sp.rich = "Recruit", Prop = prop.rec.rich)
  plot.year.rec.shan <- data.frame(Year = Year[j], Plot = Plots, Treat = Treat, sp.rich = "Recruit", Prop = prop.rec.shan)
  
  plot.year.rec.size <- data.frame(plot.year.size[plot.year.size$size == "Recruit" & !is.na(plot.year.size$size),], Prop = prop.rec.size)
  
  data.sp <- rbind(data.sp, plot.year.sp)

  data.rec.sp <- rbind(data.rec.sp, plot.year.rec.sp)
  data.rec.rich <- rbind(data.rec.rich, plot.year.rec.rich)
  data.rec.shan <- rbind(data.rec.shan, plot.year.rec.shan)
  
  data.size <- rbind(data.size, plot.year.size)
  data.rec.size <- rbind(data.rec.size, plot.year.rec.size)
  
  if (j != length(Year)) {
    
    plot.year.rec.fate <- data.frame(plot.year.fate[plot.year.fate$size == "Recruit",], Prop = prop.rec.fate)
    data.rec.fate <- rbind(data.rec.fate, plot.year.rec.fate)
    
  }
  
}

summary.size <- data.size %>% group_by(Year, size, Treat) %>% summarise(sum = sum(n, na.rm=T), stdv = sd(n, na.rm=T), se = sd(n, na.rm=T)/sqrt(n()))
summary.size$size <- ifelse (summary.size$size == "Recruits", "Recruit", as.character(summary.size$size))

data.sp[data.sp$size == "Recruit",] %>% group_by(Treat, sp) %>% summarise(n = sum(sum, na.rm=T), stdv = sd(sum, na.rm=T), se = sd(sum, na.rm=T)/sqrt(n()))
summary.size %>% group_by(size) %>% summarise(mean = mean(sum, na.rm=T), stdv = sd(sum, na.rm=T), se = sd(sum, na.rm=T)/sqrt(n()))

summary.sp <- data.sp %>% group_by(Year, size, Treat, sp) %>% summarise(sum = sum(n, na.rm=T), stdv = sd(n, na.rm=T), se = sd(n, na.rm=T)/sqrt(n()))
summary.sp %>% group_by(size, sp) %>% summarise(mean = mean(sum, na.rm=T), stdv = sd(sum, na.rm=T), se = sd(sum, na.rm=T)/sqrt(n()))

summary.rec.sp <- data.rec.sp %>% group_by(Year, Treat, sp) %>% summarise(m = mean(n, na.rm=T), stdv = sd(n, na.rm=T), se = sd(n, na.rm=T)/sqrt(n()))
summary.rec.rich <- data.rec.rich %>% group_by(Year, Treat) %>% summarise(m = mean(Prop, na.rm=T), stdv = sd(Prop, na.rm=T), se = sd(Prop, na.rm=T)/sqrt(n()))
summary.rec.shan <- data.rec.shan %>% group_by(Year, Treat) %>% summarise(m = mean(Prop, na.rm=T), stdv = sd(Prop, na.rm=T), se = sd(Prop, na.rm=T)/sqrt(n()))
summary.rec.size <- data.rec.size %>% group_by(Year, Treat) %>% summarise(m = mean(Prop, na.rm=T), stdv = sd(Prop, na.rm=T), se = sd(Prop, na.rm=T)/sqrt(n()))
summary.rec.size[,3:5] <- summary.rec.size[,3:5] / 4 #para pasarlo a metro cuadrado

summary.rec.fate <- data.rec.fate %>% group_by(Year, Treat, fate) %>% summarise(m = mean(Prop, na.rm=T), stdv = sd(Prop, na.rm=T), se = sd(Prop, na.rm=T)/sqrt(n()))

summary.rec.sp %>% group_by(Treat) %>% summarise(sum = mean(m, na.rm=T), stdv = sd(m, na.rm=T), se = sd(m, na.rm=T)/sqrt(n()))
summary.rec.sp[,4:6] <- summary.rec.sp[,4:6] / (30*40) #para pasarlo a metro cuadrado

summary.rec.rich %>% group_by(Treat) %>% summarise(sum = mean(m, na.rm=T), stdv = sd(m, na.rm=T), se = sd(m, na.rm=T)/sqrt(n()))
summary.rec.rich[,3:5] <- summary.rec.rich[,3:5] / 4 #para pasarlo a metro cuadrado

summary.rec.shan %>% group_by(Treat) %>% summarise(sum = mean(m, na.rm=T), stdv = sd(m, na.rm=T), se = sd(m, na.rm=T)/sqrt(n()))


p.rich <- ggplot(summary.rec.rich, aes(Year, m, fill = Treat)) + labs(y = expression(paste("Richness of recruits (", m^2, ")"))) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  #scale_y_continuous(limits=c(0,0.1)) + 
  geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.2, position=position_dodge(.9)) +  theme_classic()

p.shan <- ggplot(summary.rec.shan, aes(Year, m, fill = Treat)) + labs(y = "Shannon by quadrat") + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  #scale_y_continuous(limits=c(0.6,0.8)) + 
  geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.2, position=position_dodge(.9)) +  theme_classic()


p.rec <- ggplot(summary.rec.size, aes(Year, m, fill = Treat)) + labs(y = expression(paste("Density of recruits (", m^2, ")"))) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  #scale_y_continuous(limits=c(0,0.5)) + 
  geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.2, position=position_dodge(.9)) +  theme_classic()

p.fate <- ggplot(summary.rec.fate[summary.rec.fate$fate=="1",], aes(Year, m, fill = Treat)) + labs(y = "Proportion of recruit's survival") +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  #scale_y_continuous(limits=c(0,0.5)) + 
  geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.2, position=position_dodge(.9)) +  theme_classic()

ggarrange(p.rec, p.rich, p.fate, #+ rremove("x.text")
          labels = c("a)", "b)", "c)"), #"d)", "e)", "f)", "g)"
          ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")



p.fs <- ggplot(summary.rec.sp[summary.rec.sp$sp=="Fs",], aes(Year, m, fill = Treat)) + labs(y = expression(paste("Density of F. sylvatica recruits (", m^2, ")"))) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  #scale_y_continuous(limits=c(0,0.1)) +
  geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.2, position=position_dodge(.9)) +  theme_classic()

p.qh <- ggplot(summary.rec.sp[summary.rec.sp$sp=="Qh",], aes(Year, m, fill = Treat)) + labs(y = expression(paste("Density of Q. humilis recruits (", m^2, ")"))) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  #scale_y_continuous(limits=c(0,1)) +
  geom_errorbar(aes(ymin=m-se, ymax=m+se) , width=.2, position=position_dodge(.9)) +  theme_classic()

p.qi <- ggplot(summary.rec.sp[summary.rec.sp$sp=="Qi",], aes(Year, m, fill = Treat)) + labs(y = expression(paste("Density of  Q. ilex recruits (", m^2, ")"))) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  #scale_y_continuous(limits=c(0,0.1)) +
  geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.2, position=position_dodge(.9)) +  theme_classic()

ggarrange(p.qh, p.qi, p.fs, # + rremove("x.text") 
          labels = c("a)", "b)", "c)"), #"d)", "e)", "f)", "g)"
          ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")

if (save.output == T) dev.off()