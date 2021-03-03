##########################################
####   SUMMARIZE COUNTS             ######
##########################################

library(ggplot2)
library(dplyr)
library(ggpubr)
library(vegan)

Year <- 2007+1:4
Plots <- c(2,3,4,5,7,9)
Treat <- c("20%", "30%", "0%", "0%", "30%", "20%", "30%", "20%", "0%")[Plots]

data.sp <- NULL; data.size <- NULL
data.rec.sp <- NULL; data.rec.size <- NULL; data.rec.fate <- NULL
data.rec.rich <- NULL; data.rec.shan <- NULL
summary.rec.Q.yr <- NULL; summary.ad.Q.yr <- NULL

save.output <- FALSE

if (save.output == T) pdf("~/Documentos/Datos NO publicados/BioIntForest/PPA_Results/Figures/Count_data.pdf", width=10, height=4)  

for (j in 1:length(Year)) {
  
  setwd("~/Documentos/Datos NO publicados/BioIntForest/PPA_Results/Reports")
  summary.plot <- read.csv(paste0("summary.Plot_", Year[j],".csv"))
  summary.rec.Q <- read.csv(paste0("summary.recruit.Q_", Year[j],".csv"))
  summary.rec.Q <- cbind(summary.rec.Q, sp.shan = apply(summary.rec.Q[,c(5,6,7)], 1, diversity),
                   sp.rich = rowSums(ifelse(as.matrix(summary.rec.Q[,c(5,6,7)]), 1, 0)),
                   Dens = rowSums(as.matrix(summary.rec.Q[,c(5,6,7)])))
  summary.rec.Q.yr <- rbind(summary.rec.Q.yr, summary.rec.Q)
  
  summary.ad.Q <- read.csv(paste0("summary.ad.Q_", Year[j],".csv"))
  summary.ad.Q.yr <- rbind(summary.ad.Q.yr, summary.ad.Q)
  
  summary.plot.sp <- summary.plot %>% group_by(Year, Plot, Treat, sp, size) %>% summarise(n = sum(n))
  
  summary.plot.size <- summary.plot %>% group_by(Year, Plot, Treat, size) %>% summarise(n = sum(n))
  
  prop.rec.sp <- NULL; prop.rec.size <- NULL; prop.rec.fate <- NULL
  prop.rec.rich <- NULL; prop.rec.shan <- NULL
  
  for (i in Plots) {
    
    prop.rec.rich <-  c(prop.rec.rich, mean(summary.rec.Q[summary.rec.Q$Plot == i,]$sp.rich[summary.rec.Q[summary.rec.Q$Plot == i,]$sp.rich > 0], na.rm=T))
    prop.rec.shan <- c(prop.rec.shan, mean(summary.rec.Q[summary.rec.Q$Plot == i,]$sp.shan[summary.rec.Q[summary.rec.Q$Plot == i,]$sp.shan > 0], na.rm=T))
    
    plot.rec.sp <- summary.plot.sp[summary.plot.sp$Plot == i & summary.plot.sp$size == "Recruit",]
    prop.rec.sp <- c(prop.rec.sp, plot.rec.sp$n / sum(plot.rec.sp$n))
    
    prop.rec.size <- c(prop.rec.size, mean(summary.rec.Q[summary.rec.Q$Plot == i,]$Dens, na.rm=T))
    #plot.size <- summary.plot.size[summary.plot.size$Plot == i,]
    #prop.size <- c(prop.size, plot.size$Freq / sum(plot.size$Freq))
    
    if (j != length(Year)) {
      
      summary.plot.fate <- summary.plot %>% group_by(Year, Plot, Treat, fate, size) %>% summarise(n = sum(n))
      plot.rec.fate <- summary.plot.fate[summary.plot.fate$Plot == i & summary.plot.fate$size == "Recruit",]
      prop.rec.fate<- c(prop.rec.fate, plot.rec.fate$n / sum(plot.rec.fate$n))  
      
    }
    
  } 
  
  summary.plot.rec.sp <- data.frame(summary.plot.sp[summary.plot.sp$size == "Recruit",], Prop = prop.rec.sp)
  summary.plot.rec.rich <- data.frame(Year = Year[j], Plot = Plots, Treat = Treat, sp.rich = "Recruit", Prop = prop.rec.rich)
  summary.plot.rec.shan <- data.frame(Year = Year[j], Plot = Plots, Treat = Treat, sp.rich = "Recruit", Prop = prop.rec.shan)
  
  summary.plot.rec.size <- data.frame(summary.plot.size[summary.plot.size$size == "Recruit" & !is.na(summary.plot.size$size),], Prop = prop.rec.size)
  
  data.sp <- rbind(data.sp, summary.plot.sp)

  data.rec.sp <- rbind(data.rec.sp, summary.plot.rec.sp)
  data.rec.rich <- rbind(data.rec.rich, summary.plot.rec.rich)
  data.rec.shan <- rbind(data.rec.shan, summary.plot.rec.shan)
  
  data.size <- rbind(data.size, summary.plot.size)
  data.rec.size <- rbind(data.rec.size, summary.plot.rec.size)
  
  if (j != length(Year)) {
    
    summary.plot.rec.fate <- data.frame(summary.plot.fate[summary.plot.fate$size == "Recruit",], Prop = prop.rec.fate)
    data.rec.fate <- rbind(data.rec.fate, summary.plot.rec.fate)
    
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
summary.rec.size <- data.rec.size %>% group_by(Year, Treat) %>% summarise(m = mean(Prop, na.rm=T), stdv = sd(n, na.rm=T), se = sd(Prop, na.rm=T)/sqrt(n()))
summary.ad.size <- data.size[data.size$size == "Adult",] %>% group_by(Year, Treat) %>% summarise(m = mean(n, na.rm=T), stdv = sd(n, na.rm=T), se = sd(n, na.rm=T)/sqrt(n()))
summary.rec.size[,3:5] <- summary.rec.size[,3:5] / 4 #para pasarlo a metro cuadrado
summary.ad.size[,3:5] <- (summary.ad.size[,3:5] / (30*40))


summary.rec.fate <- data.rec.fate %>% group_by(Year, Treat, fate) %>% summarise(m = mean(Prop, na.rm=T), stdv = sd(Prop, na.rm=T), se = sd(Prop, na.rm=T)/sqrt(n()))

summary.rec.sp %>% group_by(Treat) %>% summarise(sum = mean(m, na.rm=T), stdv = sd(m, na.rm=T), se = sd(m, na.rm=T)/sqrt(n()))
summary.rec.sp[,4:6] <- summary.rec.sp[,4:6] / (30*40) #para pasarlo a metro cuadrado

summary.rec.rich %>% group_by(Treat) %>% summarise(sum = mean(m, na.rm=T), stdv = sd(m, na.rm=T), se = sd(m, na.rm=T)/sqrt(n()))
summary.rec.rich[,3:5] <- summary.rec.rich[,3:5] / 4 #para pasarlo a metro cuadrado

summary.rec.shan %>% group_by(Treat) %>% summarise(sum = mean(m, na.rm=T), stdv = sd(m, na.rm=T), se = sd(m, na.rm=T)/sqrt(n()))


signif.num <- function(x) {
  symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
         cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
         symbols = c("***", "**", "*", ".", "n.s"))
}


summary.ad.Q.yr$Year <- as.factor(summary.ad.Q.yr$Year)  
mod.ad.size <- aov(Size ~ Treat:Year, data = summary.ad.Q.yr)
tukey.ad <- data.frame(TukeyHSD(mod.ad.size)$`Treat:Year`, label = array(signif.num(TukeyHSD(mod.ad.size)$`Treat:Year`[,4])))
tukey.ad <- tukey.ad[c(1,14,23,28),]
max.ad <- max(summary.ad.size$m + summary.ad.size$se)


p.ad <- ggplot(summary.ad.size, aes(Year, m, fill = Treat)) + labs(y = expression(paste("Density of trees (", m^2, ")"))) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  #scale_y_continuous(limits=c(0,0.5)) + 
  geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.2, position=position_dodge(.9)) +  
  geom_text(x = c(2008, 2008, 2009, 2009, 2010, 2010, 2011, 2011), y = max.ad, label = c(tukey.ad$label[1], "", tukey.ad$label[2], "", tukey.ad$label[3], "", tukey.ad$label[4], "")) +
  theme_classic()


data.rec.size <- data.rec.size %>% mutate(Treat = ifelse( Treat == "30%", "Thnn", "Ctrl" ))
summary.rec.Q.yr$Year <- as.factor(summary.rec.Q.yr$Year)  
mod.rec.size <- aov(Size ~ Treat:Year, data = summary.rec.Q.yr)
#mod.rec.size <- lme(Size ~ Treat+Year +Treat:Year, random=~1|Plot, data=summary.rec.Q.yr)
#emmeans(mod.rec.size, pairwise ~ Treat:Year)
tukey.rec <- data.frame(TukeyHSD(mod.rec.size)$`Treat:Year`, label = array(signif.num(TukeyHSD(mod.rec.size)$`Treat:Year`[,4])))
tukey.rec <- tukey.rec[c(1,14,23,28),]
max.rec <- max(summary.rec.size$m + summary.rec.size$se)

p.rec <- ggplot(summary.rec.size, aes(Year, m, fill = Treat)) + labs(y = expression(paste("Density of recruits (", m^2, ")"))) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  #scale_y_continuous(limits=c(0,0.5)) + 
  geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.2, position=position_dodge(.9)) +  
  #geom_signif(annotations = c("",""), y_position = max.rec-0.01, xmin=c(1,1), xmax=c(2,2)) +
  geom_text(x = c(2008, 2008, 2009, 2009, 2010, 2010, 2011, 2011), y = max.rec, label = c(tukey.rec$label[1], "", tukey.rec$label[2], "", tukey.rec$label[3], "", tukey.rec$label[4], "")) +
  theme_classic()


data.rec.fate <- data.rec.fate %>% mutate(Treat = ifelse( Treat == "30%", "Thnn", "Ctrl" ))

summary.rec.Q.yr$Year <- as.factor(summary.rec.Q.yr$Year)  
summary.rec.Q.yr$Fate <- ifelse(summary.rec.Q.yr$Fate > 0, 1, 0) 
mod.rec.fate <- aov(Fate ~ Treat:Year, data = summary.rec.Q.yr, family = "binomial")
mod.rec.fate <- lme(Fate ~ Treat:Year, random=~1|Plot, data=summary.rec.Q.yr)
tukey.fate <- data.frame(TukeyHSD(mod.rec.fate)$`Treat:Year`, label = array(signif.num(TukeyHSD(mod.rec.fate)$`Treat:Year`[,4])))
tukey.fate <- tukey.fate[c(1,14,23,28),]
max.fate <- max(summary.rec.fate[summary.rec.fate$fate=="1",]$m + summary.rec.fate[summary.rec.fate$fate=="1",]$se)

p.fate <- ggplot(summary.rec.fate[summary.rec.fate$fate=="1",], aes(Year, m, fill = Treat)) + labs(y = "Proportion of recruit's survival") +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  #scale_y_continuous(limits=c(0,0.5)) + 
  geom_text(x = c(2008, 2008, 2009, 2009, 2010, 2010), y = max.fate, label = c(tukey.fate$label[1], "", tukey.fate$label[2], "", tukey.fate$label[3], "")) +
  geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.2, position=position_dodge(.9)) +  theme_classic()


p.rich <- ggplot(summary.rec.rich, aes(Year, m, fill = Treat)) + labs(y = expression(paste("Richness of recruits (", m^2, ")"))) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  #scale_y_continuous(limits=c(0,0.1)) + 
  geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.2, position=position_dodge(.9)) +  theme_classic()

p.shan <- ggplot(summary.rec.shan, aes(Year, m, fill = Treat)) + labs(y = "Shannon by quadrat") + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  #scale_y_continuous(limits=c(0.6,0.8)) + 
  geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.2, position=position_dodge(.9)) +  theme_classic()


ggarrange(p.ad, p.rec, p.fate, #+ rremove("x.text")
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