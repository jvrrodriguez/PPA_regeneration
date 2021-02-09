##########################################
####   SUMMARIZE COUNTS             ######
##########################################


library(ggplot2)
library(dplyr)
library(ggpubr)

setwd("~/Documentos/Datos NO publicados/BioIntForest/Results")

Year <- c(2008, 2009, 2010, 2011)

data.sp <- data.frame(); data.size <- data.frame(); data.fate <-data.frame()

for (j in 1:length(Year)) {
  
  data.year <- read.csv(paste0("summary.Freq_", Year[j],".csv"))
  
  Plots <- unique(data.year$Plot)
  
  data.sp.year <- data.year[data.year$Var == "Sp",]
  data.size.year <- data.year[data.year$Var == "Size",]
  data.fate.year <- data.year[data.year$Var == "Fate",]
  
  prop.sp <- NULL; prop.size <- NULL; prop.fate <- NULL
  
  for (i in Plots) {
    
    data.sp.plot <- data.sp.year[data.sp.year$Plot == i,]
    prop.sp <- c(prop.sp, data.sp.plot$Freq / sum(data.sp.plot$Freq))
    
    data.size.plot <- data.size.year[data.size.year$Plot == i,]
    prop.size <- c(prop.size, data.size.plot$Freq / sum(data.size.plot$Freq))
    
    if (j != length(Year)) {
      
      data.fate.plot <- data.fate.year[data.fate.year$Plot == i,]
      prop.fate <- c(prop.fate, data.fate.plot$Freq / sum(data.fate.plot$Freq))  
      
    }

  } 

  data.sp.year <- data.frame(data.sp.year, Prop = prop.sp)
  data.size.year <- data.frame(data.size.year, Prop = prop.size)
  
  data.sp <- rbind(data.sp, data.sp.year)
  data.size <- rbind(data.size, data.size.year)
  
  if (j != length(Year)) {
    
    data.fate.year <- data.frame(data.fate.year, Prop = prop.fate)
    data.fate <- rbind(data.fate, data.fate.year)
  
  }
  
}

summary.sp <- data.sp %>% group_by(Year, Var, Var1, Treat) %>% summarise(m = mean(Prop), stdv = sd(Prop))
summary.size <- data.size %>% group_by(Year, Var, Var1, Treat) %>% summarise(m = mean(Prop), stdv = sd(Prop))
summary.fate <- data.fate %>% group_by(Year, Var, Var1, Treat) %>% summarise(m = mean(Prop), stdv = sd(Prop))


p.fs <- ggplot(summary.sp[summary.sp$Var1=="Fs",], aes(Year, m, fill = Treat)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  scale_y_continuous(limits=c(0,0.1)) + 
  geom_errorbar(aes(ymin=m-stdv, ymax=m+stdv), width=.2, position=position_dodge(.9)) +  theme_bw()

p.qh <- ggplot(summary.sp[summary.sp$Var1=="Qh",], aes(Year, m, fill = Treat)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  scale_y_continuous(limits=c(0,0.5)) + 
  geom_errorbar(aes(ymin=m-stdv, ymax=m+stdv), width=.2, position=position_dodge(.9)) +  theme_bw()

p.qi <- ggplot(summary.sp[summary.sp$Var1=="Qi",], aes(Year, m, fill = Treat)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  scale_y_continuous(limits=c(0,0.4)) + 
  geom_errorbar(aes(ymin=m-stdv, ymax=m+stdv), width=.2, position=position_dodge(.9)) +  theme_bw()

p.rec <- ggplot(summary.size[summary.size$Var1=="Recruit",], aes(Year, m, fill = Treat)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  scale_y_continuous(limits=c(0,0.5)) + 
  geom_errorbar(aes(ymin=m-stdv, ymax=m+stdv), width=.2, position=position_dodge(.9)) +  theme_bw()

p.fate <- ggplot(summary.fate[summary.fate$Var1=="1",], aes(Year, m, fill = Treat)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  scale_y_continuous(limits=c(0,0.5)) + 
  geom_errorbar(aes(ymin=m-stdv, ymax=m+stdv), width=.2, position=position_dodge(.9)) +  theme_bw()

ggarrange(p.fs, p.qh, p.qi, p.rec, p.fate + rremove("x.text"), 
          labels = c("a)", "b)", "c)", "d)", "e)"),
          ncol = 3, nrow = 2, common.legend = TRUE, legend = "right")
