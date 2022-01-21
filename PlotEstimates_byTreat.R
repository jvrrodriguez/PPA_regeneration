library(ggpubr)
# Data preparation
housetasks <- read.delim(
  system.file("demo-data/housetasks.txt", package = "ggpubr"),
  row.names = 1
)
head(housetasks, 4)

# Visualization
ggballoonplot(housetasks, fill = "value")+
  scale_fill_viridis_c(option = "C")


library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)

source("~/Documentos/Datos NO publicados/BioIntForest/PPA_regeneration/function quantumPlot.R")

Year <- c(2008, 2009, 2010, 2011)
Plot <- c(2,3,4,5,7,9)
save.output <- TRUE
Var.dens <- c("(Intercept)", "Treat30%",  "Dens.adult", "Dens.rec", "Dens.size.rec", "Dens.rich.rec", "Dens.shan.rec")
Var.env <- c("(Intercept)", "Treat30%", "Canopy", "CanOpen", "LAI", "DiffBelow", "DiffBelow.Yr", "N.Sunflecks", "Mdn.Sunflecks", "Max.Sunflecks", "Fs", "Hed", "Pter", "Rub", "Scl") #"MDT",  "MDS",  "Lad.var", "Lad.cv", "Lad.shan",

save.output <- TRUE

#Encontrar una interaccion entre especies y ver si esto varia con el tiempo (variables ambientales
#alfun factor que tenga una tendencia temporal

estimate.env <- NULL
estimate.dens <- NULL
p.env <- NULL
p.dens <- NULL

data.pred <- data.frame()

list.env.sp <- list()
list.dens.sp <- list()

list.env.fate <- list()
list.dens.fate <- list()

estimate.env.plot <- NULL
estimate.dens.plot <- NULL

p.env.plot <- NULL
p.dens.plot <- NULL

for (i in 1:length(Year)){
  
  load(file=paste0("~/Documentos/Datos NO publicados/BioIntForest/PPA_Results/Reports/PPA_", Year[i],".RData"))
  
  data.pred <- rbind(data.pred, data.frame(Year = as.factor(Year[i]), pred.dens.sp[,1:3], fate = pred.dens.fate[,3], pred.dens=pred.dens.sp[,4], pred.env=pred.env.sp[,4]))
  
  estimate.env.tmp <- array(NA, length(Var.env))
  p.env.tmp <- array(NA, length(Var.env))
  
  estimate.env.tmp[match(estimates.env.red$Variable, Var.env)] <- as.numeric(estimates.env.red$Estimate)
  p.env.tmp[match(estimates.env.red$Variable, Var.env)] <- ifelse(as.numeric(estimates.env.red$Pr...t..) < 0.05, 1, 0.5)
  
  estimate.env <- c(estimate.env, estimate.env.tmp)
  p.env <- c(p.env, p.env.tmp)
   
  estimate.dens.tmp <- array(NA, length(Var.dens))
  p.dens.tmp <- array(NA, length(Var.dens))
  
  estimate.dens.tmp[match(estimates.dens.red$Variable, Var.dens)] <- as.numeric(estimates.dens.red$Estimate)
  p.dens.tmp[match(estimates.dens.red$Variable, Var.dens)] <- ifelse(as.numeric(estimates.dens.red$Pr...t..) < 0.05, 1, 0.5)
  
  estimate.dens <- c(estimate.dens, estimate.dens.tmp)
  p.dens <- c(p.dens, p.dens.tmp)
  
  a_ctrl <- aov(pred~Sp, data=pred.env.sp[pred.env.sp$Treat=="0%",])
  tHSD_ctrl <- TukeyHSD(a_ctrl, ordered = FALSE, conf.level = 0.95)
  
  a_thnn <- aov(pred~Sp, data=pred.env.sp[pred.env.sp$Treat=="30%",])
  tHSD_thnn <- TukeyHSD(a_thnn, ordered = FALSE, conf.level = 0.95)
  
  for (j in Plot) {
    
    estimate.env.tmp <- array(NA, length(Var.env))
    p.env.tmp <- array(NA, length(Var.env))
    
    estimate.env.tmp[match(rownames(var.env.Plot[[j]]), Var.env)] <- as.numeric(var.env.Plot[[j]]$Estimate)
    p.env.tmp[match(rownames(var.env.Plot[[j]]), Var.env)] <- var.env.Plot[[j]]$Ztest
    names(estimate.env.tmp) <- Var.env
    
    estimate.env.plot <- rbind(estimate.env.plot, data.frame(Year = Year[i], Plot = j, data.frame(t(estimate.env.tmp))))
    p.env.plot <- rbind(p.env.plot, data.frame(Year = Year[i], Plot = j, data.frame(t(p.env.tmp))))
    
    estimate.dens.tmp <- array(NA, length(Var.dens))
    p.dens.tmp <- array(NA, length(Var.dens))
    
    estimate.dens.tmp[match(rownames(var.dens.Plot[[j]]), Var.dens)] <- as.numeric(var.dens.Plot[[j]]$Estimate)
    p.dens.tmp[match(rownames(var.dens.Plot[[j]]), Var.dens)] <- var.dens.Plot[[j]]$Ztest
    names(estimate.dens.tmp) <- Var.dens
    
    estimate.dens.plot <- rbind(estimate.dens.plot, data.frame(Year = Year[i], Plot = j, data.frame(t(estimate.dens.tmp))))
    p.dens.plot <- rbind(p.dens.plot, data.frame(Year = Year[i], Plot = j, data.frame(t(p.dens.tmp))))
    
  }
  
  
  #https://stackoverflow.com/questions/18771516/is-there-a-function-to-add-aov-post-hoc-testing-results-to-ggplot2-boxplot
  
  list.env.sp[[i]] <- ggplot(sum.env.sp, aes(x=Treat, y=mean, fill=Sp)) + 
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    labs(x="Treatment", y="Probability") + 
    geom_errorbar( aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) 
    #geom_text(data = generate_label_df(pred.env.sp, tHSD_ctrl, 'Sp'), aes(x = plot.labels, y = V1, label = labels))
    #title("Best environmental model for each plot", line = -1, cex.main = 1, outer = TRUE)
  
  list.dens.sp[[i]] <- ggplot(sum.dens.sp, aes(x=Treat, y=mean, fill=Sp)) + 
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    labs(x="Treatment", y="Probability") + 
    geom_errorbar( aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9))
    #title("Best environmental model for each plot", line = -1, cex.main = 1, outer = TRUE)
  
  
  if (i != length(Year)) {
    
    signif.num <- function(x) {
      symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
             cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
             symbols = c("***", "**", "*", ".", "n.s"))
    }
    
    sum.env.fate <- sum.env.fate %>% mutate(Treat = ifelse( Treat == "30%", "Thnn", "Ctrl" ))
    mod.env <- aov(pred.env ~ Treat: fate, data = data.pred[data.pred$Year==Year[i],])
    tukey.env <- data.frame(TukeyHSD(mod.env)$`Treat:fate`, label = array(signif.num(TukeyHSD(mod.env)$`Treat:fate`[,4])))
    tukey.env <- tukey.env[c(2,5),]
    max.env <- max(sum.env.fate$mean + sum.env.fate$sd)
    
    list.env.fate[[i]] <- ggplot(sum.env.fate, aes(x=sum.env.fate$fate, y=mean, fill=sum.env.fate$fate)) + 
      geom_bar(stat="identity", color="black", position=position_dodge()) +
      facet_grid(.~Treat) + 
      labs(x="", y=expression(paste("Pred. intensity ", lambda))) + 
      scale_y_continuous(limits=c(0,max.env)) + 
      geom_errorbar( aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) +
      geom_signif(annotations = c("",""), y_position = max.env-0.01, xmin=c(1,1), xmax=c(2,2)) +
      geom_text(x = c(1.5, 2.5, 1.5, 2.5), y = rep(max.env, 4), aes(label = c(tukey.env$label[1], "",tukey.env$label[2], "")), data = tukey.env) + 
      ggtitle(paste0("Environmental factors // ", Year[i])) +
      theme(plot.title = element_text(size=10), axis.text.x = element_blank())
      #title("Best environmental model for each plot", line = -1, cex.main = 1, outer = TRUE)
    
    sum.dens.fate <- sum.dens.fate %>% mutate(Treat = ifelse( Treat == "30%", "Thnn", "Ctrl" ))
    mod.dens <- aov(pred.dens ~ Treat: fate, data = data.pred[data.pred$Year==Year[i],])
    tukey.dens <- data.frame(TukeyHSD(mod.dens)$`Treat:fate`, label = array(signif.num(TukeyHSD(mod.dens)$`Treat:fate`[,4])))
    tukey.dens <- tukey.dens[c(2,5),]
    max.dens <- max(sum.dens.fate$mean + sum.dens.fate$sd/4)
    
    list.dens.fate[[i]] <- ggplot(sum.dens.fate, aes(x=sum.dens.fate$fate, y=mean, fill=sum.dens.fate$fate)) + 
      geom_bar(stat="identity", color="black", position=position_dodge()) +
      facet_wrap(.~Treat) + 
      labs(x="", y=expression(paste("Pred. intensity ", lambda))) + 
      scale_y_continuous(limits=c(0, max.dens)) + 
      geom_errorbar( aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) +
      geom_signif(annotations = c("",""), y_position = max.dens-0.03, xmin=c(1,1), xmax=c(2,2)) +
      geom_text(x = c(1.5, 2.5, 1.5, 2.5), y = rep(max.dens, 4), aes(label = c(tukey.dens$label[1], "",tukey.dens$label[2], "")), data = tukey.dens) + 
      ggtitle(paste0("Density-dependent factors // ", Year[i])) +
      theme(plot.title = element_text(size=10), axis.text.x = element_blank())
      #title("Best Aggregation model for each plot", line = -1, cex.main = 1, outer = TRUE)
    
  }
  
}


estimate.env.plot <- estimate.env.plot[order(estimate.env.plot$Plot),]
estimate.dens.plot <- estimate.dens.plot[order(estimate.dens.plot$Plot),]

p.env.plot <- p.env.plot[order(p.env.plot$Plot),]
p.dens.plot <- p.dens.plot[order(p.dens.plot$Plot),]

p.env.plot <-  ifelse(p.env.plot[,3:17] == 1, levels(var.env.Plot[[j]]$Ztest)[[1]],
       ifelse(p.env.plot[,3:17] == 2, levels(var.env.Plot[[j]]$Ztest)[[2]],
             ifelse(p.env.plot[,3:17] == 3, levels(var.env.Plot[[j]]$Ztest)[[3]], 
                    ifelse(is.na(p.env.plot[,3:17]), "", " "))))

p.dens.plot <-  ifelse(p.dens.plot[,3:9] == 1, levels(var.dens.Plot[[j]]$Ztest)[[1]],
                      ifelse(p.dens.plot[,3:9] == 2, levels(var.dens.Plot[[j]]$Ztest)[[2]],
                             ifelse(p.dens.plot[,3:9] == 3, levels(var.dens.Plot[[j]]$Ztest)[[3]], 
                                    ifelse(is.na(p.dens.plot[,3:9]), "", " "))))

colnames(estimate.env.plot)[3:17] <- Var.env
colnames(estimate.dens.plot)[3:9] <- Var.dens

estimate.env.plot <- format(estimate.env.plot, trim = TRUE, digits=3, scientific= 5)
estimate.dens.plot <- format(estimate.dens.plot, trim = TRUE, digits=3, scientific= 5)

estimate.env.plot2 <- matrix(paste(as.matrix(estimate.env.plot[,3:17]), as.matrix(p.env.plot)), ncol = 15)
estimate.dens.plot2 <- matrix(paste(as.matrix(estimate.dens.plot[,3:9]), as.matrix(p.dens.plot)), ncol = 7)

colnames(estimate.env.plot2) <- Var.env
colnames(estimate.dens.plot2) <- Var.dens

estimate.env.plot2[estimate.env.plot2 == "NA NA"] <- ""
estimate.dens.plot2[estimate.dens.plot2 == "NA NA"] <- ""

estimate.env.plot2 <- cbind(estimate.env.plot[,2:1], estimate.env.plot2)
estimate.dens.plot2 <- cbind(estimate.dens.plot[,2:1], estimate.dens.plot2)


if (save.output == T) {
  
  write.table(estimate.env.plot2, sep = ",",
            paste0("~/Documentos/Datos NO publicados/BioIntForest/PPA_Results/Figures/estimates.env.Plot.txt"))
  
  write.table(estimate.dens.plot2, sep = ",",
            paste0("~/Documentos/Datos NO publicados/BioIntForest/PPA_Results/Figures/estimates.dens.Plot.txt"))
  
}
  
sum.env <- data.frame(Variable = Var.env, Year = rep(Year, each=length(Var.env)), estimate = estimate.env, p = p.env)

sum.dens <- data.frame(Variable = Var.dens, Year = rep(Year, each=length(Var.dens)), estimate = estimate.dens, p = p.dens)





if (save.output == T) pdf("~/Documentos/Datos NO publicados/BioIntForest/PPA_Results/Figures/PPA_estimates.pdf", width=5, height=8)  

ggballoonplot(sum.env, x = "Year", y = "Variable", fill = "estimate", size = "p") + 
  scale_fill_viridis_c(option = "C")

ggballoonplot(sum.dens, x = "Year", y = "Variable", fill = "estimate", size = "p") + 
  scale_fill_viridis_c(option = "C")

if (save.output == T) dev.off()


if (save.output == T) pdf("~/Documentos/Datos NO publicados/BioIntForest/PPA_Results/Figures/PPA_estimates.sp.pdf", width=6, height=9)  

prow <- plot_grid(list.env.sp[[1]] + theme(legend.position="none"), list.dens.sp[[1]] + theme(legend.position="none"), list.env.sp[[2]]+ theme(legend.position="none"), list.dens.sp[[2]]+ theme(legend.position="none"),
                  list.env.sp[[3]] + theme(legend.position="none"),  list.dens.sp[[3]] + theme(legend.position="none"), list.env.sp[[4]] + theme(legend.position="none"), list.dens.sp[[4]] + theme(legend.position="none"),
                  align = 'vh', ncol = 2,
                  labels = c("a)", "b)", "c)", "d)", "e)", "f)", "h)", "i)"),
                  hjust = -1)

# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  list.env.sp[[1]] + theme(legend.box.margin = margin(0, 0, 0, 12))
)

# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
print(
  plot_grid(prow, legend, rel_widths = c(3, .4)))


if (save.output == T) dev.off()


if (save.output == T) pdf("~/Documentos/Datos NO publicados/BioIntForest/PPA_Results/Figures/PPA_estimates.fate.pdf", width=6, height=9)  

prow <- plot_grid(list.env.fate[[1]] + theme(legend.position="none"), list.dens.fate[[1]] + theme(legend.position="none"),
                  list.env.fate[[2]] + theme(legend.position="none"),  list.dens.fate[[2]] + theme(legend.position="none"), 
                  list.env.fate[[3]] + theme(legend.position="none"), list.dens.fate[[3]] + theme(legend.position="none"),
                  align = 'vh', ncol = 2,
                  labels = c("a)", "b)", "c)", "d)", "e)", "f)"),
                  hjust = -1)

# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  list.env.fate[[1]] + theme(legend.box.margin = margin(0, 0, 0, 12))
)

# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
# print(
#   plot_grid(prow, legend, rel_widths = c(3, .4)))

plot(prow)

if (save.output == T) dev.off()

data.pred$pred.dens <- ifelse (data.pred$pred.dens > 1, 1, data.pred$pred.dens)




library(lme4)
library(lmerTest)
library(MuMIn)
library(sjPlot)
library(interactions)

model.dens <- glmer(pred.env ~ pred.env + Treat + Year + Sp + pred.env:Sp + pred.env:Treat + pred.env:Year + (1 | Plot), family=gaussian, data = data.pred)
summary(model.dens)
(aov.dens <- anova(model.dens))

model.env <- glmer(pred.dens ~ pred.dens +  Treat + Year + Sp + pred.dens:Sp + pred.dens:Treat + pred.dens:Year + (1 | Plot), family=gaussian, data = data.pred)
summary(model.env)
(aov.env <- anova(model.env))

model.fate <- glmer(fate ~ pred.dens + pred.env + Year + Treat + Sp + pred.dens:pred.env + pred.dens:Sp + pred.dens:Treat + pred.dens:Year + pred.env:Sp + pred.env:Treat + pred.env:Year + (1 | Plot), family=binomial, data = data.pred[!is.na(data.pred$fate),])
(aov <- anova(model.fate))
summary(model.fate)
lmerTest::rand(model.fate)

show_tests(model.fate, fractions = TRUE)$Product
step(kk)

# https://cran.r-project.org/web/packages/interactions/vignettes/interactions.html
interact_plot(model.fate, pred = pred.env, modx = pred.dens)
interact_plot(model.fate, pred = pred.dens, modx = Treat, interval = TRUE, int.width = 0.8)
interact_plot(model.fate, pred = pred.env, modx = Treat, interval = TRUE, int.width = 0.8)
