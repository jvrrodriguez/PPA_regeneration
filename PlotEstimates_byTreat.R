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

source("~/Documentos/Datos NO publicados/BioIntForest/Analysis/function quantumPlot.R")

Year <- c(2008, 2009, 2010, 2011)
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

for (i in 1:length(Year)){
  
  load(file=paste0("~/Documentos/Datos NO publicados/BioIntForest/Results/PPA_", Year[i],".RData"))
  
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
  
  #https://stackoverflow.com/questions/18771516/is-there-a-function-to-add-aov-post-hoc-testing-results-to-ggplot2-boxplot
  #esta pensada para trabajar con un boxplot y no para la grafica de barras
  generate_label_df <- function(d, HSD, flev){
    # Extract labels and factor levels from Tukey post-hoc 
    Tukey.levels <- HSD[[flev]][,4]
    Tukey.labels <- multcompLetters(Tukey.levels)['Letters']
    plot.labels <- names(Tukey.labels[['Letters']])
    
    # Get highest quantile for Tukey's 5 number summary and add a bit of space to buffer between    
    # upper quantile and label placement
    boxplot.df <- ddply(d, flev, function (x) max(fivenum(x$pred)) + 0.2)
    
    # Create a data frame out of the factor levels and Tukey's homogenous group letters
    plot.levels <- data.frame(plot.labels, labels = Tukey.labels[['Letters']],
                              stringsAsFactors = FALSE)
    
    # Merge it with the labels
    labels.df <- merge(plot.levels, boxplot.df, by.x = 'plot.labels', by.y = flev, sort = FALSE)
    
    return(labels.df)
  }
  
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
    
    list.env.fate[[i]] <- ggplot(sum.env.fate, aes(x=Treat, y=mean, fill=fate)) + 
      geom_bar(stat="identity", color="black", position=position_dodge()) +
      labs(x="Treatment", y="Probability") + 
      geom_errorbar( aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) + 
      geom_signif(comparisons = list(c("Fs", "Qh","Qi")), map_signif_level=TRUE)
    #title("Best environmental model for each plot", line = -1, cex.main = 1, outer = TRUE)
    
    list.dens.fate[[i]] <- ggplot(sum.dens.fate, aes(x=Treat, y=mean, fill=fate)) + 
      geom_bar(stat="identity", color="black", position=position_dodge()) +
      labs(x="Treatment", y="Probability") + 
      geom_errorbar( aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9))
    #title("Best Aggregation model for each plot", line = -1, cex.main = 1, outer = TRUE)
    
  }
  
}

sum.env <- data.frame(Variable = Var.env, Year = rep(Year, each=length(Var.env)), estimate = estimate.env, p = p.env)

sum.dens <- data.frame(Variable = Var.dens, Year = rep(Year, each=length(Var.dens)), estimate = estimate.dens, p = p.dens)

if (save.output == T) pdf("~/Documentos/Datos NO publicados/BioIntForest/Results/PPA_estimates.pdf", width=5, height=8)  

ggballoonplot(sum.env, x = "Year", y = "Variable", fill = "estimate", size = "p") + 
  scale_fill_viridis_c(option = "C")

ggballoonplot(sum.dens, x = "Year", y = "Variable", fill = "estimate", size = "p") + 
  scale_fill_viridis_c(option = "C")

if (save.output == T) dev.off()


if (save.output == T) pdf("~/Documentos/Datos NO publicados/BioIntForest/Results/PPA_estimates.sp.pdf", width=6, height=9)  

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


if (save.output == T) pdf("~/Documentos/Datos NO publicados/BioIntForest/Results/PPA_estimates.fate.pdf", width=6, height=9)  

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
print(
  plot_grid(prow, legend, rel_widths = c(3, .4)))

if (save.output == T) dev.off()

data.pred$pred.dens <- ifelse (data.pred$pred.dens > 1, 1, data.pred$pred.dens)

library(lme4)
library(lmerTest)
library(MuMIn)
library(sjPlot)
library(sjmisc)
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
