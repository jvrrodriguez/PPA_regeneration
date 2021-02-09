
folder <- paste("~/Documentos/Datos NO publicados/BioIntForest/Data/Fotos_hemisfericas/")

Year <- 2008

hemiphot.data <- read.table(paste0(folder, "Hemiphot_Plots_", Year, ".csv"), sep = ",", header = T)

hemiview.data <- read.table(paste0(folder, "Luz Aspurz ", Year, ".csv"), sep = ",", header = T)

dim(hemiphot.data)
dim(hemiview.data)

hemiphot.data <- hemiphot.data[match(hemiview.data$Estaca, hemiphot.data$File), ]
hemiphot.data <- hemiphot.data[,-16]

cor(hemiphot.data[,5:16], hemiview.data[,2:12])