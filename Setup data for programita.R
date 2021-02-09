##########################################
####   PREPARE DATA FOR PROGRAMITA  ######
##########################################
library(shapefiles)

#capa con toda la información
setwd("~/Documentos/Datos NO publicados/BioIntForest/Data")
Data<- data.frame(read.table(file="cruce_sombras_plantulas_2008_2011_JVR.csv", header=T, sep=";",dec=","))
#DataB$MAPX<- round(DataB$MAPX,3); DataB$MAPY<- round(DataB$MAPY,3)

Sc<-1#00000#factor to convert to metes

#parcelas
PLOT <- unique(Data$PARCELA) #3, 4 y 9= 0%; 2, 5, 7= 30% 
#unique(Data$SP)
Data$SP2 <- ifelse(Data$SP=="Qi" | Data$SP=="Qir" | Data$SP=="Qic" | Data$SP=="Qi?", "Qi",
                   ifelse(Data$SP=="Qh" | Data$SP=="Qhr", "Qh",
                   ifelse(Data$SP=="Fs" | Data$SP=="Fsr", "Fs",NA)))
Data$SP2<-as.factor(Data$SP2)
Data$RES <- ifelse(Data$SP=="Qir" | Data$SP=="Qhr" | Data$SP=="Fsr", 1, 0)
SP <- unique(Data$SP2)[1:3] #Eliminar NA
sizejuv <- 50 #umbral para discriminar juvenil de adulto
sizead <- 200 #umbral para discriminar adulto

for (i in 1:length(PLOT)){ #se para porque en algunos casos no hay mortalidad (la seleccion de lineas es cero)
  DataB<-subset(Data, PARCELA==PLOT[i] & SP2!="<NA>") #& RES==0 
  #plot(DataB$MAPX,DataB$MAPY,col=as.integer(DataB$SP2)) Todos
  #calculate header (limits and counts) based on area sample area
  XMin <- min(DataB$MAPX); XMax <- max(DataB$MAPX)
  YMin <- min(DataB$MAPY); YMax <- max(DataB$MAPY) #pay attention we remove negative values (Y axis flipped)
  
  #select data
  selJuv08<-which((DataB$ALT2008<sizejuv & DataB$ALT2008>0))
  selJuv09<-which((DataB$ALT2009<sizejuv & DataB$ALT2009>0))
  selJuv10<-which((DataB$ALT2010<sizejuv & DataB$ALT2010>0))
  selJuv11<-which((DataB$ALT2011<sizejuv & DataB$ALT2011>0))
  selAd<-which((DataB$ALT2008==-9999.0 | DataB$ALT2008>=sizead) | (DataB$ALT2009==-9999.0 | DataB$ALT2009>=sizead) |
               (DataB$ALT2010==-9999.0 | DataB$ALT2010>=sizead) | (DataB$ALT2011==-9999.0 | DataB$ALT2011>=sizead)) 
    
  hist(DataB$ALT2008[selJuv08]);hist(DataB$ALT2009[selJuv09])
  hist(DataB$ALT2010[selJuv10]);hist(DataB$ALT2011[selJuv11])
  
  plot(DataB$MAPX[selJuv10],DataB$MAPY[selJuv10],col=as.integer(DataB$SP2[selJuv10]))
  points(DataB$MAPX[selAd],DataB$MAPY[selAd],pch=21,col=1,bg=as.integer(DataB$SP2[selAd]))
  table(DataB$SP2[selJuv10]);table(DataB$SP2[selAd])
  
  #construct tables
  Data08<-rbind(cbind((DataB$MAPX[selAd]-XMin)*Sc,(DataB$MAPY[selAd]-YMin)*Sc,1,0),
                      cbind((DataB$MAPX[selJuv08]-XMin)*Sc,(DataB$MAPY[selJuv08]-YMin)*Sc,0,1))
  Data09<-rbind(cbind((DataB$MAPX[selAd]-XMin)*Sc,(DataB$MAPY[selAd]-YMin)*Sc,1,0),
                cbind((DataB$MAPX[selJuv09]-XMin)*Sc,(DataB$MAPY[selJuv09]-YMin)*Sc,0,1)) 
  Data10<-rbind(cbind((DataB$MAPX[selAd]-XMin)*Sc,(DataB$MAPY[selAd]-YMin)*Sc,1,0),
                cbind((DataB$MAPX[selJuv10]-XMin)*Sc,(DataB$MAPY[selJuv10]-YMin)*Sc,0,1)) 
  Data11<-rbind(cbind((DataB$MAPX[selAd]-XMin)*Sc,(DataB$MAPY[selAd]-YMin)*Sc,1,0),
                cbind((DataB$MAPX[selJuv11]-XMin)*Sc,(DataB$MAPY[selJuv11]-YMin)*Sc,0,1)) 
      
  #prepare headers
  nr08<-length(Data08[,1]); nr09<-length(Data09[,1]) #number of points
  nr10<-length(Data10[,1]); nr11<-length(Data11[,1]) 
  header08 <-paste("0.0",(XMax-XMin)*Sc,"0.0",(YMax-YMin)*Sc,nr08, sep=" ")
  header09 <-paste("0.0",(XMax-XMin)*Sc,"0.0",(YMax-YMin)*Sc,nr09, sep=" ")
  header10 <-paste("0.0",(XMax-XMin)*Sc,"0.0",(YMax-YMin)*Sc,nr10, sep=" ")
  header11 <-paste("0.0",(XMax-XMin)*Sc,"0.0",(YMax-YMin)*Sc,nr11, sep=" ")
  colnames(Data08)<-c(header08," "," "," "); colnames(Data09)<-c(header09," "," "," ") 
  colnames(Data10)<-c(header10," "," "," "); colnames(Data11)<-c(header11," "," "," ")     
  
  #save stuff
  setwd("~/Documentos/Datos NO publicados/PP Aspurz/PPA") #change folder if necessary
  write.table(Data08,paste("Plot",PLOT[i],"_JuvAll08.dat",sep=""),sep="\t",row.names=F,quote=F,eol = "\r\n") #programita data names the same
  write.table(Data09,paste("Plot",PLOT[i],"_JuvAll09.dat",sep=""),sep="\t",row.names=F,quote=F,eol = "\r\n") 
  write.table(Data10,paste("Plot",PLOT[i],"_JuvAll10.dat",sep=""),sep="\t",row.names=F,quote=F,eol = "\r\n") 
  write.table(Data11,paste("Plot",PLOT[i],"_JuvAll11.dat",sep=""),sep="\t",row.names=F,quote=F,eol = "\r\n")
  
  for (j in 1:length(SP)){
    #selec adults and select alive (a) and dead (d) juveniles
    selAda<-which(DataB$SP2==SP[j] & ((DataB$ALT2008==-9999.0 | DataB$ALT2008>=sizead) | (DataB$ALT2009==-9999.0 | DataB$ALT2009>=sizead) |
                   (DataB$ALT2010==-9999.0 | DataB$ALT2010>=sizead) | (DataB$ALT2011==-9999.0 | DataB$ALT2011>=sizead))) 
    selJuv08a<-which((DataB$ALT2008<sizejuv & DataB$ALT2008>0 & DataB$SP2==SP[j] & DataB$ALT2009!="NA"))
    selJuv08d<-which((DataB$ALT2008<sizejuv & DataB$ALT2008>0 & DataB$SP2==SP[j] & is.na(DataB$ALT2009)))
    selJuv09a<-which((DataB$ALT2009<sizejuv & DataB$ALT2009>0 & DataB$SP2==SP[j] & DataB$ALT2010!="NA"))
    selJuv09d<-which((DataB$ALT2009<sizejuv & DataB$ALT2009>0 & DataB$SP2==SP[j] & is.na(DataB$ALT2010)))
    selJuv10a<-which((DataB$ALT2010<sizejuv & DataB$ALT2010>0 & DataB$SP2==SP[j] & DataB$ALT2011!="NA"))
    selJuv10d<-which((DataB$ALT2010<sizejuv & DataB$ALT2010>0 & DataB$SP2==SP[j] & is.na(DataB$ALT2011)))
    selJuv11a<-which((DataB$ALT2011<sizejuv & DataB$ALT2011>0 & DataB$SP2==SP[j])) #aqui no se tienen datos del siguiente año
    
    length(is.na(selJuv08d))
    
    #construct tables
    #selecciona los adultos y juveniles de cada especie
    Data08s<-rbind(cbind((DataB$MAPX[selAda]-XMin)*Sc,(DataB$MAPY[selAda]-YMin)*Sc,1,0),
                   cbind((DataB$MAPX[selJuv08a]-XMin)*Sc,(DataB$MAPY[selJuv08a]-YMin)*Sc,0,1),
                   cbind((DataB$MAPX[selJuv08d]-XMin)*Sc,(DataB$MAPY[selJuv08d]-YMin)*Sc,0,1))
    Data09s<-rbind(cbind((DataB$MAPX[selAda]-XMin)*Sc,(DataB$MAPY[selAda]-YMin)*Sc,1,0),
                   cbind((DataB$MAPX[selJuv09a]-XMin)*Sc,(DataB$MAPY[selJuv09a]-YMin)*Sc,0,1),
                   cbind((DataB$MAPX[selJuv09d]-XMin)*Sc,(DataB$MAPY[selJuv09d]-YMin)*Sc,0,1))
    Data10s<-rbind(cbind((DataB$MAPX[selAda]-XMin)*Sc,(DataB$MAPY[selAda]-YMin)*Sc,1,0),
                   cbind((DataB$MAPX[selJuv10a]-XMin)*Sc,(DataB$MAPY[selJuv10a]-YMin)*Sc,0,1),
                   cbind((DataB$MAPX[selJuv10d]-XMin)*Sc,(DataB$MAPY[selJuv10d]-YMin)*Sc,0,1)) 
    Data11s<-rbind(cbind((DataB$MAPX[selAda]-XMin)*Sc,(DataB$MAPY[selAda]-YMin)*Sc,1,0),
                   cbind((DataB$MAPX[selJuv11a]-XMin)*Sc,(DataB$MAPY[selJuv11a]-YMin)*Sc,0,1))
    
    #selecciona todos los adultos para calcular la mortalidad
    Data08d<-rbind(cbind((DataB$MAPX[selAd]-XMin)*Sc,(DataB$MAPY[selAd]-YMin)*Sc,1,0),
                   cbind((DataB$MAPX[selJuv08a]-XMin)*Sc,(DataB$MAPY[selJuv08a]-YMin)*Sc,0,1),
                   cbind((DataB$MAPX[selJuv08d]-XMin)*Sc,(DataB$MAPY[selJuv08d]-YMin)*Sc,0,0)) 
    Data09d<-rbind(cbind((DataB$MAPX[selAd]-XMin)*Sc,(DataB$MAPY[selAd]-YMin)*Sc,1,0),
                   cbind((DataB$MAPX[selJuv09a]-XMin)*Sc,(DataB$MAPY[selJuv09a]-YMin)*Sc,0,1),
                   cbind((DataB$MAPX[selJuv09d]-XMin)*Sc,(DataB$MAPY[selJuv09d]-YMin)*Sc,0,0))
    Data10d<-rbind(cbind((DataB$MAPX[selAd]-XMin)*Sc,(DataB$MAPY[selAd]-YMin)*Sc,1,0),
                   cbind((DataB$MAPX[selJuv10a]-XMin)*Sc,(DataB$MAPY[selJuv10a]-YMin)*Sc,0,1),
                   cbind((DataB$MAPX[selJuv10d]-XMin)*Sc,(DataB$MAPY[selJuv10d]-YMin)*Sc,0,0)) 
        
    #prepare headers
    nr08s<-length(Data08s[,1]); nr09s<-length(Data09s[,1]) #number of points
    nr10s<-length(Data10s[,1]); nr11s<-length(Data11s[,1]) 
    header08s <-paste("0.0",(XMax-XMin)*Sc,"0.0",(YMax-YMin)*Sc,nr08s, sep=" ")
    header09s <-paste("0.0",(XMax-XMin)*Sc,"0.0",(YMax-YMin)*Sc,nr09s, sep=" ")
    header10s <-paste("0.0",(XMax-XMin)*Sc,"0.0",(YMax-YMin)*Sc,nr10s, sep=" ")
    header11s <-paste("0.0",(XMax-XMin)*Sc,"0.0",(YMax-YMin)*Sc,nr11s, sep=" ")
    colnames(Data08s)<-c(header08s," "," "," "); colnames(Data09s)<-c(header09s," "," "," ") 
    colnames(Data10s)<-c(header10s," "," "," "); colnames(Data11s)<-c(header11s," "," "," ") 
    
    nr08d<-length(Data08d[,1]); nr09d<-length(Data09d[,1]) #number of points
    nr10d<-length(Data10d[,1])
    header08d <-paste("0.0",(XMax-XMin)*Sc,"0.0",(YMax-YMin)*Sc,nr08d, sep=" ")
    header09d <-paste("0.0",(XMax-XMin)*Sc,"0.0",(YMax-YMin)*Sc,nr09d, sep=" ")
    header10d <-paste("0.0",(XMax-XMin)*Sc,"0.0",(YMax-YMin)*Sc,nr10d, sep=" ")
    colnames(Data08d)<-c(header08d," "," "," "); colnames(Data09d)<-c(header09d," "," "," ") 
    colnames(Data10d)<-c(header10d," "," "," ")
        
    #save stuff
    write.table(Data08s,paste("Plot",PLOT[i],"_Juv",SP[j],"08.dat",sep=""),sep="\t",row.names=F,quote=F,eol = "\r\n") #programita data names the same
    write.table(Data09s,paste("Plot",PLOT[i],"_Juv",SP[j],"09.dat",sep=""),sep="\t",row.names=F,quote=F,eol = "\r\n") 
    write.table(Data10s,paste("Plot",PLOT[i],"_Juv",SP[j],"10.dat",sep=""),sep="\t",row.names=F,quote=F,eol = "\r\n") 
    write.table(Data11s,paste("Plot",PLOT[i],"_Juv",SP[j],"11.dat",sep=""),sep="\t",row.names=F,quote=F,eol = "\r\n")
    write.table(Data08d,paste("Plot",PLOT[i],"_JuvD",SP[j],"08.dat",sep=""),sep="\t",row.names=F,quote=F,eol = "\r\n")
    write.table(Data09d,paste("Plot",PLOT[i],"_JuvD",SP[j],"09.dat",sep=""),sep="\t",row.names=F,quote=F,eol = "\r\n") 
    write.table(Data10d,paste("Plot",PLOT[i],"_JuvD",SP[j],"10.dat",sep=""),sep="\t",row.names=F,quote=F,eol = "\r\n") 
  }
}