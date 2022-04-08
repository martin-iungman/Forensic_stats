#Este script genera archivo excel de tipo "legajos", para los jovenes seleccionados. Falta incluirle el # de mismatch. 

library(pedtools)
library(pedmut)
library(forrel)
library(dvir)
library(writexl)
library(tidyverse)
library(matrixStats)
library(doparallel)
library(foreach)

#Usar un archivo que incluya los jovenes 
dvi = readFam("BNDG-DVI.fam", useDVI = TRUE, verbose = TRUE)
#Restrinjo a las muestras de interes. En este caso los primeros 10
#Tambien se puede escribir con los nombres con dvi[[1]][c("D79790", "M39999")]. Se puede ver los nombres de los jovenes con names(dvi)[[1]]
dvi[[1]]<-dvi[[1]][1:10]

#Forrel

missing = "Missing person"

#df vacio
df=tibble()

#loop PI

n.cores <- parallel::detectCores() - 1

my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()


loops_GF<-function(dvi, umbral){
  pm = dvi$'Unidentified persons'[[i]]
  pi = c(pm$ID)
  df<-data.frame()
  for (j in seq(2,length(dvi))){
    pm2=pm
    markers<-match(name(dvi[[j]]$'Reference pedigree'), name(pm2))
    pm2$MARKERS<-pm2$MARKERS[markers]
    
    #calculo de LR
    lr = missingPersonLR(
      reference = dvi[[j]]$'Reference pedigree',
      poi=pm2,
      missing)$LRtotal
    vector<-c(names(dvi)[j], pi, lr, length(markers[!is.na(markers)]))
    #apendeo el LR en un vector  
    df<-rbind(df,vector)
    names(df)<-c("Grupo Familiar", "Muestra","LR","#Marcadores")
    
  }
  df$LR<-as.numeric(df$LR)
  df<-df%>%filter(LR>umbral)
  return(df)
}

tic()
df_list<-foreach(i=1:length(dvi[[1]]),.packages = c("forrel","pedtools", "dplyr","magrittr"), .inorder=T)%dopar%{
  loops_GF(dvi = dvi,umbral= 5)
}
toc()

names(df_list)<-map(dvi[[1]], ~.x$ID)
walk2(df_list, names(df_list),~write_xlsx(.x, paste0("tablas_legajos/",.y, "_comparaciones.xlsx")))



df<-as.data.frame(df)

#transforma en numeros a los PI
for (m in 2:ncol(df)){
  df[,m]=as.numeric(df[,m])
}

#LR maximos arriba
maxdf=rowMaxs(as.matrix(df[,2:ncol(df)]), na.rm=T)
df2=df[order(maxdf, decreasing = T),]

lista<-map(2:ncol(df), ~df[,c(1,.x )]%>%mutate(GF=names(df)[.x])%>%rename("lr"=names(df)[.x]))
lista2<-map(lista, ~if(length(.x$lr[.x$lr>10000])>3){ .x=filter(.x, lr>10000)}
            else{ top_n(.x, 3, lr)})
output2<-bind_rows(lista2)%>%arrange(desc(lr) )


#exporto
write_xlsx(df2,"BNDG-DVI (LR).xlsx")