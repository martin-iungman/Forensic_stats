#Este script es para comparar un joven contra una familia
#se requiere de un archivo .fam con las flias y los jovenes
#Inputs necesarios:
#id_joven
#flia_id

#carga o instalacion de paquetes (segun necesario)
packages<-c("pedtools","forrel", "pedmut", "pedprobr", "tidyverse")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))


#primero setear el working directory (funcion setwd o desde el panel de Files, yendo a la carpeta y setearla con "More" y "Set As Working Directory")
#cargo dvi:
fams = list.files()[grep("*.fam$",list.files())]
#(si hay mas de un archivo fam, seleccionar cual es el deseado con readFam(fams[[1]]) o el numero que corresponda) 
dvi = readFam(fams, useDVI = TRUE, verbose = TRUE)

# ID del joven, entre comillas
joven_id= "M0001"
  
# Nombre familia ( como esta cargado en fams, entre comillas - se puede mirar todos los nombre con names(fams))  
flia_id="Flia1"   
  
#edito los nombres:
modifyID<-function(ped){
  ped[['Reference pedigree']]$ID<-ped[['Reference pedigree']]$ID%>%gsub("\\[[^][]*]","",.)
  return(ped)
}

dvi<-map(dvi, ~.x%>%modifyID())

#SELECCION MUESTRA A COMPARAR:
jovenes=dvi$'Unidentified persons'
jovenes['Reference pedigree']<-NULL
ids<-map_chr(jovenes, ~.x$ID)
pi=jovenes[[names(ids[ids==joven_id])]]
familia=dvi[[flia_id]]$'Reference pedigree'

#plot hipotesis
missingPersonPlot(familia, missing = "Missing person", newdev = F)

#calculo de LR
lr = missingPersonLR(familia,'Missing person',pi)

#LR por marcador
lr$LRperMarker%>%as.data.frame()%>%
  rownames_to_column(var="Marker")%>%rename(LR=".")%>%
  ggplot(aes(Marker,log10(LR),fill=LR>1))+geom_bar(stat="identity")+
  #scale_y_log10()+
  theme(axis.text.x = element_text(angle=45), legend.position = "none")

#sensitivity plot: como se ve afectado el LR con la tasa mutacional. 
rvec = seq(0.000001, 
           0.1, 
           length= 100)

lrvec= sapply(rvec, function(r) {
H2 = setMutationModel(familia, model= "equal", rate = r)  
lr_aux = missingPersonLR(H2, "Missing person",pi)
lr_aux$LRtotal[1]
})

plot(log10(rvec),
     log10(lrvec),
     type = "l",
     xlab= "log10(rate)",
     ylab= " log10(LR)",
     main="Mutation Rate")
abline(v=log10(0.02), col="red")