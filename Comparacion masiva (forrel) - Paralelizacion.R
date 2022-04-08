#cargo paquetes
packages<-c("pedtools","pedmut", "forrel", "dvir", "writexl","tidyverse","matrixStats","doParallel", "foreach", 
            "kinship2", "Matrix")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))

modifyID<-function(ped){
  ped$ID<-ped$ID%>%gsub("-.*?$","",.)
  return(ped)
}

fecha = format(Sys.Date(), format="%d-%m-%Y")
# si hay mas de un .fam agregar a la siguiente linea si es el [[1]], el [[2]], etc
fams = list.files()[grep("*.fam$",list.files())]
#cargo dvis. Se necesita que tenga los Unidentified person en el primer elemento
dvi = readFam(fams, useDVI = TRUE, verbose = TRUE)
xref<-map(dvi[2:length(dvi)], ~.x[['Reference pedigree']]%>%modifyID())
for (i in seq(1:length(xref))) {
  if (length(xref[[i]])!=9) {
    xref[[i]]=NA
  }
}

print(paste("Grupos con familias inconexas : ", paste(names (dvi[is.na(xref)]), collapse=" ;  ")))
xref=xref%>%discard(is.na(xref))
dvi2=append(list('Unidentified persons'=dvi[[1]]),xref)

################################################################################
#Forrel

missing = "Missing person"

#loop PI
n.cores <- parallel::detectCores() - 2

my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

#funcion
loops_GF<-function(dvi, pi){
  pm = dvi$'Unidentified persons'[[i]]
  pi = c(pm$ID)
  vector<-append(pi,vector("numeric",length=length(dvi)-1))
  for (j in seq(2,5)){
    pm2=pm
    gf<-dvi[[j]]
    markers<-match(name(gf), name(pm2))
    markers2<-match(name(pm2),name(gf))
    
    pm2$MARKERS<-pm2$MARKERS[markers[!is.na(markers)]]
    gf$MARKERS<-gf$MARKERS[markers2[!is.na(markers2)]]
    
    
    #calculo de LR
    lr = missingPersonLR(
      reference = gf,
      poi=pm2,
      missing)$LRtotal
    
    #apendeo el LR en un vector  
    vector[[j]] = lr
    
  }
  return(vector)
}

#paralelizacion
df<-foreach(i=1:2, .combine='rbind', .packages = c("forrel","pedtools"), .inorder=T)%do%{
  loops_GF(dvi = dvi2)
}



#nombre de los PI
names(df)=  c("PI", names(dvi)[2:ncol(df)])

df<-as.data.frame(df)

#transforma en numeros a los PI
for (m in 2:ncol(df)){
  df[,m]=as.numeric(df[,m])
}

#LR maximos arriba (matriz completa)
#maxdf=rowMaxs(as.matrix(df[,2:ncol(df)]), na.rm=T)
#df2=df[order(maxdf, decreasing = T),]

lista<-map(2:ncol(df), ~df[,c(1,.x )]%>%mutate(GF=names(df)[.x])%>%rename("LR"=names(df)[.x]))

#agrega nombres de los GF
lista<-map2(lista,names(dvi)[2:length(dvi)], ~.x%>%mutate(GF=.y) )


lista2<-map(lista, ~if(length(.x$LR[.x$LR>10000])>3){ .x=filter(.x, LR>10000)}
            else{ top_n(.x, 3, LR)})

output2<-bind_rows(lista2)%>%arrange(desc(LR) )
  
#exporto
write_xlsx(output2,paste0("BNDG-DVI (LR)","_", fecha, ".xlsx"))
#write_tsv(output2,"BNDG-DVI (LR).tsv")
