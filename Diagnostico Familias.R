packages<-c( "ibdsim2","ribd","furrr","tidyverse", "pedprobr", 
            "pedtools", "pedmut", "forrel", "poibin", "mispitools", "doParallel", "plotly", "htmlwidgets", "progressr","writexl")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

Trates_modif<-function (datasim, threshold) 
{
  nsims = nrow(datasim)
  TPED = datasim$Related
  RPED = datasim$Unrelated
  FPR = sum(RPED > threshold)/nsims
  FNR = sum(TPED < threshold)/nsims
  TPR = sum(RPED < threshold)/nsims
  TNR = sum(TPED > threshold)/nsims
  MCC = (TPR * TNR - FPR * FNR)/(sqrt(TPR + FPR) * sqrt(TPR + 
                                                          FNR) * sqrt(TNR + FPR) * sqrt(TNR + FNR))
  list(FNR=FNR, FPR=FPR, MCC=MCC)
}


DTsim_modif<-function (datasim, weight) 
{
  nsims = nrow(datasim)
  TPED = datasim$Related
  RPED = datasim$Unrelated
  ValoresLR = seq(1, nsims, length.out = nsims)
  FPs = 0
  FNs = 0
  for (i in 1:10000) {
    FPs[i] = sum(RPED > ValoresLR[i])
    FNs[i] = sum(TPED < ValoresLR[i])
  }
  Dis = sqrt(((FNs/nsims))^2 + ((weight * FPs/nsims)^2))
  Dis=Dis[!is.na(Dis)]
  Tabla = base::data.frame(x = ValoresLR, y = Dis)
  DT = which.min(Tabla$y)
  return(DT)
}


#Elegir Nombre de Directorio donde guardar toda la informacion. Si no existe se va a crear
file<-list.files()[grep("*.fam$",list.files())]

#cargar archivo .fam con las relaciones (sale del Familias3)
x = readFam(file, useDVI = T, verbose = T)
sexo_MP<-read_delim("sexo_MP.csv", delim= ";")[1:10,]

#creamos subcarpeta para guardar el output y lo seteamos como directorio de trabajo 
Dir<-tools::file_path_sans_ext(file)
Dir<-paste0(Dir,"_",Sys.Date())
if (!dir.exists(Dir)) {dir.create(Dir)}
setwd(Dir)

#Me quedo solo con el pedigree de referencia
xref<-map(x, ~.x[['Reference pedigree']])

#modifico nombres de la fllias
lista_familias<-ifelse(grepl("RM$",names(xref)), gsub("-.+","_RM",names(xref)),ifelse(
  grepl("\\)$", names(xref)), gsub("-.+\\(","\\(",names(xref)),gsub("-.+","",names(xref))))

names(xref)=lista_familias

#Creo df vacio
stats=matrix(nrow=length(xref), ncol=9)%>%as_tibble()
names(stats)=c("Familias", "EPtotal", "ExpMM", "log_LRmedio", "IP10000","LRdt", "FNR","FPR","MCC")
stats$Familias=lista_familias

#Numero de cores a usar
ncores=detectCores()-1


with_progress({
  p <- progressor(steps = length(xref))
  plan(multisession, workers = ncores)
  EP<-future_map(xref, ~{
    p()
    possibly(missingPersonEP, otherwise = NA_real_)(.x, missing="Missing person", verbose=F)
  })
})

stats$EPtotal<-map(EP, ~.x$EPtotal%>%round(2))%>%unlist()
stats$ExpMM<-map(EP, ~.x$expectedMismatch%>%round(2))%>%unlist()


#Simulaciones
correct_files=c()
for (flia in length(lista_familias)) {
  print(lista_familias[[flia]])
  print("Simulacion relacionados vs. no relacionados")
    timeStart<-Sys.time()
    datasim_sep<-list()
    Unrelated<-list()
    Related<-list()
    #La idea es separar las simulaciones en los ncores cada uno con un seed distinto
    seeds=seq(1900, 1900+ncores-1)
    #es importante que el nsim sea divisible por el ncores
    nsim=round(1000/ncores)*ncores
    
    plan(multisession, workers=ncores)
    datasim_sep<-future_map(seeds, ~makeLRsims(xref[[flia]], 
                                               missing ="Missing person", 
                                               nsim/ncores, 
                                               .x), 
                            future.rng.onMisuse="ignore", .options = furrr_options(seed=NULL))
    
    Unrelated<-map(datasim_sep, ~.x$Unrelated)%>%unlist()
    Related<-map(datasim_sep, ~.x$Related)%>%unlist()
    datasim <- data.frame(Unrelated,Related)
    
    timeEnd<-Sys.time()
    print(difftime(timeEnd, timeStart, units='mins')%>%round(2))
    

#Diagnostico general
    stats$log_LRmedio[flia]<-mean(log10(datasim$Related))%>%round(2)
    stats$IP10000[flia]<-sum(datasim$Related>10000)/nrow(datasim)
    LRdt<-DTsim_modif(datasim, 10)
    stats$LRdt[flia]<-LRdt
    rates<-Trates_modif(datasim, LRdt)
    stats$FNR[flia]<-rates$FNR
    stats$FPR[flia]<-rates$FPR
    stats$MCC[[flia]]<-rates$MCC
    ped<-xref[[flia]]
    ped$SEX<-unlist(replace(ped$SEX,ped$ID=="Missing person",ifelse(
      names(xref)[flia]%in%sexo_MP$`Grupo Familiar`, sexo_MP[which(sexo_MP$`Grupo Familiar`%in%names(xref)[flia]),3],3)))

  #Graficos
   try({ png(paste0(lista_familias[flia],"_diagnostico.png"), width=1920, height=1080)
    sub<-str_c(Sys.Date())
    layout(rbind(1:2))
    
    #pedigree
    ped=plot(
      ped,
      hatched = typedMembers, 
      margins = c(1,1,4,1),
      marker = xref[[lista_familias[flia]]]$MARKERS, 
      sep = "-", 
      missing = "", 
      cex=0.7
      )
    title(paste0(lista_familias[flia],"   ", sub), cex.main=1)
    
    #IBD
    if(length(typedMembers(xref[flia]))>2){
      ibd=checkPairwise(xref[[flia]], 
                    labels = TRUE, mar=c(1,.1,.1,.1)
                    )
    }
    dev.off()
    correct_files<-append(correct_files, lista_familias[[flia]])
   })
    datasim_long<-datasim%>%pivot_longer(cols=c(Related,Unrelated), names_to = "Simulacion", values_to = "LR")
    ggplot(datasim_long,aes(log10(LR), fill=Simulacion))+geom_density(alpha=0.7)+xlim(log10(min(datasim_long$LR))-2, log10(max(datasim_long$LR))+2)+geom_vline(xintercept=LRdt, linetype="dashed", col="darkred")+
      ggtitle(str_glue(lista_familias[flia], "   -   ",sub,"   -    nsim=", nsim), 
              subtitle=str_glue("\n","log10(LRmedio)=",as.numeric(stats$log_LRmedio[flia]), "\n", "Expected Mismatch=", as.numeric(stats$ExpMM[flia]), "\n","LRdt=", LRdt, "     ", "FPR=", rates$FPR, "     ", "FNR=", rates$FNR,"     ", "MCC=", rates$MCC))+
      theme(plot.title = element_text(hjust = 0.5, face="bold"), plot.subtitle = element_text(hjust = 0.5))
    ggsave(paste0(lista_familias[flia], "_density.png"))
    rm(datasim, datasim_long, ped, datasim_sep, Related, Unrelated)
    gc()
    #Distribucion de LR
    #htmlwidgets::saveWidget(as_widget(LRdist(datasim)),paste0(lista_familias[flia],"_LRdist.html"))
    #htmlwidgets::saveWidget(as_widget(deplot(datasim)),paste0(lista_familias[flia],"_deplot.html"))
    
}

print(paste0("Estas familas no se graficaron:", lista_familias[!lista_familias%in%correct_files]))

#Exporto
write_xlsx(stats, paste0("DIAGNOSTICO GENERAL - ",sub,".xlsx"))