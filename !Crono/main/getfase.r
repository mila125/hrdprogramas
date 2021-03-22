
CEP_csv  <- read.csv2("C:\\Users\\milax\\Documents\\!Crono_November_2020\\Data base\\ActualDatabase.csv")



library(gWidgets2)
library(tcltk)
library(gWidgets2tcltk)
library(dplyr)        #for data frame functions
library(import::here("!Crono.r"))        #for data frame functions
library(stringr)      #for scan / print functions
library(data.table)   # For like function (%like%)

Main<- function()numeric();
search<-function(i,elements = c(j) )data_frame();
gibbs_to_eachr<-function(TK,TO,Pbar,p)numeric(); #funcion gibbs to each one
Build_interface<- function()numeric();  



 #esto es main
Main<- function()
{
  
  Z<- 0;
  V<- 0;
  w<-0; 
  X<- 0 ;
  Y<-0;
  U<-0;
  A<- 0 ;
  B<-0;  
  C<-0;
  D<- 0;
  E<-0;
  
  
  #para energia libre de gibbs
  mG<-0;
  mGu<-0;
  aG<-0;
  aGu<-0;
  gG<-0;
  gGu<-0;
  
  #leo los valores para las variables de entorno
  Tk<-scan("", what = numeric(),1);
  Pbar<-scan("", what = numeric(),1);
  p<-scan("", what = numeric(),1);
  TC<-scan("", what = numeric(),1);
  gH2O<-scan("", what = numeric(),1);
  eTP<- scan("", what = numeric(),1);

  
  
  
  
  
  
  i<-0;
  i<-scan("", what = numeric(),1);
  elements <- c(i);  #arreglo para Min , 
  j<-0;
  
  for (j in 1:i) 
  {
    elements[j]<-scan("", what = character(),1);   #leo los  elementos que deben ser incluidos en el compuesto
    if (elements[j]=="/") #fin de la cadena
    {
      break;
    }
    CEP <- search_name(i,elements<- c(i)); 
  }
 
  

  
  gibbs_to_eachr(TK,TO,Pbar,p,CEP);
  
  return(0);
  
}

search_name<-function(i,elements = c(i) )  #la funcion recibe array de elementos requeridos como parametro
{
  j<-0;
  
  #hago una busqueada sucesssiva - primeramente busco por compuestos que contengan primer elemento del array , ya despues entre ellos por los que contengan tambien el segundo elemento e.t.c.
  CEP <- CEP_csv %>% dplyr::filter(Name %like% ""); #inicializo el cep
  
  for(j in 1:i)
  {
    if(elements[j]=="/")  #si ya no hay mas elementos requiridos , salir
    {
      break;
    }
    
    CEP <- CEP %>% dplyr::filter(Name %like% elements[j]); #si el nombre existe calcular e.g.
    
   }
  #View(CEP) #Visualización de la base
  
  return(CEP); #devuelve el CEP
}





search_form<-function(i,elements = c(i) )  #la funcion recibe array de elementos requeridos como parametro
{
  j<-0;
  
  #hago una busqueada sucesssiva - primeramente busco por compuestos que contengan primer elemento del array , ya despues entre ellos por los que contengan tambien el segundo elemento e.t.c.
  CEP <- CEP_csv %>% dplyr::filter(Formula %like% ""); #inicializo el cep
  
  for(j in 1:i)
  {
    if(elements[j]=="/")  #si ya no hay mas elementos requiridos , salir
    {
      break;
    }
    
    CEP <- CEP %>% dplyr::filter(Formula %like% elements[j]);
    
  }
  #View(CEP) #Visualización de la base
  
  return(CEP); #devuelve el CEP
}


gibbs_to_eachr<-function(TK,TO,Pbar,p,CEP) #funcion gibbs to each one

{
  
  j<-0;
  
  
  wG<-waterGibbsFreeEnergy(Tk,Pbar,p,i);  
  
  cols<-count(CEP,NULL,NULL);
  
  #calcular energia libre de gibbs para cada fila de las seleccionadas
  for(j in 1:cols$n)
  {
    if(CEP$Category==Min)   #numero de elementos en CEP
    { 
      mG<- MinGibbsFreeEnergy(CEP$Name[j], Tk,Pbar);
      mGu<- MinGibbsFreeEnergyUncertainty(CEP$Name[j], Tk)  ;
      valuesGibbs[i] = mG;
      valuesGibbsUni[i] = mGu;
    }
    if(CEP$Category==Gas)
    {
      gG-GasibbsFreeEnergy(CEP$Name[j], Tk,Pbar);
      gGu<-GasGibbsFreeEnergyUncertainty(CEP$Name[j], Tk);  
      valuesGibbs[i] = gG;
      valuesGibbsUni[i] = gGu;
    }
    if(CEP$Category==Aq)
    {
      aG<- AqGibbsFreeEnergy(CEP$Name[j], Tk,Pbar);
      aGu<-  AqGibbsFreeEnergyUncertainty(CEP$Name[j], Tk);  
      valuesGibbs[i] = aG;
      valuesGibbsUni[i] = aGu;  #guardo los resultados de energia libre de Gibbs normal /  uncertanity de cada componente de la data base en un array destinado
    }
    #los 2  arrays se declararan como variables globales
    
  }
  return(0);
}

 




end

