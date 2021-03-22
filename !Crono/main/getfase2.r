
CEP_csv  <- read.csv2("C:\\Users\\milax\\Documents\\!Crono_November_2020\\Data base\\ActualDatabase.csv")


library(dplyr)        #for data frame functions
library(stringr)      #for scan / print functions
library(data.table)   # For like function (%like%)

Main<- function()numeric();
search_name<-function(i,elements = c(i) )character();#Devuelve un array
gibbs_to_eachr<-function(TK,TO,Pbar,p,Category,elements=c(i),i)numeric(); #funcion gibbs to each one




 #esto es main
Main<- function()
{
  
  aux<-0;   #leo los  elementos que deben ser incluidos en el compuesto
  
  for (j in 1:i) 
  {
    aux<-scan("", what = character(),1);   #leo los  elementos que deben ser incluidos en el compuesto
    if(!(CEP_csv %>% dplyr::filter(Name == aux)))  #Averiguo sieste elemento existe
    {
      printf("try again!");
      j<-j-1; #si el elemento no existe en la db , el Crono lo pide ingressar nuevamente
    }
  }
  #para energia libre de gibbs
  #valores unicos
  wG<-0;
  aqIG<-0;
  hG<-0;  
  #arrays a eleccion del querido usuario
  MinAr<-0;
  MinGAr<-0;
  AqAr<-0;
  AqGar<-0;
  GasAr<-0;
  GasGAr<-0;
  #Categoria que cambia para cada array
  Category<-0;
  #leo los valores para las variables de entorno
  Tk<-scan("", what = numeric(),1);
  Pbar<-scan("", what = numeric(),1);
  p<-scan("", what = numeric(),1);
  TC<-scan("", what = numeric(),1);
  gH2O<-scan("", what = numeric(),1);
  eTP<- scan("", what = numeric(),1);

  
  
  
  
  
  
  i<-0;  #numero de elementos para cada array
  
  i<-scan("", what = numeric(),1);
  
  
  MinAr<- c(i);  #arreglo para Min , 
  
  
  MinAr<-search_name(i,MinAr );#Funcion para leer el Array y a la vez averiguar la existencia  de sus elementos en la db

  MinGAr<-c(i);   
  Category<-'Min';
  gibbs_to_eachr(TK,TO,Pbar,p,Category,MinGAr,i);
  
  i<-scan("", what = numeric(),1);
  AqAr<- c(i);  #arreglo para Aq, 
  AqAr<-search_name(i,AqAr );#Funcion para leer el Array y a la vez averiguar la existencia  de sus elementos en la db
  Category<-'Aq';
  AqGAr<-c(i);
  gibbs_to_eachr(TK,TO,Pbar,p,Category,AqGAr,i);
  
  i<-scan("", what = numeric(),1);
  GasAr<- c(i);  #arreglo para Gas ,
  Category<-'Gas';
  GasAr<-search_name(i,GasAr );#Funcion para leer el Array y a la vez averiguar la existencia  de sus elementos en la db
  
  GasGAr<-c(i);
  gibbs_to_eachr(TK,TO,Pbar,p,Category,GasGAr,i);
  
 return(0); 
}


  search_name<-function(i )
{
    aux<-0;
    elements <- c(i);
  j<-0;
  
  for (j in 1:i) 
  {
    aux<-scan("", what = character(),1);   #leo los  elementos que deben ser incluidos en el compuesto
    elements[j]<-aux;
    if (elements[j]=="/") #fin de la cadena
    {
      break;
    }
    if(!(CEP_csv %>% dplyr::filter(Name == elements[j])))  #Averiguo sieste elemento existe
    {
      printf("try again!");
      j<-j-1; #si el elemento no existe en la db , el Crono lo pide ingressar nuevamente
    }
     
  }
 
 return(elements );
  
}



gibbs_to_eachr<-function(TK,TO,Pbar,p,Category,elements=c(i),i) #funcion gibbs to each one

{
  
  j<-0;
  Aux1<-0; #variables auxiliares para GibbsEnergy y GibbsEnergy uncerntainity respectivamente
  Aux2<-0;
  
  wG<-waterGibbsFreeEnergy('H2O',Tk,Pbar);  #Por default se calculan los valores para H2O , OH y H
  aqIG<-AqGibbsFreeEnergy('OH',Tk,Pbar); 
  hG<-GasGibbsFreeEnergy('H',Tk,Pbar); 
  
  
  #calcular energia libre de gibbs para cada elemento del array . 
  for(j in 1:i)
  {
    if(Category==Min)   #Cada array tiene solo una categoria . O es Min o es Aq o es Gas
    { 
      Aux1<- MinGibbsFreeEnergy(elements[j], Tk,Pbar);
      Aux2<- MinGibbsFreeEnergyUncertainty(elements[j], Tk)  ;
      MinGibbsAr[j] = Aux1;
      MinGibbsUniAr[j] = Aux2;
    }
    if(Category==Gas)
    {
      Aux1<-GasibbsFreeEnergy(elements[j], Tk,Pbar);
      Aux2<-GasGibbsFreeEnergyUncertainty(elements[j], Tk);  
      GasGibbsAr[j] = Aux1;
      GasGibbsUniAr[j] = Aux2;
    }
    if(Category==Aq)
    {
      Aux1<- AqGibbsFreeEnergy(elements[j], Tk,Pbar);
      Aux2<-  AqGibbsFreeEnergyUncertainty(elements[j], Tk);  
      AqGibbsAr[j] = Aux1;
      AqGibbsUniAr[j] = Aux2;  #guardo los resultados de energia libre de Gibbs normal /  uncertanity de cada componente de la data base en un array destinado
    }
    #los 2  arrays se declararan como variables globales
    
  }
  return(0);
}

 




end

