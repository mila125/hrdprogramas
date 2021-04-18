import::from("C:\\Users\\Mila\\Documents\\R-projects\\!Crono\\thermodynamics\\!Crono.R","waterGibbsFreeEnergy","MinGibbsFreeEnergy")
import::from("C:\\Users\\Mila\\Documents\\R-projects\\!Crono\\thermodynamics\\!Crono.R","GasGibbsFreeEnergyUncertainty","GasGibbsFreeEnergy");
import::from("C:\\Users\\Mila\\Documents\\R-projects\\!Crono\\thermodynamics\\!Crono.R","AqGibbsFreeEnergyUncertainty","AqGibbsFreeEnergy");
library(dplyr)        #for data frame functions
library(stringr)      #for scan / print functions
library(stringi)      #for scan / print functions
library(data.table)   # For like function (%like%)

Main<- function()numeric();
search_name<-function(i)data.frame();
gibbs_to_eachr<-function(TK,T0,Pbar,p,CEP_2,i)data.frame(); 
get_element_coefic<-function()data.frame();   #las cuatro funciones devuelven data frame , actualizan el CEP_2


CEP_csv<-read.csv("C:\\Users\\Mila\\Documents\\R-projects\\!Crono\\data_base\\ActualDatabase.csv");
CEP_2 <-read.csv("C:\\Users\\Mila\\Documents\\R-projects\\!Crono\\data_base\\SelectedDatabase.csv");


 #esT0 es main
Main<- function()
{
  



  #para energia libre de gibbs
  #valores unicos
  T0<-0;
  wG<-0;
  aqIG<-0;
  hG<-0;  
  #arrays a eleccion del querido usuario


  Ar<-0;



  #leo los valores para las variables de enT0rno
  
  print("Escribe las  variables de enT0rno del sistema \n");
  Tk<-scan("", what = numeric(),1);
  Pbar<-scan("", what = numeric(),1);
  p<-scan("", what = numeric(),1);
  
  TC<-scan("", what = numeric(),1);
  gH2O<-scan("", what = numeric(),1);
  eTP<- scan("", what = numeric(),1);

  
  
  
  
  
  
  print("CuanT0s elemenT0s tendra el array?   \n");
 
  
  i<-scan("", what = numeric(),1);
  
  

 
  CEP <- read.csv("C:\\Users\\Mila\\Documents\\R-projects\\!Crono\\data_base\\ActualDatabase.csv");
  
  CEP_2<-search_name(i);#Funcion para leer el Array y a la vez averiguar la existencia  de sus elemenT0s en la db
 
  CEP_2<-get_element_coefic();#Funcion para leer el Array y a la vez averiguar la existencia  de sus elemenT0s en la db

 
  gibbs_to_eachr(TK,T0,Pbar,p,i); #Relleno el campo del G
  
  
  
 return(0); 
}


  search_name<-function(i)
{
  
    CEP_2  <- read.csv( "C:\\Users\\Mila\\Documents\\R-projects\\!Crono\\data_base\\SelectedDataBase.csv") ;  #lo leo otra vez porque se pierde

    aux<-0;   #leo los  elemenT0s que deben ser incluidos en el compuesT0
    
    for (j in 1:i) 
    {
      print("Proximo compuesT0  \n");
      aux<-scan("", what = character(),1);   #leo los  elemenT0s que deben ser incluidos en el compuesT0
      
    
     
      
  
        CEP_2<- rbind(CEP_2,CEP_csv %>% dplyr::filter(Name==aux));
    

      
      
      
      
    }



 return(CEP_2);
   }
  
get_element_coefic<-function() #la llamo durante cada ciclo del for each en la funcion element coefic
  {
  formula<-0;
  row<-1;
  i<-0; #variable de referencia
  element <- substr(formula, i, i); 
  i<-1;
  j<-0;
  coefic <- substr(formula, j, j); 
  

   
  
  
   for(row in 1 : count(CEP_2)[[1]][[1]])
    {
i<-1;
     
     formula<-CEP_2$Formula[row];
     
     #la formula del compuesT0 en la fila acc_row
     
     while(i <= stri_length(formula))   #i no cambia , pero se ira cambiando en los dos for
     {
       for(i in i:stri_length(formula))  # el i va de un en uno hasta llegar a (  ,  que es inicio del numero que es la cantid de dicho elemenT0 en el compuesT0
       {
         
         
         if(substr(formula, i, i)=="(")  #termina el elemenT0
         {
           
           element <-gsub("0", "",element) ; 
           print(element);
           if(!(element %in% colnames(CEP_2) ) )
           {
             CEP_2[element] <- NA;
            
           }
           break;     #solo el numero
         }
         
         element <-   gsub(" ", "", paste(element , substr(formula, i, i))) ; 
         
         
         
         
       }
       
       j<-i+1;   #el j da un salT0 y empieza desde la posicion del numero
       
       for(j in  j:stri_length(formula))  # el j va de un en uno hasta llegar a )  ,  que es inicio del proximo  elemenT0 presente en el compuesT0
       {
         
         
         if(substr(formula, j, j)==")")  #termina el numero
         {
           print(coefic);
           CEP_2[row,as.character(element)]<-coefic;
           element<-0;
           coefic<-0;
           break;     #solo el numero
         }
         coefic <- gsub(" ", "",paste(coefic ,substr(formula, j, j)));
         
         
       }
       
       i<-j+1;#in next position #el i da un salT0 y empieza desde la posicion del elemenT0
       #CEP_2[acc_row][as.character(element)] <-coefic;
       #la escribe coef en la columna element  y fila acc_row
     }
     #proximo compuesT0 de la lsta   
     
     
      
    }
    
      
      
      
      
  
  
    
   
  return (CEP_2);   #actualizo CEP_2
}
    


 gibbs_to_eachr<-function(Tk,T0,Pbar,p,CEP_2,i) #funcion gibbs T0 each one

{
   Gibbs<-"Gibbs";
   Gibbs_un<-"Gibbs_un";#col names T0 Gibbs/Gibbs unsertainty 
   CEP_2[Gibbs] <- NA;
   CEP_2[Gibbs_un] <- NA;
  j<-0;
  Aux1<-0; #variables auxiliares para GibbsEnergy y GibbsEnergy uncerntainity respectivamente
  Aux2<-0;

 
  
  #calcular energia libre de gibbs para cada elemenT0 del array . 
  for(j in 1:i)
  {
    if(CEP_2$Category[[j]]=="Min")   #Si la  categoria  es Min o es Aq o es Gas
    { 
      Aux1<-CEP_2$Name[[j]];
      Aux2<-Aux1;
      Aux1<- MinGibbsFreeEnergy(Aux1, Tk,T0,Pbar);
      Aux2<- MinGibbsFreeEnergyUncertainty(elements[j], Tk)  ;
      
      CEP_2[j,as.character(Gibbs)]<-Aux1;
      CEP_2[j,as.character(Gibbs_un)]<-Aux2;
    }
    if(CEP_2$Category[[j]]=="Gas")
    {
      Aux1<-CEP_2$Name[[j]];
      Aux2<-Aux1;
      Aux1<-GasGibbsFreeEnergy(Aux1, Tk,Pbar);
      Aux2<-GasGibbsFreeEnergyUncertainty(Aux2, Tk,Pbar);  

      CEP_2[j,as.character(Gibbs)]<-Aux1;
      CEP_2[j,as.character(Gibbs_un)]<-Aux2;
    } 
    if(CEP_2$Category[[j]]=="Aq")
    {
      Aux1<-CEP_2$Name[[j]];
      Aux2<-Aux1;
      Aux1<- AqGibbsFreeEnergy(Aux1, Tk,Pbar, gH2O,  eTP);
      Aux2<-  AqGibbsFreeEnergyUncertainty(elements[j], Tk,Pbar);  
   
      
      CEP_2[j,as.character(Gibbs)]<-Aux1;
      CEP_2[j,as.character(Gibbs_un)]<-Aux2;  #guardo los resultados de energia libre de Gibbs normal /  uncertanity de cada componente de la data base en un array destinado
    }
    #los 2  arrays se declararan como variables globales
    
  }

  return(CEP_2);
}

 




end

