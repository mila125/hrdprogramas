
library(dplyr)        #for data frame functions
library(stringr)      #for scan / print functions
library(stringi)      #for scan / print functions
library(data.table)   # For like function (%like%)

Main<- function()numeric();
search_name<-function(total_rows)data.frame();
gibbs_to_eachr<-function(Tk,T0,Pbar,p,CEP_2,total_rows)data.frame(); 
                  
get_element_coefic<-function(total_rows)data.frame();   #las cuatro funciones devuelven data frame , actualizan el CEP_2
#thermodynamics

Born<- function( w, charge, Tk, Pbar, gH2O,  eTP,CEP_2)
SolvFunc <- function(p, TC, Pbar,CEP_2) 
waterDensity<- function(Tk,Pbar,p,CEP_2)
DielConst<- function(p, TC, Tk, Pbar,CEP_2)
  
waterGibbsFreeEnergy<- function(Tk,Pbar,p,CEP_2)   

GasGibbsFreeEnergy<- function(Gas,Tk,Pbar,CEP_2,row)
GasGibbsFreeEnergyUncertainty<- function(Gas, Tk,CEP_2,row)
  
AqGibbsFreeEnergy<- function(Aq,Tk, Pbar, gH2O,  eTP,CEP_2,row)
AqGibbsFreeEnergyUncertainty<- function(Aq, Tk, Pbar,CEP_2,row)



MinGibbsFreeEnergy<- function(Min,Tk,T0,Pbar,CEP_2,row)
MinGibbsFreeEnergyUncertainty<- function(Min, Tk,T0,CEP_2,row) 

CEP_csv  <- read.csv("C:\\Users\\Mila\\Documents\\R-projects\\!Crono\\data_base\\ActualDatabase.csv");


CEP_2<-0;

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
  T0<-scan("", what = numeric(),1);
  print("Escribe las  variables de enT0rno del sistema \n");
  Tk<-scan("", what = numeric(),1);
  Pbar<-scan("", what = numeric(),1);
  p<-scan("", what = numeric(),1);
  
  TC<-scan("", what = numeric(),1);
  gH2O<-scan("", what = numeric(),1);
  eTP<- scan("", what = numeric(),1);

  
  
  
  
  
  
  print("CuanT0s elemenT0s tendra el array?   \n");
 
  
  total_rows<-scan("", what = numeric(),1);
  
  

 
 
  
  CEP_2<-search_name(total_rows);#Funcion para leer el Array y a la vez averiguar la existencia  de sus elemenT0s en la db
  total_rows<-count(CEP_2);
  total_rows<-total_rows$n;
  
  CEP_2<-get_element_coefic(total_rows);#Funcion para leer el Array y a la vez averiguar la existencia  de sus elemenT0s en la db
  

  Gibbs<-"Gibbs";
  Gibbs_un<-"Gibbs_un";#col names T0 Gibbs/Gibbs unsertainty    
  
  CEP_2[Gibbs] <- NA;
  
  CEP_2[Gibbs_un] <- NA;  #creation of cols to store the gibbs energy comon/uncertainty
  
  total_rows<-count(CEP_2);
  total_rows<-total_rows$n;
  total_rows<-total_rows+1;
  
  water<-c(64) ; #num of cols include the elements and gibbs/gibbs_uncertainty
  
  
  rbind(CEP_2, water);
  
  
  #CEP_2$Name[total_rows]<-'water';  #inicialize 
 # CEP_2$Formula[total_rows]<-'H(2)O(1)';
#  CEP_2$H[total_rows]<-2;
#  CEP_2$O[total_rows]<-'1';
  #water<-waterGibbsFreeEnergy(Tk,Pbar,p,CEP_2,total_rows);
 # if(is.na(water))
#  {water<-0;}
 # CEP_2[total_rows,as.character(Gibbs)]<-water;
  gibbs_to_eachr(Tk,T0,Pbar,p,CEP_2,total_rows); #Relleno el campo del G
  

  
 return(0); 
}



  search_name<-function(total_rows)
{




    aux<-0;   
    
    for (j in 1:total_rows) 
    {
      print("Proximo compuesto  \n");
      aux<-scan("", what = character(),1);   #leo los  elemenT0s que deben ser incluidos en el compuesT0
      print("Cuanto quieres  \n");

     
   
  
        CEP_2<- rbind(CEP_2,CEP_csv %>% dplyr::filter(Name==aux));
    

      
      
      
      
    }



 return(CEP_2);
   }
  
get_element_coefic<-function(total_rows) #la llamo durante cada ciclo del for each en la funcion element coefic
  {
  formula<-0;
  row<-1;
  i<-0; #variable de referencia
  element <- substr(formula, i, i); 
  i<-1;
  j<-0;
  coefic <- substr(formula, j, j); 
  



   for(row in 1 : total_rows)
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
         
         
         if(substr(formula, j, j)==")") #termina el numero
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
    


 gibbs_to_eachr<-function(Tk,T0,Pbar,p,CEP_2,total_rows) #funcion gibbs T0 each one

{



  
   
   
   


  j<-0;
  Aux1<-0; #variables auxiliares para GibbsEnergy y GibbsEnergy uncerntainity respectivamente
  Aux2<-0;

 
  
  #calcular energia libre de gibbs para cada elemenT0 del array . 
  for(j in 1:total_rows)
  {
  if( isTRUE(CEP_2$Category[[j]]=="min") )   #Si la  categoria  es Min o es Aq o es Gas
    { 
      Aux1<-CEP_2$Name[[j]];
      Aux2<-Aux1;

      
      Aux1<- MinGibbsFreeEnergy(Aux1,Tk,T0,Pbar,CEP_2,j);
      Aux2<- MinGibbsFreeEnergyUncertainty(Aux2, Tk,T0,CEP_2,j)  ;
      if(is.na(Aux1))
      {
        Aux1<-0;
      }
      if(is.na(Aux2))
      {
        Aux2<-0;
      }
      CEP_2[j,as.character(Gibbs)]<-Aux1;
      CEP_2[j,as.character(Gibbs_un)]<-Aux2;
    }
   if( isTRUE(CEP_2$Category[[j]]=="gas"))
    {
      Aux1<-CEP_2$Name[[j]];
      Aux2<-Aux1;
      Aux1<-as.double(GasGibbsFreeEnergy(Aux1, Tk,Pbar,CEP_2,j));
      Aux2<-GasGibbsFreeEnergyUncertainty(Aux2, Tk,CEP_2,j);  

      if(is.na(Aux1))
      {
        Aux1<-0;
      }
      if(is.na(Aux2))
      {
        Aux2<-0;
      }
      
      CEP_2[j,as.character(Gibbs)]<-Aux1;
      CEP_2[j,as.character(Gibbs_un)]<-Aux2;
    } 
  if( isTRUE(CEP_2$Category[[j]]=="aq") )
    {
      Aux1<-CEP_2$Name[[j]];
      Aux2<-Aux1;

      Aux1<- AqGibbsFreeEnergy(Aux1, Tk,Pbar, gH2O,  eTP,CEP_2,j);
      Aux2<-  AqGibbsFreeEnergyUncertainty( Aux2,Tk,Pbar,CEP_2,j);  
                 
      if(is.na(Aux1))
      {
        Aux1<-0;
      }
      if(is.na(Aux2))
      {
        Aux2<-0;
      }
      CEP_2[j,as.character(Gibbs)]<-Aux1;
      
     CEP_2[j,as.character(Gibbs_un)]<-Aux2;  #guardo los resultados de energia libre de Gibbs normal /  uncertanity de cada componente de la data base en un array destinado
     }
    #los 2  arrays se declararan como variables globales
    
  }

  return(CEP_2);
}

 

#imported functions

 MinGibbsFreeEnergyUncertainty <- function(Min, Tk,T0,CEP_2,row)    #Calculates uncertainties of minerals
 {
   oG<-0;
   oS<-0; 
   oV<-0; 
   oa<-0; 
   ob<-0;

   
   
   


   
    oG<-as.numeric(CEP_2$SD..cal.mol.1[row]);
    oS<-as.numeric(CEP_2$S..cal.mol.1.?.K.1[row]);
  oV<-as.numeric(CEP_2$SD..cal.mol.1.?.K.1[row]);
   oa<-as.numeric(CEP_2$V..cm3.mol.1[row]);
   ob<-as.numeric(CEP_2$SD..cm3.mol.1[row]);
   oc<-as.numeric(CEP_2$a..cal.mol.1.?.K.1[row] );

     Result<-sqrt (oG^2) + (oS*(Tk-T0)^2)+ (oa*(Tk-T0-Tk));#*log(Tk/T0))^2) ; #+ (-ob/1000*(Tk-T0)*(Tk-T0)/2^2)+ (oc*100000*(Tk-T0)*(Tk-T0)/2/T0/T0/Tk^2)+ (oV*(Pbar-1)/41.84^2)); 

   
   
   
 }
 
 waterDensity <- function(Tk, Pbar,CEP_2,row)          #Calculates a density ofwater in g/cm3
 {  
   dyn.load("C:\\Users\\Mila\\Documents\\R-projects\\!Crono\\okawsp6.dll",TRUE,TRUE); 

   #from CronoMF 
   #  dyn.load("C:\\Users\\Mila\\Documents\\Rprojects\\Crono\\okawsp6.dll",TRUE,TRUE);
   wspDPT <- function (P, T)as.double; 
   
   
   
   if (Pbar<=1000) 
   {
     Result<- wspDPT(Pbar*1E+5, Tk)/1000;
   }
   else
   {
     Result<-(( wspDPT(99900000, Tk)*99900000- wspDPT(100000000, Tk)*99900000)/(-9990000000000)*Pbar*1E+5+( wspDPT(100000000, Tk)*99900000- wspDPT(99900000, Tk)*100000000)/(-100000))/1000;    #   {g/cm3} //Linear extrapolation 
   }  
   
 }
 
 SolvFunc <- function(p, TC, Pbar,CEP_2)           #Solvent function (g), Shock et al , 1992; Sue et al , 2002 - valid at ??>=0.20 and <=1, ??>=0.35 - Shock et al., 1992 (p. 809, Fig. 6)
 {
   g<-0;
   ag<-0;
   bg<-0;
   p_a<-0; 
   ag_i = -2.037662;         #??    Shock et al., 1992 (Table 3, p. 807)
   ag_ii = 0.005747;         #?? ??k     Shock et al., 1992 (Table 3, p. 807)
   ag_iii = -0.000006557892; #?? ??k^2  Shock et al., 1992 (Table 3, p. 807)
   bg_i = 6.107361;          #         Shock et al., 1992 (Table 3, p. 807)
   bg_ii = -0.01074377;      #1 ??C    Shock et al., 1992 (Table 3, p. 807)
   bg_iii = 0.00001268348;   #1 ??C^2  Shock et al., 1992 (Table 3, p. 807)
   ag_1 = 3.66666E-16;       #??k       Shock et al., 1992 (Table 4, p. 809)
   ag_2 = -1.504956E-10;     #?? bar^3  Shock et al., 1992 (Table 4, p. 809)
   ag_3 = 5.01799E-14;       #?? bar^4  Shock et al., 1992 (Table 4, p. 809)
   cg_1 = 0.18359;           #        Sue et al., 2002 (p. 3304)
   cg_2 = -0.18632;          #         Sue et al., 2002 (p. 3304)
   cg_3 = 0.11531;           #         Sue et al., 2002 (p. 3304)
   
   p_a<-p/ waterDensity(T0,1);
   p_a<-1;
   
   if (p_a>1 )
   { 
     p_a<-1;                                         #//Need T0 avoid the calculation collapse at T>25??C and provide some minor error on the Born function
     ag<-ag_i+ag_ii*TC+ag_iii*TC*TC;                               #//Shock et al., 1992 (Eq. 25, p. 807)
     bg<-bg_i+bg_ii*TC+bg_iii*TC*TC;                               #//Shock et al., 1992 (Eq. 26, p. 807)
     g<-ag* (1-p_a^ bg);  
   }
   #//T-P-dependent solvent function, Tanger and Helgeson, 1988; Shock et al., 1992 (Eq. 32, p. 809)
   if(((TC>=150) && (Pbar<=1000))&&((TC<=350) || (Pbar>=500)))  
   {  
     #//Shock et al., 1992 (p. 809, Fig. 6)
     
     #Shock et al., 1992 (p. 809, Fig. 6)
     
     g<-g-( ((TC-155)/300^ 4.8)+ag_1* ((TC-155)/300^ 16))*(ag_2* (1000-Pbar^3)+ag_3* (1000-Pbar^4)) ; #//Correction at T-P region II, Shock et al., 1992 (Eq. 32-33, p. 809)
   }
   else
   {
     g<-g-(cg_1*p_a+cg_2* (p_a^2)+cg_3*log(p_a));             #Sue et al., 2002 (Eq. 25, p. 3304) - valid at T - from 250 T0 600??C, Pbar - from 223 T0 998 bar, ?? - from 0.20 T0 0.81 g cm-3 (p. 3303)
   }
   
   
   Result<-g;   
   
 }  
 
 DielConst <- function(p, TC, Tk, Pbar,CEP_2)                      #//Dielectric constant ofwater ( e) at specified T and P
 {
   y<-0;
   c<-0;
   s<-0;
   d<-0;
   
   a_1 <- -1.576377E-3;    #1/??C                //Sverjensky et al., 2014 (Table 1, p. 131)
   a_2 <- 6.810288E-2;     #1/sqrt(??C)         //Sverjensky et al., 2014 (Table 1, p. 131)
   a_3 <- 7.548755E-1;     #                    //Sverjensky et al., 2014 (Table 1, p. 131)
   b_1 <- -8.016651E-5;    #1/??C                //Sverjensky et al., 2014 (Table 1, p. 131)
   b_2 <- -6.871618E-2;   #1/sqrt(??C)          //Sverjensky et al., 2014 (Table 1, p. 131)
   b_3 <- 4.747973;        #                    //Sverjensky et al., 2014 (Table 1, p. 131)
   c_1 <- 0.4117;                                #//Marshall, 2008 (p. 5); Dolejs, 2013 (p. 50)
   c_2 <- 366.6;                                 #//Marshall, 2008 (p. 5); Dolejs, 2013 (p. 50)
   c_3 <- -1.491E+5;                             #//Marshall, 2008 (p. 5); Dolejs, 2013 (p. 50)
   c_4 <- 9.190E+6;                              #//Marshall, 2008 (p. 5); Dolejs, 2013 (p. 50)
   d_1 <- 0.290;                                 #//Marshall, 2008 (p. 5); Dolejs, 2013 (p. 50)
   d_2 <- 275.4;                                 #//Marshall, 2008 (p. 5); Dolejs, 2013 (p. 50)
   d_3 <- 0.3245E+5;                             #//Marshall, 2008 (p. 5); Dolejs, 2013 (p. 50)
   s_1 <- 1.667;                                 #//Marshall, 2008 (p. 5); Dolejs, 2013 (p. 50)
   s_2 <- -11.41;                              #//Marshall, 2008 (p. 5); Dolejs, 2013 (p. 50)
   s_3 <- -3.526E+4;                             #//Marshall, 2008 (p. 5); Dolejs, 2013 (p. 50)
   
   if ((TC>100)&&(Pbar>1000)) 
   {
     Result<-exp(b_1*TC+b_2*sqrt(TC)+b_3)* (p^a_1*TC+a_2*sqrt(TC)+a_3) ;   #//Sverjensky et al., 2014 (Eq. 8, p. 131)
     
   }
   else
   {
     y<-1/(1+0.0012/p/p);                                       #//Marshall, 2008 (p. 4); Dolejs, 2013 (Eq. 64, p. 50)
     c<-c_1+c_2/Tk+c_3/Tk/Tk+c_4/Tk/Tk/Tk;                      #//Marshall, 2008 (p. 5); Dolejs, 2013 (Eq. 65, p. 50)
     d<-d_1+d_2/Tk+d_3/Tk/Tk;                                   #//Marshall, 2008 (p. 5); Dolejs, 2013 (Eq. 66, p. 50)
     s<-s_1+s_2/Tk+s_3/Tk/Tk;                                   #//Marshall, 2008 (p. 5); Dolejs, 2013 (Eq. 67, p. 50)
     Result<-1+ (10^ y*(c+(s-1)*log10(p))+d+log10(p));      #//Marshall, 2008 (Eq. 2, p. 5); Dolejs, 2013 (Eq. 63, p. 50)  
   }
 }
 
 
 
 
 MinGibbsFreeEnergy <- function(Min,Tk,T0,Pbar,CEP_2,row)  
 {  
   

   
   G<-as.numeric(CEP_2$G..cal.mol.1[row]);
   S<- as.numeric(CEP_2$SD..cal.mol.1[row]);
   
   
   V<- as.numeric(CEP_2$S..cal.mol.1.?.K.1[row]);
   
   a<- as.numeric(CEP_2$SD..cal.mol.1.?.K.1[row]);
   
   
   b<- as.numeric(CEP_2$V..cm3.mol.1[row]); 
   c<- as.numeric(CEP_2$SD..cm3.mol.1[row]);  
   
   A<-10;
   B<-40;
   
   
   Result <- G-S*(Tk-T0)+A*(Tk-T0-Tk*log10(Tk/T0))-(c*100000+B/1000*T0*T0*Tk)*(Tk-T0)*(Tk-T0)/2/T0/T0/Tk+V*(Pbar-1)/41.84; 
  return(Result);
 } 
 
 #other functions
 waterGibbsFreeEnergy<- function(Tk,Pbar,p,CEP_2,row) #Calculates the Gibbs energy ofwater from triple point Ptr and Ttr T0 10 kbar and 1273.15 ??k (Helgeson and kirkham, 1974; Hill, 1990; Johnson and NorT0n, 1991; Han, 2008)
 { 
   
   
   Rcal<-0;
   HH2O<-0;
   SH2O<-0;
   T<-0;
   k<-0;
   k<-0; 
   k0<-0;
   k1<-0;
   E<-0; 
   G<-0;
   H<-0;
   DT<-0;
   Dp<-0;
   UnitMinusE<-0;
   k0p<-0;
   kT<-0;
   k0T<-0;
   Z<-0;
   h<-0;
   F<-0;
   Ep<-0;
   FT<-0;
   Fh<-0;
   Fp<-0;
   Dkn<-0;
   DkT<-0;
   tmp<-0;
   w1<-0;
   w2<-0;
   w3<-0;
   w4<-0;
   w1T<-0;
   w2T<-0;
   w3T<-0;
   w4T<-0;
   R1<-0;
   R2<-0;
   R3<-0;
   R1p<-0;
   R2p<-0;
   R3p<-0;
   w1p<-0;
   w2p<-0;
   w3p<-0;
   w4p<-0;
   Gp<-0;
   GT<-0;
   HT<-0;
   dPdro<-0;
   I<-0;
   J<-0;
   A1<-matrix(0, nrow = 7, ncol = 7);  
   A2<-matrix(0, nrow = 7, ncol = 12);
   A3<-matrix(0, nrow = 5, ncol = 5);
   A4<-matrix(0, nrow = 5, ncol = 10); 
   
   #const
   MH2O <- 18.015268;
   Tcr <- 647.067; # {??k}      #critical temperature ofwater from (Johnson and NorT0n, 1991),  Tcr = 647.096??k in http://en. wikipedia.org/ wiki/Critical_point_(thermodynamics)
   Ttr <- 273.16;  # {??k}      #temperature at the triple point
   pcr <- 0.322778; #{g/cm3}   #critical density ofwater, Johnson and NorT0n, 1991
   e0 <- 1.028667;    #Hill, 1990 (p. 1235)
   Dp0 <- 0.23;       #Hill, 1990 (p. 1235)
   DT0 <- 0.05;       #Hill, 1990 (p. 1235)
   O <- 80;      #Hill, 1990 (p. 1236)
   B <- 1;       #Hill, 1990 (p. 1236)
   Y <- 130;     #Hill, 1990 (p. 1236)
   e <- 12;      #Hill, 1990 (p. 1236)
   v <- 4;       #Hill, 1990 (p. 1236) 
   C <- c(8);
   C<- list(c(7.07501275112,-8.34240569963,-0.364601380,-0.036897043,0.003033815,0.000390109,0.113592870,2.413178500));    #Hill, 1990 (Appendix A, p. 1271)
   
   C001 = -0.000034631815;      #Hill, 1990 (p. 1235)
   C002 = -0.000030378112;      #Hill, 1990 (p. 1235)
   GH2Otr = -56290;  # {cal/mol}     #Johnson and NorT0n, 1991
   HH2Otr = -68767;  # {cal/mol}     #Johnson and NorT0n, 1991
   SH2Otr = 15.132;   #{cal/mol}     #Johnson and NorT0n, 1991
   #begin
   p <- p  /pcr;         # {g/cm3}       #Hill, 1990 (p. 1235, Appendix C, p. 1272)
   
   T<- -Tcr / Tk;                      #Hill, 1990 (Appendix C, p. 1272)
   DT<-1+T;                         #Hill, 1990 (p. 1235, Appendix C, p. 1272)
   Dp<-p-1;                         #Hill, 1990 (Appendix C, p. 1272)
   Dkn<-C001+C002*DT;               #Hill, 1990 (p. 1235), Han, 2008 (p. 360, 402)
   DkT<-C002;                       #Han, 2008 (p. 360, 402)
   E<-exp(-p*p);                                      #Hill, 1990 (p. 1236)
   Ep<- -2*p*E;                                        #Hill, 1990 (Appendix E, p. 1273)
   if (p<0.2E-4)
   {
     UnitMinusE<- p*p;
   }
   
   if (!(p<0.2E-4))              #Hill, 1990 (p. 1236); Han, 2008 (p. 361, 402)
   { 
     UnitMinusE<- 1-E;                             #Hill, 1990 (p. 1236)
     tmp<- (-O*DT-B*Dp-Y*DT*DT-e*Dp*Dp);                   #Hill, 1990 (p. 1236, Appendix E, p. 1274); Han, 2008 (p. 361, 402)
     if (tmp>88) 
     {G<-1.01E+38;}                        #Han, 2008 (p. 361, 402)
     if (tmp==88) 
     {G<--1.01E-38;}                      #Han, 2008 (p. 361, 402)
     if ((tmp>=-88) && (tmp<=88)) 
     {
       G<-exp(tmp);      #Hill, 1990 (p. 1236, Appendix E, p. 1274)
       Gp<-G*(-B-2*e*Dp);                                 #Hill, 1990 (Appendix E, p. 1274)
       H<-exp(-v*(T+3));                                  #Hill, 1990 (p. 1236, Appendix E, p. 1274)
       GT<-G*(-O-2*Y*DT);                                 #Hill, 1990 (Appendix E, p. 1274)
       HT<--V*H;                                          #Hill, 1990 (Appendix E, p. 1274)
       h<-sqrt((Dp/Dp0)*(Dp/Dp0)+(DT/DT0)*(DT/DT0));      #Hill, 1990 (p. 1235)
       tmp<- (h/e0^4);
     }
     #Hill, 1990 (p. 1235); Han, 2008 (p. 360, 402)
     if (tmp>88) 
     {
       Z<-1E+38;
     }
     if(tmp<0.01)                        #Han, 2008 (p. 360, 402) 
     {
       Z<-tmp;
     }                               #Han, 2008 (p. 360, 402)
     if(!(tmp<0.01) )
     {
       Z<-exp(tmp)-1;                              #Han, 2008 (p. 360, 402)
       tmp<--1/Z;
     }                                  #Hill, 1990 (p. 1235); Han, 2008 (p. 360, 402)
     if (tmp==88 )
     {
       F<-1 ;
     }                           #Han, 2008 (p. 360, 402)
     if (!(tmp==88) )
     {
       F<-1-exp(tmp);
     }                                         #Hill, 1990 (p. 1235, Appendix E, p. 1273); Han, 2008 (p. 360, 402)
     if (F==0) 
     {
       FT<- 0; 
       Fp<-0;
     }
     if (!(F==0)) 
     {
       
       Fh<-as.complex(-4/e0*exp(-1/Z)*(1/Z+1/Z/Z)*(log10(1+Z)^(3/4)));         #Hill, 1990 (Appendix E, p. 1274)  
       
       Fp<-as.complex(Fh*Dp/h/Dp0/Dp0);        #as complex T0 prevent posible errors    #Hill, 1990 (Appendix E, p. 1274)
       FT<-as.complex(Fh*DT/h/DT0/DT0);                                  
     }              #Han, 2008 (p. 360, 402)
     #Hill, 1990 (Appendix E, p. 1274)
   }
   
   k0<-0;
   k0T<-0;
   for (i in 1:length(6)) 
     
   {
     a<-i;
     
     
     k0<-k0+i*(-T^2-i);               #Hill, 1990 (Appendix A,D, pp. 1271,1273)
     k0T<-k0T-i*(2-i)*(-T^1-i);
   }
   #Hill, 1990 (Appendix D, p. 1273)}
   
   C_Num<-as.numeric(unlist(C)); #Convert the array T0 numeric
   k0<-(C_Num[7]*T+C_Num[8])*log10(-T);               #Hill, 1990 (Appendix A,D, pp. 1271,1273)
   k0T<-(k0T+C_Num[7])*(1+log10(-T))+C_Num[8]/T;           #Hill, 1990 (Appendix D, p. 1273)
   w1<-0;
   w1T<-0;
   w1p<-0;
   
   for (I in 1:length(7)) 
   {  for (J in 1:length(7))
   { 
     A1[I,J]<-0;
   } 
   }
   
   A1[4,1]<-0.3384249125E+00;       #Hill, 1990 (Appendix B, p. 1271)
   A1[5,1]<--0.7153393406E-01;       #Hill, 1990 (Appendix B, p. 1271)
   A1[7,1]<- 0.5493680814E-03;       #Hill, 1990 (Appendix B, p. 1271)
   A1[3,2]<- 0.4933218501E-01;       #Hill, 1990 (Appendix B, p. 1271)
   A1[6,2]<--0.2328491212E-01;       #Hill, 1990 (Appendix B, p. 1271)
   A1[7,2]<- 0.2402095181E-02;       #Hill, 1990 (Appendix B, p. 1271)
   A1[1,3]<- 0.7529422956E+00;       #Hill, 1990 (Appendix B, p. 1271)
   A1[3,3]<--0.2280260070E+01;       #Hill, 1990 (Appendix B, p. 1272)
   A1[2,4]<- 0.1142004144E+01;       #Hill, 1990 (Appendix B, p. 1272)
   A1[3,4]<--0.2619059624E+01;       #Hill, 1990 (Appendix B, p. 1272)
   A1[5,4]<-0.4395237702E+00;       #Hill, 1990 (Appendix B, p. 1272)
   A1[6,4]<--0.3161046646E-01;       #Hill, 1990 (Appendix B, p. 1272)
   A1[7,4]<- 0.6814467692E-03;       #Hill, 1990 (Appendix B, p. 1272)
   A1[1,5]<--0.3924227294E+00;       #Hill, 1990 (Appendix B, p. 1272)
   A1[3,5]<--0.2738770648E+00;       #Hill, 1990 (Appendix B, p. 1272)
   A1[4,6]<--0.1943443857E-01;       #Hill, 1990 (Appendix B, p. 1272)
   A1[5,6]<- 0.3048860434E-02;       #Hill, 1990 (Appendix B, p. 1272)
   A1[3,7]<- 0.3946510403E-02;       #Hill, 1990 (Appendix B, p. 1272)
   
   for (I in 1:length(7)) 
   {
     if (I==2 )
       
     {
       R1<-UnitMinusE*log10(p)-p*p*log10(p)+p*p/2;           #Hill, 1990 (p. 1236; Appendix D, p. 1273)
       R1p<-UnitMinusE/p-Ep*log10(p)-2*p*log10(p);
     }                        #Hill, 1990 (p. 1236; Appendix D, p. 1273)
     else 
     {
       R1<-UnitMinusE* (p^I-2);                                              #Hill, 1990 (p. 1236; Appendix D, p. 1273)
       R1p<-UnitMinusE*(I-2)* (p^I-3)-Ep* (p^I-2);                 #Hill, 1990 (p. 1236; Appendix D, p. 1273)
       for (J in 1:length(7)) 
       {
         w1<- w1+A1[I,J]*R1* (T^J-1);                                           #Hill, 1990 (p. 1236; Appendix D, p. 1273)
         w1T<- w1T+A1[I,J]*R1*(J-1)* (T^J-2);                                   #Hill, 1990 (Appendix D, p. 1273)
         w1p<- w1p+A1[I,J]*R1p* (T^J-1);
         
       }
     }
     #Hill, 1990 (p. 1236; Appendix D, p. 1273)
     
   }
   
   
   w2<-0;
   w2T<-0;
   w2p<-0;
   
   for (I in 1:length(7)) 
   {  for (J in 1:length(12)) 
   {A2[I,J]<-0;}
   }
   A2[1,1]<-  0.2243610314E+00;        #Hill, 1990 (Appendix B, p. 1272)
   A2[2,1]<-  0.1193250201E+00;        #Hill, 1990 (Appendix B, p. 1272)
   A2[5,1]<-  0.6582959348E-01;        #Hill, 1990 (Appendix B, p. 1272)
   A2[5,2]<- 0.1651430628E+00;        #Hill, 1990 (Appendix B, p. 1272)
   A2[1,3]<- -0.2178969357E+01;        #Hill, 1990 (Appendix B, p. 1272)
   A2[2,3]<-  0.2674090542E+00;        #Hill, 1990 (Appendix B, p. 1272)
   A2[3,3]<-  0.8647490995E+00;        #Hill, 1990 (Appendix B, p. 1272)
   A2[1,4]<- -0.1530432257E+00;        #Hill, 1990 (Appendix B, p. 1272)
   A2[3,4]<-  0.2059881454E+01;        #Hill, 1990 (Appendix B, p. 1272)
   A2[6,4]<- -0.4888628703E+00;        #Hill, 1990 (Appendix B, p. 1272)
   A2[7,4]<-  0.1375328753E+00;        #Hill, 1990 (Appendix B, p. 1272)
   A2[5,5]<- -0.9015180666E+00;        #Hill, 1990 (Appendix B, p. 1272)
   A2[6,5]<- -0.1444258609E+00;        #Hill, 1990 (Appendix B, p. 1272)
   A2[7,5]<-  0.1558046279E+00;        #Hill, 1990 (Appendix B, p. 1272)
   A2[4,6]<- -0.2740652563E+01;        #Hill, 1990 (Appendix B, p. 1272)
   A2[6,6]<-  0.4983771706E+00;        #Hill, 1990 (Appendix B, p. 1272)
   A2[4,7]<- -0.3261978564E+01;        #Hill, 1990 (Appendix B, p. 1272)
   A2[5,7]<-  0.1609338784E+01;        #Hill, 1990 (Appendix B, p. 1272)
   A2[1,8]<-  0.3484674963E-01;        #Hill, 1990 (Appendix B, p. 1272)
   A2[2,8]<- -0.1537646434E+01;        #Hill, 1990 (Appendix B, p. 1272)
   A2[5,8]<-  0.2316225257E+00;        #Hill, 1990 (Appendix B, p. 1272)
   A2[2,9]<- -0.1419249232E+01;        #Hill, 1990 (Appendix B, p. 1272)
   A2[3,9]<-  0.7969984635E+00;        #Hill, 1990 (Appendix B, p. 1272)
   A2[5,10]<- 0.7510544627E-02;        #Hill, 1990 (Appendix B, p. 1272)
   A2[1,12]<- 0.5364384732E-03;        #Hill, 1990 (Appendix B, p. 1272)
   
   for (I in 1:length(7)) 
   {
     R2<-  (p^I);                                    #Hill, 1990 (p. 1236; Appendix D, p. 1273)
     R2p<-I* (p^I-1);                                #Hill, 1990 (p. 1236; Appendix D, p. 1273)
     for (J in 1:length(12)) 
     { w2<- w2+A2[I,J]*R2* (T^J-1);             #Hill, 1990 (p. 1236; Appendix D, p. 1273)
     w2T<- w2T+A2[I,J]*R2*(J-1)* (T^J-2);            #Hill, 1990 (Appendix D, p. 1273)
     w2p<- w2p+A2[I,J]*R2p* (T^J-1);                #Hill, 1990 (p. 1236; Appendix D, p. 1273)
     }                 
     
   }
   
   
   
   w3<-0;
   w3T<-0;
   w3p<-0;
   
   for (I in 1:length(5)) 
   { for (J in 1:length(5))
   { A3[I,J]<-0;}
   }
   
   A3[1,1]<- 0.6109381296E+00;                #Hill, 1990 (Appendix B, p. 1272)
   A3[3,1]<--0.1906644459E-01;                #Hill, 1990 (Appendix B, p. 1272)
   A3[5,1]<- 0.7976092188E-02;                #Hill, 1990 (Appendix B, p. 1272)
   A3[1,2]<- 0.1934466766E+01;                #Hill, 1990 (Appendix B, p. 1272)
   A3[1,3]<- 0.1921820547E+01;                #Hill, 1990 (Appendix B, p. 1272)
   A3[3,3]<--0.4410105919E-01;                #Hill, 1990 (Appendix B, p. 1272)
   A3[1,4]<- 0.6130354419E+00;                #Hill, 1990 (Appendix B, p. 1272)
   A3[2,4]<--0.2855258689E+00;                #Hill, 1990 (Appendix B, p. 1272)
   A3[5,4]<- 0.2526137080E-01;                #Hill, 1990 (Appendix B, p. 1272)
   A3[2,5]<--0.2374074642E+00;                #Hill, 1990 (Appendix B, p. 1272)
   A3[4,5]<- 0.3855866402E-01;                #Hill, 1990 (Appendix B, p. 1272)
   A3[5,5]<- 0.8041672150E-02;                 #Hill, 1990 (Appendix B, p. 1272)
   
   
   
   for (I in 1:length(5)) 
   {R3<- (p^I+1);                             #Hill, 1990 (p. 1236; Appendix D, p. 1273)
   R3p<-(I+1)* (p^I);                         #Hill, 1990 (p. 1236; Appendix D, p. 1273)
   for (J in 1:length(10)) 
   { 
     w3<- w3+A3[I,J]*R3* (T^J+1);               #Hill, 1990 (p. 1236; Appendix D, p. 1273)
     w3T<- w3T+A3[I,J]*R3*(J+1)* (T^J);        #Hill, 1990 (Appendix D, p. 1273)
     w3p<- w3p+A3[I,J]*R3p* (T^J+1);            #Hill, 1990 (p. 1236; Appendix D, p. 1273)
   }
   
   }
   
   
   
   
   w4<-0;
   w4T<-0;
   w4p<-0;
   for (I in 1:length(5))
   {  for (J in 1:length(10)) 
   {A4[I,J]<-0;} 
   }
   
   A4[3,1]<- -0.1635439033E+02;                  #Hill, 1990 (Appendix B, p. 1272)
   A4[1,2]<- -0.5025818675E+02;                  #Hill, 1990 (Appendix B, p. 1272)
   A4[2,4]<-  0.1649003040E+00;                  #Hill, 1990 (Appendix B, p. 1272)
   A4[1,5]<- -0.8499893502E+00;                  #Hill, 1990 (Appendix B, p. 1272)
   A4[1,9]<-  0.8314382544E-02;                #Hill, 1990 (Appendix B, p. 1272)
   A4[2,9]<-  0.8781327858E-03;                  #Hill, 1990 (Appendix B, p. 1272)
   A4[2,10]<- 0.1537391213E-02;                  #Hill, 1990 (Appendix B, p. 1272)
   A4[3,10]<--0.9016873786E-03;                  #Hill, 1990 (Appendix B, p. 1272)
   A4[5,10]<- 0.3326628664E-03;                  #Hill, 1990 (Appendix B, p. 1272)
   
   for (I in 1:length(5)) 
   {
     if (I == 2)
     {R1<-UnitMinusE*log10(p)-p*p*log10(p)+p*p/2;                            #Hill, 1990 (p. 1236; Appendix D, p. 1273)
     R1p<-UnitMinusE/p+2*p*exp(-p*p)*log10(p)-2*p*log10(p);                 #Hill, 1990 (p. 1236; Appendix D, p. 1273)
     } else {
       R1<-UnitMinusE* (p^I-2);                                     #Hill, 1990 (p. 1236; Appendix D, p. 1273)
       R1p<-UnitMinusE*(I-2)* (p^I-3)+2*p*exp(-p*p)* (p^I-2);
     } 
     #Hill, 1990 (p. 1236; Appendix D, p. 1273)
     
     
     for (J in 1:length(10)) 
     {   w4<- w4+A4[I,J]*R1* (T^J-1);                                  #Hill, 1990 (p. 1236; Appendix D, p. 1273)
     w4T<- w4T+A4[I,J]*R1*(J-1)* (T^J-2);                          #Hill, 1990 (Appendix D, p. 1273)
     w4p<- w4p+A4[I,J]*R1p* (T^J-1);                               #Hill, 1990 (p. 1236; Appendix D, p. 1273)
     }
     
     
   }
   
   
   
   k1<- w1+E* w2+G* w3+H* w4;                                       #the "far-field" Helmholtz function, Hill, 1990 (p. 1236, Appendix D, p. 1272)
   
   kf<-log10(p)+(k0+k1+F*Dkn);         #adicional () T0 prevent errors when k1 is inf                           #Helmholtz function, Hill, 1990 (p. 1235, Appendix D, p. 1272)
   
   k<-kf*Tk*Rcal;    #  {cal/mol}   #Helmholtz free energy ofwater, Hill, 1990 (p. 1235)
   
   if(kf==Inf)  #case kf is Inf
   {k<-Inf;}
   if(kf==-Inf)
   {k<--Inf;}
   
   k0p<-1/p+ w1p+Ep* w2+E* w2p+Gp* w3+G* w3p+H* w4p+Dkn*Fp;              #Hill, 1990 (Appendix D, p. 1272)
   if(Gp==Inf||G==Inf)  #case Gp or G is Inf
   {k0p<-Inf;}
   if(Gp==-Inf||G==-Inf)  #case Gp or G is Inf
   {k0p<--Inf;}
   kT<-k0T+ w1T+E* w2T+GT* w3+G* w3T+HT* w4+H* w4T+FT*Dkn+F*DkT;         #Hill, 1990 (Appendix D, p. 1273)
   
   HH2O<-k+Rcal*Tcr*(-kT-p/T*k0p)+HH2Otr;   #{cal/mol}                #Enthalpy ofwater, Johnson and NorT0n, 1991 (p. 586)
   
   SH2O<-Rcal*(T*kT)+SH2Otr;     #{cal/??k/mol}                          #Entropy ofwater, Johnson and NorT0n, 1991 (p. 586)
   
   
   Result<-HH2O-Tk*SH2O+GH2Otr-HH2Otr+Ttr*SH2Otr;# {cal/mol}        #Gibbs free energy ofwater, Helgeson and kirkham, 1974 (p. 1096)
   
 }
 #okey
 Born<- function( w, charge, Tk, Pbar, gH2O,  eTP,T0)     #//Calculates the Gibbs free energies of solvation
 {
   
   wTP<-0;
   re<-0;
   
   Y = -5.8495398E-5;     #{1/??k}                //Shock et al., 1992 (Appendix D, p. 824)
   n = 166027;            #{??*cal/mol}           // n=NAv*e*e/2 (e=4.80298E-10 {cm^1.5*g^0.5/s} - electronic charge, Shock et al., 1992 (p. 803, Appendix A, p. 818)
   e = 78.38092765;       #{}                    //Dielectric constant ofwater at T = 298.15??k and P = 1 bar
   kz = 0.94;             #{??}                   //Constant for cations, Shock et al., 1992 (p. 805)
   
   if (!(charge==0))                         #//Born coefficient of the ion, Shock et al., 1992 (p. 803-805, Appendix D, p. 824)
   { 
     re<-charge*charge/( w*100000/ n+charge/3.082);   #//Effective electroctatic radii of ions, 3.082 corresponds T0 the effective electrostatic radius of H+ at 1 bar and 298.15??k}
     
     if (charge>0)                          #//Cations
     {
       re<-re-charge*kz;                             #//Crystallographic radii of cations, Shock et al., 1992 (p. 804)
       wTP<- n*charge*(charge/(re+charge*(kz+gH2O))-1/(3.082+gH2O));                #//Shock et al., 1992 (Eq. D4, p. 824) 
       
     }
     else                                                                     #//Anions
     {  
       wTP<- n*charge*(charge/(re-charge*gH2O)-1/(3.082+gH2O));                     #//Shock et al., 1992 (Eq. D4, p. 824)
       
     } 
     
     #//Neutral molecules
     
   }
   else
   {
     wTP<- w*100000; 
   }  
   Result<- wTP*(1/ eTP-1)+ w*100000*(1-1/ e)+ w*100000*Y*(Tk-T0);                    #//Shock et al., 1992 (Appendix B, Eq. B25, p. 821)} )
   CEP_2$Gibbs[row]<-Result;
 }
 AqGibbsFreeEnergy<- function(Aq,Tk, Pbar, gH2O,  eTP,CEP_2,row)     #//Calculates Gibbs Free energies of aqueous species
 {
   
   G<-0;
   S<-0; 
   a1<-0;
   a2<-0;
   a3<-0;
   a4<-0;
   c1<-0;
   c2<-0;
   w<-0; 
   charge<-0;
   theta = 228; #{??k}   //Solvent constant, Shock et al., 1992 (Appendix A-B, pp. 819-820)
   psi = 2600;  #{bar}  //Solvent constant, Shock et al., 1992 (Appendix A-B, pp. 819-820)
   
   
   

   
   

  
    G<- as.numeric(CEP_2$G..cal.mol.1[[row]]) ;

    S<- as.numeric(CEP_2$S..cal.mol.1.?.K.1[[row]]) ;
      
  

     a1<- as.numeric(CEP_2$a..cal.mol.1.?.K.1[[row]]) ;
  
     a2<- as.numeric(CEP_2$SD..cal.mol.1.?.K.1[row])  ;
  
    a3<- as.numeric(CEP_2$c?.10.5..cal.mol.1.?.K.1[row]) ;
  
    a4<- as.numeric(CEP_2$X[row]);
   
    c1<- as.numeric(CEP_2$X.1[row]) ;   
    
    c2<- as.numeric(CEP_2$X.3[row]) ;
    
    w<- as.numeric(CEP_2$X.5[row]);
    
     if ( !(is.na( CEP_2$X.7[row])))
     {
       charge<- as.numeric(CEP_2$X.7[row]) ;
     } 
     else 
     { 
       charge<-10;
     }
   d<-log10((psi+Pbar)/(psi+1))-c2*10000*(1/(Tk-theta)-1/(T0-theta));#*(theta-Tk)/theta-Tk/theta/theta*log10(T0*(Tk-theta)/Tk/70.15))+(a3*(Pbar-1)+a4*10000*log10((psi+Pbar)/(psi+1)));#/(Tk-theta)#+Born( w,charge,Tk,Pbar,gH2O, eTP);         #//Tanger and Helgeson, 1988; Shock et al., 1989; Shock et al., 1992; Johnson et al., 1992; Sverjensky et al., 1997
   Result<- G -S*(Tk-T0)+c1*(Tk-T0-Tk*log10(Tk/T0))+a1/10*(Pbar-1)+a2*100;  

     
   }

   
 
 #-------------------------------------------------
 AqGibbsFreeEnergyUncertainty<- function(Aq, Tk, Pbar,CEP_2,row)    #//Calculates uncertainties of aqueous species
 {  
   oG<-0;
   oS<-0;
   oa1<-0;
   oa2<-0;
   oa3<-0;
   oa4<-0; 
   oc1<-0; 
   oc2<-0; 
   
   theta = 228; #{??k}   //Solvent constant, Shock et al., 1992 (Appendix A-B, pp. 819-820)
   psi = 2600;  #{bar}  //Solvent constant, Shock et al., 1992 (Appendix A-B, pp. 819-820)
   


   
     
       oG<- as.numeric(CEP_2$G..cal.mol.1[row]) ;
     
       oS<- as.numeric(CEP_2$S..cal.mol.1.?.K.1[row])  ;
     
  
       oa1<- as.numeric(CEP_2$a..cal.mol.1.?.K.1[row] );
  
       oa2<- as.numeric(CEP_2$b?.1000..cal.mol.1.?.K.1[row] ) ;
    
  
       oa3<- as.numeric(CEP_2$c?.10.5..cal.mol.1.?.K.1[row]) ;
    
  #  oa4<- as.numeric(CEP_2$X[row]) ;  
       
    
  #  oc1<- as.numeric(CEP_2$X.2[row]) ;
   
     # oc2<- as.numeric(CEP_2$X.4[row]) ; 

     Result<-sqrt((oG^2)+ (oS*(Tk-T0)^2)+ (oc1*(Tk-T0-Tk*log10(Tk/T0))^2)+ (oa1/10*(Pbar-1)^2)+ (oa2*100*log10((psi+Pbar)/(psi+1))^2));    #//+ (o w*100000*((1-1/epsilon)+Y*(Tk-T0)),2));  Result 
 }    
   
                            
 GasGibbsFreeEnergy<-function(Gas,Tk,Pbar,CEP_2,row)                #//Calculates Gibbs Free energies of gas species
 {
   

   
   G<-0;
   S<-0;
   a<-0;
   b<-0
   c<-0;
   
  
   #View(CEP) #Visualizaci??n de la base
   
   
   
   
   G<- as.numeric(CEP_2$G..cal.mol.1[row]) ;
   S<- as.numeric(CEP_2$S..cal.mol.1.?.K.1[row]) ;
   a<- as.numeric(CEP_2 $a..cal.mol.1.?.K.1[row] ) ;
   b<- as.numeric(CEP_2 $b?.1000..cal.mol.1.?.K.1 [row]) ;
   c<- as.numeric(CEP_2 $c?.10.5..cal.mol.1.?.K.1[row]) ;
   
   
   
   
   
   
   
   Result<-G-S*(Tk-T0)+a*(Tk-T0-Tk*log10(Tk/T0))-(c*100000+b/1000*T0*T0*Tk)*(Tk-T0)*(Tk-T0)/2/T0/T0/Tk;         #//Helgeson et al., 1978
   

 }
 
 GasGibbsFreeEnergyUncertainty<-function(Gas, Tk,CEP_2,row)       #//Calculates uncertainties of gas species
 {            
   
   oG<-0;
   oS<-0;
   oa<-0;
   ob<-0;
   oc<-0;
   T0<-2;
   
   library(dplyr) #Cargar paquete, si no est?? cargado desde antes.
   

   
   # as.numeric(SD..cal.mol.1);#[4]
   #as.numeric(CEP$SD..cal.mol.1.?.K.1);#[6]
   #as.numeric(CEP$SD..cal.mol.1.?.K.1.1);#[10]
   #as.numeric(CEP$SD.1000..cal.mol.1..K.1);#[12]
   #as.numeric(CEP$SD.10.5..cal.mol.1..K.1);#[14]
  
       oG<-as.numeric(CEP_2$G..cal.mol.1[row]);#[4];
    
       oS<- as.numeric(CEP_2$S..cal.mol.1.?.K.1[row]);#[6]
     
    
       oa<- as.numeric(CEP_2 $a..cal.mol.1.?.K.1[row] );
    
    
       ob<- as.numeric(CEP_2 $b?.1000..cal.mol.1.?.K.1[row] ) ;
    
       oc<- as.numeric(CEP_2 $c?.10.5..cal.mol.1.?.K.1[row]); 
    
     
     Result<-sqrt( (oG^2)   + (oS*(Tk-T0)^2) + (oa*(Tk-T0-Tk*log10(Tk/T0))^2)+ (-ob/1000*(Tk-T0)*(Tk-T0)/2^2)+ (oc*100000*(Tk-T0)*(Tk-T0)/2/T0/T0/Tk^2));
     
     
   
 }


end

