#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

library(dplyr) #Cargar paquete, si no está cargado desde antes.
CEP_csv  <- read.csv2("C:\\Users\\milax\\Documents\\!Crono_November_2020\\Data base\\ActualDatabase.csv")


waterGibbsFreeEnergy<- function(Tk,Pbar,p) numeric();
Born<- function( w, charge, Tk, Pbar, gH2O,  eTP)numeric();
SolvFunc <- function(p, TC, Pbar)numeric(); 
waterDensity<- function(Tk,Pbar,p) numeric();
GasGibbsFreeEnergy<- function(Gas,Tk,p) numeric();
GasGibbsFreeEnergyUncertainty<- function(Gas, Tk,Pbar)
AqGibbsFreeEnergyUncertainty<- function(Aq, Tk, Pbar)
AqGibbsFreeEnergy<- function(Aq,Tk, Pbar, gH2O,  eTP)
MinGibbsFreeEnergyUncertainty<- function(Min, Tk) 
MinGibbsFreeEnergy<- function(Min,Tk,T0,Pbar)
DielConst<- function(p, TC, Tk, Pbar) 

  
MinGibbsFreeEnergyUncertainty <- function(Min, Tk)    #Calculates uncertainties of minerals
{
  oG<-0;
  oS<-0; 
  oV<-0; 
  oa<-0; 
  ob<-0;
  CEP <-filter(CEP_csv, Name == Min);

  

  if (!(  CEP$Name == Min))
  {
    Result<-0     #Find existence of the mineral in the Data base
  }
  else 
  {
    if ( !(is.na( CEP$SD..cal.mol.1)))
  {oG<-CEP$SD..cal.mol.1}else {oG<-0;} 
    if ( !(is.na(CEP$S..cal.mol.1..K.1)) )
   {oS<-CEP$SD..cal.mol.1..K.1} else {oS<-0;}
    if( !(is.na(CEP$SD..cm3.mol.1)))
   {oV<-CEP$SD..cm3.mol.1} else {oV<-0;}
    if (!(is.na(CEP$SD..cal.mol.1..K.1.1)))
    {oa<-CEP$SD..cal.mol.1..K.1.1} else {oa<-0;}
    if (!(is.na(CEP$SD.1000..cal.mol.1..K.1)))
    {ob<-CEP$SD.1000..cal.mol.1..K.1} else {ob<-0;}
    if (!(is.na(CEP$SD.1000..cal.mol.1..K.1)))
    {oc<-CEP$SD.10.5..cal.mol.1..K.1} else {oc<-0;}
    Result<-sqrt( (oG^2)+ (oS*(Tk-T0)^2)+ (oa*(Tk-T0-Tk*log(Tk/T0))^2)+ (-ob/1000*(Tk-T0)*(Tk-T0)/2^2)+ (oc*100000*(Tk-T0)*(Tk-T0)/2/T0/T0/Tk^2)+ (oV*(Pbar-1)/41.84^2)); 
    
  }
  
  
}

 waterDensity <- function(Tk, Pbar)          #Calculates a density ofwater in g/cm3
{
   #from CronoMF
   dyn.load ('okawsp6.dll');
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

SolvFunc <- function(p, TC, Pbar)           #Solvent function (g), Shock et al , 1992; Sue et al , 2002 - valid at ρ>=0.20 and <=1, ρ>=0.35 - Shock et al., 1992 (p. 809, Fig. 6)
{
g<-0;
ag<-0;
bg<-0;
p_a<-0; 
ag_i = -2.037662;         #Å    Shock et al., 1992 (Table 3, p. 807)
ag_ii = 0.005747;         #Å °k     Shock et al., 1992 (Table 3, p. 807)
ag_iii = -0.000006557892; #Å °k^2  Shock et al., 1992 (Table 3, p. 807)
bg_i = 6.107361;          #         Shock et al., 1992 (Table 3, p. 807)
bg_ii = -0.01074377;      #1 °C    Shock et al., 1992 (Table 3, p. 807)
bg_iii = 0.00001268348;   #1 °C^2  Shock et al., 1992 (Table 3, p. 807)
ag_1 = 3.66666E-16;       #°k       Shock et al., 1992 (Table 4, p. 809)
ag_2 = -1.504956E-10;     #Å bar^3  Shock et al., 1992 (Table 4, p. 809)
ag_3 = 5.01799E-14;       #Å bar^4  Shock et al., 1992 (Table 4, p. 809)
cg_1 = 0.18359;           #        Sue et al., 2002 (p. 3304)
cg_2 = -0.18632;          #         Sue et al., 2002 (p. 3304)
cg_3 = 0.11531;           #         Sue et al., 2002 (p. 3304)

p_a<-p/ waterDensity(T0,1);
if (p_a>1 )
{ 
p_a<-1;                                         #//Need to avoid the calculation collapse at T>25°C and provide some minor error on the Born function
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
  g<-g-(cg_1*p_a+cg_2* (p_a^2)+cg_3*log(p_a));             #Sue et al., 2002 (Eq. 25, p. 3304) - valid at T - from 250 to 600°C, Pbar - from 223 to 998 bar, ρ - from 0.20 to 0.81 g cm-3 (p. 3303)
}
 
   
Result<-g;   
 
}  
  
DielConst <- function(p, TC, Tk, Pbar)                      #//Dielectric constant ofwater ( e) at specified T and P
{
  y<-0;
  c<-0;
  s<-0;
  d<-0;
  
  a_1 <- -1.576377E-3;    #1/°C                //Sverjensky et al., 2014 (Table 1, p. 131)
  a_2 <- 6.810288E-2;     #1/sqrt(°C)         //Sverjensky et al., 2014 (Table 1, p. 131)
  a_3 <- 7.548755E-1;     #                    //Sverjensky et al., 2014 (Table 1, p. 131)
  b_1 <- -8.016651E-5;    #1/°C                //Sverjensky et al., 2014 (Table 1, p. 131)
  b_2 <- -6.871618E-2;   #1/sqrt(°C)          //Sverjensky et al., 2014 (Table 1, p. 131)
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
  



MinGibbsFreeEnergy <- function(Min,Tk,T0,Pbar)  
{ 

 
  CEP_csv  <- read.csv2("C:\\Users\\milax\\Documents\\!Crono_November_2020\\Data base\\ActualDatabase.csv")
  library(dplyr) #Cargar paquete, si no está cargado desde antes.

CEP <-filter(CEP_csv, Name == Min)
View(CEP) #Visualización de la base

G<-as.numeric(CEP$G..cal.mol.1);
S<- as.numeric(CEP$SD..cal.mol.1);


V<- as.numeric(CEP$S..cal.mol.1..K.1);

a<- as.numeric(CEP$SD..cal.mol.1..K.1);
b<- as.double(CEP$V..cm3.mol.1); 
c<- as.numeric(CEP$SD..cm3.mol.1);

A<-10;
B<-40;


Result <- G-S*(Tk-T0)+A*(Tk-T0-Tk*log10(Tk/T0))-(c*100000+B/1000*T0*T0*Tk)*(Tk-T0)*(Tk-T0)/2/T0/T0/Tk+V*(Pbar-1)/41.84;  
} 

#other functions
 waterGibbsFreeEnergy<- function(Tk,Pbar,p) #Calculates the Gibbs energy ofwater from triple point Ptr and Ttr to 10 kbar and 1273.15 °k (Helgeson and kirkham, 1974; Hill, 1990; Johnson and Norton, 1991; Han, 2008)
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
Tcr <- 647.067; # {°k}      #critical temperature ofwater from (Johnson and Norton, 1991),  Tcr = 647.096°k in http://en. wikipedia.org/ wiki/Critical_point_(thermodynamics)
Ttr <- 273.16;  # {°k}      #temperature at the triple point
pcr <- 0.322778; #{g/cm3}   #critical density ofwater, Johnson and Norton, 1991
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
GH2Otr = -56290;  # {cal/mol}     #Johnson and Norton, 1991
HH2Otr = -68767;  # {cal/mol}     #Johnson and Norton, 1991
SH2Otr = 15.132;   #{cal/mol}     #Johnson and Norton, 1991
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

    Fp<-as.complex(Fh*Dp/h/Dp0/Dp0);        #as complex to prevent posible errors    #Hill, 1990 (Appendix E, p. 1274)
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

C_Num<-as.numeric(unlist(C)); #Convert the array to numeric
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

kf<-log10(p)+(k0+k1+F*Dkn);         #adicional () to prevent errors when k1 is inf                           #Helmholtz function, Hill, 1990 (p. 1235, Appendix D, p. 1272)

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

HH2O<-k+Rcal*Tcr*(-kT-p/T*k0p)+HH2Otr;   #{cal/mol}                #Enthalpy ofwater, Johnson and Norton, 1991 (p. 586)

SH2O<-Rcal*(T*kT)+SH2Otr;     #{cal/°k/mol}                          #Entropy ofwater, Johnson and Norton, 1991 (p. 586)


Result<-HH2O-Tk*SH2O+GH2Otr-HH2Otr+Ttr*SH2Otr;# {cal/mol}        #Gibbs free energy ofwater, Helgeson and kirkham, 1974 (p. 1096)

}
#okey
Born<- function( w, charge, Tk, Pbar, gH2O,  eTP,T0)     #//Calculates the Gibbs free energies of solvation
{
 
   wTP<-0;
  re<-0;
  
  Y = -5.8495398E-5;     #{1/°k}                //Shock et al., 1992 (Appendix D, p. 824)
  n = 166027;            #{Å*cal/mol}           // n=NAv*e*e/2 (e=4.80298E-10 {cm^1.5*g^0.5/s} - electronic charge, Shock et al., 1992 (p. 803, Appendix A, p. 818)
  e = 78.38092765;       #{}                    //Dielectric constant ofwater at T = 298.15°k and P = 1 bar
  kz = 0.94;             #{Å}                   //Constant for cations, Shock et al., 1992 (p. 805)
  
  if (!(charge==0))                         #//Born coefficient of the ion, Shock et al., 1992 (p. 803-805, Appendix D, p. 824)
  { 
    re<-charge*charge/( w*100000/ n+charge/3.082);   #//Effective electroctatic radii of ions, 3.082 corresponds to the effective electrostatic radius of H+ at 1 bar and 298.15°k}
  
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
}
AqGibbsFreeEnergy<- function(Aq,Tk, Pbar, gH2O,  eTP)     #//Calculates Gibbs Free energies of aqueous species
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
  theta = 228; #{°k}   //Solvent constant, Shock et al., 1992 (Appendix A-B, pp. 819-820)
  psi = 2600;  #{bar}  //Solvent constant, Shock et al., 1992 (Appendix A-B, pp. 819-820)
  

  
  CEP <-filter(CEP_csv, Name == Aq)


  if (! ( CEP$Name == Aq))
  {
    printf(Aq+' is not found in the data base!'); 
    Result<-0; 
  }                          #//Looks for the aqueous species in the data base
  else 
  {
    if ( !(is.na( CEP[4])))
    {  G<- as.numeric(CEP[4]) ;}
    else
    { G<-0;}
    if( !(is.na( CEP[6])) )
     { S<- as.numeric(CEP[6]) ;}
    else 
      {S<-0;}
    if ( !(is.na(  CEP[10]))) 
      {a1<- as.numeric(CEP[10]) ;}
    else 
      {a1<-0;}
    if ( !(is.na( CEP[12]))) 
      {a2<- as.numeric(CEP[12])  ;}
    else
      {a2<-0;}
    if ( !(is.na( CEP[14])) )
    {
      a3<- as.numeric(CEP[14]) ;
    }
    else
    {
      a3<-0;
    }
    if ( !(is.na( CEP[16])))
    {
      a4<- as.numeric(CEP[16]);
    }
    else 
    {
      a4<-0;
    }
    if ( !(is.na( CEP[18])))
     {
      c1<- as.numeric(CEP[18]) ;
    }
    else 
    {
       c1<-0;
      }
    if ( !(is.na( CEP[20]))) 
    {
      c2<- as.numeric(CEP[20]) ;
    }
    else
    {  
      c2<-0;
    }
    if ( !(is.na( CEP[22]))) 
    {
     w<- as.numeric(CEP[22]);
    }
    else
    { 
      w<-0;
    }
    if ( !(is.na( CEP[24])))
    {
      charge<- as.numeric(CEP[24]) ;
    } 
    else 
    { 
      charge<-10;
    }
     

  }
  d<-log10((psi+Pbar)/(psi+1))-c2*10000*(1/(Tk-theta)-1/(T0-theta));#*(theta-Tk)/theta-Tk/theta/theta*log10(T0*(Tk-theta)/Tk/70.15))+(a3*(Pbar-1)+a4*10000*log10((psi+Pbar)/(psi+1)));#/(Tk-theta)#+Born( w,charge,Tk,Pbar,gH2O, eTP);         #//Tanger and Helgeson, 1988; Shock et al., 1989; Shock et al., 1992; Johnson et al., 1992; Sverjensky et al., 1997
  Result<- G -S*(Tk-T0)+c1*(Tk-T0-Tk*log10(Tk/T0))+a1/10*(Pbar-1)+a2*100;
 
}
#-------------------------------------------------
AqGibbsFreeEnergyUncertainty<- function(Aq, Tk, Pbar)    #//Calculates uncertainties of aqueous species
{
oG<-0;
oS<-0;
oa1<-0;
oa2<-0;
oa3<-0;
oa4<-0; 
oc1<-0; 
oc2<-0; 

theta = 228; #{°k}   //Solvent constant, Shock et al., 1992 (Appendix A-B, pp. 819-820)
psi = 2600;  #{bar}  //Solvent constant, Shock et al., 1992 (Appendix A-B, pp. 819-820)

CEP <-filter(CEP_csv, Name == Aq)
View(CEP) #Visualización de la base

if(! ( CEP$Name == Aq))   #//Find existence of the aqueous species in the data base
{
  Result<-0;
}
else
 {
  if (!(is.na( CEP[4])))
  {
  oG<- as.numeric(CEP[4]) ;
  } 
   else 
   {
    oG<-0;
    }
  
  if (!(is.na( CEP[6])))
  {
    oS<- as.numeric(CEP[6])  ;
  }
  else
  {
    oS<-0;
  }

  if (!(is.na( CEP[10]) ))
   { 
    oa1<- as.numeric(CEP[10]);
   }
  else 
  {
    oa1<-0;
  }
  if(!(is.na(CEP[12])) )  
  {
    oa2<- as.numeric(CEP[12]) ;
    } 
  else
    {
    oa2<-0;
    }
  if (!(is.na(CEP[14])) )
  {
    oa3<- as.numeric(CEP[14]) ;
   }
  else
  { 
  oa3<-0;
    }
   
  if (!(is.na(CEP[16])) )
 {
    oa4<- as.numeric(CEP[16]) ;  
    
  } 
  else 
  {
      oa4<-0;
  }
  if(!(is.na(CEP[18])) ) 
  {
    oc1<- as.numeric(CEP[18]) ;
  } 
    else 
    {
    oc1<-0;
    }
if (!(is.na(CEP[20])) )
  {
    oc2<- as.numeric(CEP[20]) ;
    }
  else
  {
      oc2<-0;
  }
Result<-sqrt( (oG^2)+ (oS*(Tk-T0)^2)+ (oc1*(Tk-T0-Tk*log10(Tk/T0))^2)+ (oa1/10*(Pbar-1)^2)+ (oa2*100*log10((psi+Pbar)/(psi+1))^2)+ (oc2*10000*((1/(Tk-theta)-1/(T0-theta))*(theta-Tk)/theta-Tk/theta/theta*log10(T0*(Tk-theta)/Tk/70.15))^2)+ (oa3*(Pbar-1)/(Tk-theta)^2)+ (oa4*10000*log10((psi+Pbar)/(psi+1))/(Tk-theta)^2));    #//+ (o w*100000*((1-1/epsilon)+Y*(Tk-T0)),2));

 }
}
 GasGibbsFreeEnergy<-function(Gas,Tk)                #//Calculates Gibbs Free energies of gas species
 {


   Tk<-2; 
   G<-0;
   S<-0;
   a<-0;
   b<-0
   c<-0;
   CEP <-filter(CEP_csv, Name == Gas)
   View(CEP) #Visualización de la base
   if (!(CEP$Name == Gas))     #//Find existence of the gas in the data base
   {
     # print(Gas ,' is not found in the data base!');
     Result<-0; 
   } 
   else
    {
       
       if (!( is.na(CEP[4])) ) 
     {
       G<- as.numeric(CEP[4]) ;
     }
     else
     {
       G<-0;
     }
     if (!( is.na(CEP[6]) ) )
     {
       S<- as.numeric(CEP[6]) ;
     }
     else 
     {
       S<-0;
     }
     if (!( is.na(CEP[10]) ) )
     {
       a<- as.numeric(CEP[10]) ;
     } 
     else 
     {
       a<-0;
     }
     if (!( is.na(CEP[12]) ) )
     {
       b<- as.numeric(CEP[12]) ;
     }
     else 
     {
       b<-0;
     }
     if (!( is.na(CEP[14] ) ))
     {
       c<- as.numeric(CEP[14]) ;
     } 
     else
     {
       c<-0;
     }
   
    }
   
    
      
   
   Result<-G-S*(Tk-T0)+a*(Tk-T0-Tk*log10(Tk/T0))-(c*100000+b/1000*T0*T0*Tk)*(Tk-T0)*(Tk-T0)/2/T0/T0/Tk;         #//Helgeson et al., 1978
   end;
 }

GasGibbsFreeEnergyUncertainty<-function(Gas, Tk)       #//Calculates uncertainties of gas species
{

  oG<-0;
  oS<-0;
  oa<-0;
  ob<-0;
  oc<-0;
  T0<-2;
  
  library(dplyr) #Cargar paquete, si no está cargado desde antes.
  CEP_csv  <- read.csv2("C:\\Users\\milax\\Documents\\!Crono_November_2020\\Data base\\ActualDatabase.csv")
  CEP <-filter(CEP_csv, Name == Gas);

 # as.numeric(CEP$SD..cal.mol.1);#[4]
  #as.numeric(CEP$SD..cal.mol.1..K.1);#[6]
  #as.numeric(CEP$SD..cal.mol.1..K.1.1);#[10]
  #as.numeric(CEP$SD.1000..cal.mol.1..K.1);#[12]
  #as.numeric(CEP$SD.10.5..cal.mol.1..K.1);#[14]
  if(is.na( CEP$Name == Gas))
   {
    Result<-0 ;
    }   #//Find existence of the mineral in the Data base
   else 
    {  
    if (!(is.na(  as.numeric(CEP[4])) ))
    {
    oG<-as.numeric(CEP[4]);#[4];
    } 
    else
    {
    oG<-0;
    } 
    if (!(is.na( as.numeric(CEP[6])) ))
      {
      oS<- as.numeric(CEP[6]);#[6]
      }
      else
      {
        oS<-0;
      }
    if (!(is.na(CEP[10]) ))
      {
      oa<- as.numeric(CEP[10]);
       }
      else
      {
        oa<-0;
      }
       
    if (!(is.na(CEP[12]) ))
      {
      ob<- as.numeric(CEP[12]) ;
      } 
      else 
      {
        ob<-0;
      }  
    if (!(is.na(CEP[14]) ))
     {
      oc<- as.numeric(CEP[14]); 
     }
      else
     {
        oc<-0;
     }

    Result<-sqrt( (oG^2)   + (oS*(Tk-T0)^2) + (oa*(Tk-T0-Tk*log10(Tk/T0))^2)+ (-ob/1000*(Tk-T0)*(Tk-T0)/2^2)+ (oc*100000*(Tk-T0)*(Tk-T0)/2/T0/T0/Tk^2));
  
  
     }
}




end

