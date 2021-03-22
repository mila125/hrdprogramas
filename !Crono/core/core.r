
#Math


#{*******************************************************************************************************
#    Core - a module for calculation of equilibrium in the minerals-solution-gas system, realize the Conjugate
#  Directions algorithm (Algorithm#1) and the method combining the Simplex linear minimization with gradient
#                        search modified from the approach of de Capitani-Brown (Algorithm#2)
#Developed by Alexey Novoselov (Institute of Geosciences, University of Campinas (UNICAMP)) in 2013-2015
#Supported by FAPESP grant No. 2011/12682-3 and by CNPq grant No. 164939/2014-8
#Last revision on 24 August 2015
#                                                                                
#LEGEND
#                                                                                
#n   - the number of independent components of the system (elements+charge)
#m   - the number of dependent components (sum of minerals, water, solution species and gases, m=min+aq+gas)
#min - the number of minerals
#aq -  the number of solution species+water
#gas - the number of gases
#H2O - the position of H2O in the list of dependent components
#G0  - initial value of the Gibbs energy at given P-T normalized to R*T, for solution species also +ln(Mw)
#Gbulk - total Gibbs energy of the system
#Gsol - Gibbs energy of each solution (phase)
#vij - stoichiometric coefficients
#b   - balance of the system in mol
#x   - content of the dependent components (solution of the task) in mol
#xmin - the minimum possible meaning of x
#Aphi, Bphi - Debye-Huckel parameters for water
#IoS - ionic strength
#phi - osmotic coefficient of water
#eps - negligibly small value, can change from 1e-6 to 1e-14 (input by user)
#rp  - penalty parameter
#Lag - Lagrange multipliers
#g   - gradient vector
#dk  - directions
#Betta - Fletcher-Reeves conjugate parameter
#MP  - mole portion of the solution (portion of 1) (analogue for x)
#SB  - solution balance (analogue of b for each solution)
#SB_2 - dublicate of SB
#SB_3 - transpose of SB_2
#Title - descriptor of each solution: (i) if Title<0, this is a mineral (i.e. -1 is a first mineral, -2 - second, ...);
#(ii) if Title=0, this is a end-member of aqueous solution, it should be saved during all calculation procedure;
#(iii) if Title>0, this is a calculated solution, it should be deleted from the system, if Title>lim (limiting number of iterations)
#N_Gauss - the maximum possible number of solutions (phases) in the system, calculated with Gauss formula
#Nmax - the current number of solutions (phases) in the system
#Penalty - parameter estimating deviation from imposed restrictions (deviation between calculated and inputed balance)
#D, DNew - optimal decisions obtained with Simplex method
                                                                                
#                                                                                
#REFERENCES
#                                                                               
#de Capitani, C., 1987. PhD Thesis. The computation of chemical equilibrium and the distribution of Fe, Mn and Mg among sites and phases in olivines and garnets. The University of British Columbia, 273 p.
#de Capitani, C., Brown, T.H., 1987. The computation of chemical equilibrium in complex systems containing non-ideal solutions. Geochemica et Cosmochemica Acta 51, 2639?2652.
#Mironenko, M.V., Grant, S.A., Marion, G.M., Farren, R.E., 1997. FREZCHEM2 ? A chemical thermodynamic model for electrolyte solutions at subzero temperatures. CRREL Report 97-5, 44 p.
#Pitzer, K.S., 1987. A thermodynamic model for aqueous solutions of liquid-like density. Reviews in Mineralogy 17, 97?142.
# Spencer, R.J., Moller, N., Weare, J.H., 1990. The prediction of mineral solubility in natural waters: A chemical equilibrium model for the Na-K-Ca-Mg-Cl-SO4-H2O system at temperatures below 25?C. Geochimica et Cosmochimica Acta 54, 575?590.
                                                                                
#********************************************************************************************************}

#const
Mw = 1/0.0180153011158497;#//55.50837;  {mol kg-1 H2O}  //from Mironenko et al., 1997
nmin = 1e-66;

  
ExApp<-0;

 MinimumSearch   <- function(vij, x , G0, b, diam, IoS, eps, Aphi, ADH, BDH, Bdot, m, n, min, aq, gas)  numeric() ;    #//Algorithm#1
MinimumSearch_NR <- function(vij,  x, G0, b, diam, IoS, eps, Aphi, ADH, BDH, Bdot, m, n, min, aq, gas) numeric();
 MinimumSearch_New2 <- function(vij,x, G0, b, IoS, eps, Aphi, Kw, m, n, min, aq, gas, H, OH) numeric();
 MinimumSearch_CB <- function(vij, x, G0, b, diam, IoS, eps, Aphi, ADH, BDH, Bdot, m, n, min, aq, gas)  numeric() ;    #//Algorithm#2
MinimumSearch_Simplex <- function(vij,  x, G0, b, diam, IoS, eps, Aphi, ADH, BDH, Bdot, m, n, min, aq, gas) numeric();
MinimumSearch_Plus <- function(vij,  x, G0, b, diam, IoS, eps, Aphi, ADH, BDH, Bdot, m, n, min, aq, gas) numeric();

IonicStrength <- function(z, x, H2O, minaq)as.numeric()        #//Ionic Strength
  
i<-0;

Result<-0;

for (  i in 1:length(minaq)) 
if z[i]<>0 {}                                         #//if the charge of the solution species is not zero
Result<-Result+x[i]*z[i]*z[i];                     #//Pitzer, 1987, p.6
Result<-0.5*Result/x[H2O]*Mw;//0.0180153011158497);#//x[H2O]*Mw/2;
                        end;
                        
                        function phi(z, x: array of extended; species, Aphi: extended; H2O, minaq: integer): extended;           #//Osmotic coefficient (Pitzer, 1987)
                          
                        Msum, Bphi, IoS, SUMPHI, SCATON, SANON: extended;
                        begin
                        Msum<-species/x[H2O]*Mw;#//0.0180153011158497);#//x[H2O]*Mw;
  IoS<-IonicStrength(z,x,H2O,minaq);           #//Ionic Strength
  Bphi<-1.2;                                      #//a universal parameter from the Debye-Huckel equation, [kg**1/2 mol**-1/2], Pitzer, 1987, p.9
  SUMPHI<-0.0;                                      #//parameters from Mironenko et al., 1997
  SCATON<-0;                                      #//
    SANON<-0;                                      # //
    Result<-1+2/Msum*(-Aphi*power(IoS,1.5)/(1+Bphi*sqrt(IoS))+SUMPHI+SCATON+SANON);         #//Water osmotic coefficient
  end;
  
  function lnGamma(x, diam, z, ADH, BDH, Bdot, IoS: extended):extended;

  sqrtI: extended;
  i: integer;
  begin
  sqrtI<-sqrt(IoS);
  if z<>0 {} begin
  if diam>0 {}
  Result<-2.303*-ADH*z*z*sqrtI/(1+diam*BDH*sqrtI)+Bdot*IoS
  else
    Result<-2.303*-ADH*z*z*sqrtI/(1+sqrtI)+Bdot*IoS;
  end else
    Result<-0.2303*IoS;                               #//Bdot=0.1 (Drummond, 1981)
  end;
  
  function GibbsDurham(vij: AAE; x, G0, diam: array of extended; Aphi, ADH, BDH, Bdot: extended; m, n, min, aq, gas: integer): extended;     #//Calculates G of solution
    
  i, H2O, minaq: integer;
  Gbulk, species, gassum, IoS, lnH2O: extended;
  z: array of extended;
  begin
  H2O<-min;
  minaq<-min+aq-1;
  lnH2O<-0;
  SetLength(z,m);
  
  for i<-0 to m-1 do for (  i in 1: m-1) 
  z[i]<-vij[n-1,i];
  Gbulk<-0;
  if min>0 {}
for (  i in 1: min-1) 
  Gbulk<-Gbulk+x[i]*G0[i];
  {ShowMessage('Diam of H+ = '+FloatToStr(diam[min+1]));
    ShowMessage('IoS = '+FloatToStr(IoS));
    ShowMessage('ADH = '+FloatToStr(ADH));
    ShowMessage('BDH = '+FloatToStr(BDH));
    ShowMessage('Bdot = '+FloatToStr(Bdot));}
  species<-0;
 for (i<-min+1 in 1:minaq ) 
  species<-species+x[i];
  if (species>0) and (x[H2O]>species) {} begin
  Gbulk<-Gbulk+x[H2O]*G0[H2O]-phi(z,x,species,Aphi,H2O,minaq)*species;       #//Water from Mironenko et al., 1997
  IoS<-IonicStrength(z,x,min,minaq);
  lnH2O<-ln(x[H2O]/Mw);//*0.0180153011158497);
//ShowMessage('IoS = '+FloatToStr(IoS));
//if IoS>0.01 {} Activity(x,diam,z,ADH,BDH,Bdot,IoS,H2O,minaq);
//ShowMessage('aH+ = '+FloatToStr(diam[min+1]));
for (i<-min+1 in 1:minaq ) 
Gbulk<-Gbulk+x[i]*(G0[i]+ln(x[i])-lnH2O+lnGamma(x[i],diam[i-H2O],z[i],ADH,BDH,Bdot,IoS));                                          //Solution species from Mironenko et al., 1997
end else begin
for (i<-min in 1:minaq ) 
Gbulk<-Gbulk+x[i]*G0[i];
end;
if gas>0 {} begin
gassum<-0;

for (i<-min+1 in 1:m-1 ) 
gassum<-gassum+x[i];
if gassum>0 {}

for (i<-minaq+1 in 1:m-1 ) 
Gbulk<-Gbulk+x[i]*(G0[i]+ln(x[i]/gassum))

  for (i<-minaq+1 in 1:m-1) 
Gbulk<-Gbulk+x[i]*G0[i];
end;
Result<-Gbulk;
end;

function Divergence(x: array of extended; vij: AAE; b: extended; l, m: integer): extended;
  
i: integer;
sum: extended;
begin
sum<-0;

for (i<-0 in 1:m-1 ) 
if vij[l,i]<>0 {} sum<-sum+x[i]*vij[l,i];
Result<-sum-b;
end;

function Lagrange(Lag: array of extended; Delta: array of extended; n: integer): extended;
  
i: integer;
begin
Result<-0;

for (i<-0 in 1:n-1 ) 
Result<-Result+Lag[i]*Delta[i];
end;

function Penalty(Delta: array of extended; rp: extended; n: integer): extended;
  
i: integer;
begin
Result<-0;

for (i<-0 in 1:n-1 ) 
Result<-Result+Delta[i]*Delta[i];
Result<-rp*Result;
end;

function GibbsDurhamDerivative(x: array of extended; G0, diam, z, ADH, BDH, Bdot, IoS: extended; i, m, H2O, aq: integer): extended;
  
j: integer;
gassum, species: extended;
begin
{species<-0;

  for (j<-H2O+1 in 1:H2O+aq-1 ) 
  species<-species+x[j];}
if i<H2O+1 {}
Result<-G0
else
  if i<H2O+aq {} begin
{if x[H2O]>species {}
  Result<-G0+ln(x[i]/x[H2O])+lnGamma(x[i],diam,z,ADH,BDH,Bdot,IoS)
  else} Result<-G0+ln(x[i]/x[H2O]){+1+lnGamma(x[i],diam,z,ADH,BDH,Bdot,IoS)}
end else begin
gassum<-0;

for (j<-H2O+aq in 1:m-1 ) 
gassum<-gassum+x[j];
Result<-G0+ln(x[i]/gassum){+1};
end;
end;

function LagrangeDerivative(Lag: array of extended; vij: AAE; l, n: integer): extended;
  
i: integer;
begin
Result<-0;

for (i<-0 in 1:n-1 ) 
if vij[i,l]<>0 {} Result<-Result+vij[i,l]*Lag[i];
end;

function PenaltyDerivative(vij: AAE; Delta: array of extended; rp: extended; l, n: integer): extended;
  
i: integer;
begin
Result<-0;

for (i<-0 in 1:n-1 )
if vij[i,l]<>0 {} Result<-Result+vij[i,l]*Delta[i];
Result<-rp*Result;
end;

function MinimumSearch(vij: AAE; out x: array of extended; G0, b, diam: array of extended; IoS, eps, Aphi, ADH, BDH, Bdot: extended; m, n, min, aq, gas: integer): extended;              //Search of equilibrium composition
  
x_1, Lag, g, dk, Delta, Delta_1, z: array of extended;
Gbulk, Gbulk_old, Gbulk_1, P_new, P_old, rp, gradx, gradx_1, gradx_old, step, Betta, eps10, ddd: extended;
f, fff, i, k, l, lll, o, H2O: integer;
GreenLight: boolean;
begin
rp<-0.1;                        # //The penalty parameter
eps10<-eps/10;
H2O<-min;
SetLength(x_1,m);
SetLength(z,m);
SetLength(g,m);                  #//Gradient
SetLength(dk,m);                 #//Conjugate directions
SetLength(Lag,n);                #//Lagrange Multipliers
SetLength(Delta,n);              #//Divergence of x from b
SetLength(Delta_1,n);

for (i<-0 in 1:n-1 )
Lag[i]<-0;

for (i<-0 in 1:m-1 )
z[i]<-vij[n-1,i];
l<-0;
repeat                           //Start of the minimization loop
o<-0;
for i<-0 to n-1 do
Delta[i]<-Divergence(x,vij,b[i],i,m);
P_new<-Penalty(Delta,rp/2,n);
//ShowMessage('Point#1');
Gbulk<-GibbsDurham(vij,x,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas)+Lagrange(Lag,Delta,n)+P_new;
//ShowMessage('Point#2');
IoS<-IonicStrength(z,x,min,min+aq-1);
//ShowMessage('I#1 = '+FloatToStr(IoS));
gradx<-0;
for i<-0 in 1:m-1 do begin                 //Calculate new gradient
if i>H2O {} ddd<-diam[i-H2O] else ddd<-0;
g[i]<-GibbsDurhamDerivative(x,G0[i],ddd,z[i],ADH,BDH,Bdot,IoS,i,m,min,aq)+LagrangeDerivative(Lag,vij,i,n)+PenaltyDerivative(vij,Delta,rp,i,n);
gradx<-gradx+g[i]*g[i];
end;
gradx_1<-sqrt(gradx);
for i<-0 in 1:m-1 do
dk[i]<--g[i]/gradx_1;
lll<-0;
fff<-0;
repeat
  step<-1e-4;
Gbulk_old<-Gbulk;
Gbulk_1<-Gbulk;
GreenLight<-true;
k<-0;
repeat                                  //One-dimentional minimization (modified from Capitani and Brown, 1987)
if GreenLight=true {} begin
for i <- 0 in 1:m-1 do begin
x_1[i]<-x[i]+step*dk[i];
if x_1[i]<=0 {}
if i>min-1 {} x_1[i]<-x[i]/2
else x_1[i]<-0;
end;
for i <- 0 in 1:n-1 do
Delta_1[i]<-Divergence(x_1,vij,b[i],i,m);
Gbulk_1<-GibbsDurham(vij,x_1,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas)+Lagrange(Lag,Delta_1,n)+Penalty(Delta_1,rp/2,n);
end;
if Gbulk_1<Gbulk {} begin
step<-step+step;
for i <- 0 in 1:m-1 do
x[i]<-x_1[i];
for i <- 0 in 1:n-1 do
Delta[i]<-Delta_1[i];
Gbulk<-Gbulk_1;
GreenLight<-true;
end else begin
f<-1;
repeat
  step<-step/(f+f);
for i <- 0 in 1:m-1 do begin
x_1[i]<-x[i]+step*dk[i];
if x_1[i]<=0 {}
if i>min-1 {} x_1[i]<-x[i]/2
else x_1[i]<-0;
end;
for i <- 0 in 1:n-1 do
Delta_1[i]<-Divergence(x_1,vij,b[i],i,m);
Gbulk_1<-GibbsDurham(vij,x_1,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas)+Lagrange(Lag,Delta_1,n)+Penalty(Delta_1,rp/2,n);
if Gbulk_1<Gbulk {} GreenLight<-false;
Inc(f);
until (GreenLight=false) or (step<eps10) or (f>100);
end;
Inc(k);
until (Gbulk_1>Gbulk) or (k>100);                    #//Finish of one-dimentional minimization
if Gbulk_old-Gbulk<nmin {} begin                #//Exit from the cycle if the meaning of the function doesn't change
    Inc(lll);
    if lll>m {} break;
   end else lll<-0;
   gradx_old<-gradx;
   IoS<-IonicStrength(z,x,min,min+aq-1);
   //ShowMessage('I#2 = '+FloatToStr(IoS));
gradx<-0;
for i<-0 in 1:m-1 do begin                             #//Calculate new gradient
if i>H2O {} ddd<-diam[i-H2O] else ddd<-0;
g[i]<-GibbsDurhamDerivative(x,G0[i],ddd,z[i],ADH,BDH,Bdot,IoS,i,m,min,aq)+LagrangeDerivative(Lag,vij,i,n)+PenaltyDerivative(vij,Delta,rp,i,n);
gradx<-gradx+g[i]*g[i];
end;
gradx_1<-sqrt(gradx);
if fff>m {} begin                                  #//Refresh of the conjugate directions every n+1 steps
for i<-0 in 1:m-1 do
dk[i]<--g[i]/gradx_1;
fff<-0; end else begin
Betta<-gradx/gradx_old;                             #//Calculate new Fletcher-Reeves conjugate parameter
for i<-0 in 1:m-1 do
dk[i]<--g[i]/gradx_1+Betta*dk[i];
Inc(fff);
end;
Inc(o);
until o>10000;
P_old<-P_new;
P_new<-Penalty(Delta,rp/2,n);
if P_old>P_new {}
rp<-rp+rp;                                           #//Increase of the penalty parameter
if rp>1e+500 {} break;
for i <- 0 in 1:n-1 do                                  #//Calculate Lagrange multipliers
Lag[i]<-Lag[i]+rp*Delta[i];
Inc(l);
{if l=10 {} ShowMessage('l = '+FloatToStr(sqrt((P_new+P_new)/rp)));
if l=20 {} ShowMessage('l = '+FloatToStr(sqrt((P_new+P_new)/rp)));
if l=30 {} ShowMessage('l = '+FloatToStr(sqrt((P_new+P_new)/rp)));
if l=40 {} ShowMessage('l = '+FloatToStr(sqrt((P_new+P_new)/rp)));
if l=50 {} ShowMessage('l = '+FloatToStr(sqrt((P_new+P_new)/rp)));}
until (l>10000) or (P_new+P_new<eps*eps*rp);           #//Finish of the minimization cycle
Result<-sqrt(2*P_new/rp);
end;

function MinimumSearch_NR(vij: AAE; out x: array of extended; G0, b, diam: array of extended; IoS, eps, Aphi, ADH, BDH, Bdot: extended; m, n, min, aq, gas: integer): extended;              //Search of equilibrium composition
  
x_1, Lag, g, dk, Delta, Delta_1, z: array of extended;
Gbulk, Gbulk_old, Gbulk_1, P_new, P_old, rp, gradx, gradx_1, gradx_old, step, Betta, eps10, ddd: extended;
f, fff, i, k, l, lll, o, H2O: integer;
GreenLight: boolean;
begin
//ShowMessage('Newton-Raphson');
rp<-0.1;                         #//The penalty parameter
eps10<-eps/10;
H2O<-min;
SetLength(x_1,m);
SetLength(z,m);
SetLength(g,m);                  #//Gradient
SetLength(dk,m);                 #//Conjugate directions
SetLength(Lag,n);                #//Lagrange Multipliers
SetLength(Delta,n);              #//Divergence of x from b
SetLength(Delta_1,n);
#//for i <- 0 to m-1 do begin
#//if i=H2O {} x[i]<-Mw else x[i]<-1e-100;
#//Showmessage('x['+IntToStr(i)+'] = '+FloatToStr(x[i]));
#//end;
for i<-0 in 1:n-1 do
Lag[i]<-0;
for i<-0 in 1:m-1 do
z[i]<-vij[n-1,i];
l<-0;
repeat                           #//Start of the minimization loop
o<-0;
for i<-0 in 1:n-1 do
Delta[i]<-Divergence(x,vij,b[i],i,m);
P_new<-Penalty(Delta,rp/2,n);
#//ShowMessage('Point#1');
Gbulk<-GibbsDurham(vij,x,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas)+Lagrange(Lag,Delta,n)+P_new;
#//ShowMessage('Point#2');
IoS<-IonicStrength(z,x,min,min+aq-1);
#//ShowMessage('I#1 = '+FloatToStr(IoS));
gradx<-0;
for i<-0 in 1:m-1 do begin                 #//Calculate new gradient
if i>H2O {} ddd<-diam[i-H2O] else ddd<-0;
g[i]<-GibbsDurhamDerivative(x,G0[i],ddd,z[i],ADH,BDH,Bdot,IoS,i,m,min,aq)+LagrangeDerivative(Lag,vij,i,n)+PenaltyDerivative(vij,Delta,rp,i,n);
gradx<-gradx+g[i]*g[i];
end;
gradx_1<-sqrt(gradx);
for i<-0 in 1:m-1 do
dk[i]<--g[i]/gradx_1;    #//-Gbulk/g[i];//
  lll<-0;
fff<-0;
repeat
  step<-1e-4;
Gbulk_old<-Gbulk;
Gbulk_1<-Gbulk;
GreenLight<-true;
k<-0;
repeat                                  #//One-dimentional minimization (modified from Capitani and Brown, 1987)
if GreenLight=true {} begin
for i <- 0 in 1:m-1 do begin
x_1[i]<-x[i]+step*dk[i];
if x_1[i]<=0 {}
if i>min-1 {} x_1[i]<-x[i]/2
else x_1[i]<-0;
end;
for i <- 0 in 1:n-1 do
Delta_1[i]<-Divergence(x_1,vij,b[i],i,m);
Gbulk_1<-GibbsDurham(vij,x_1,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas)+Lagrange(Lag,Delta_1,n)+Penalty(Delta_1,rp/2,n);
end;
if Gbulk_1<Gbulk {} begin
step<-step+step;
for i <- 0 in 1:m-1 do
x[i]<-x_1[i];
for i <- 0 in 1:n-1 do
Delta[i]<-Delta_1[i];
Gbulk<-Gbulk_1;
GreenLight<-true;
end else begin
f<-1;
repeat
  step<-step/(f+f);
for i <- 0 in 1:m-1 do begin
x_1[i]<-x[i]+step*dk[i];
if x_1[i]<=0 {}
if i>min-1 {} x_1[i]<-x[i]/2
else x_1[i]<-0;
end;
for i <- 0 in 1:n-1 do
Delta_1[i]<-Divergence(x_1,vij,b[i],i,m);
Gbulk_1<-GibbsDurham(vij,x_1,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas)+Lagrange(Lag,Delta_1,n)+Penalty(Delta_1,rp/2,n);
if Gbulk_1<Gbulk {} GreenLight<-false;
Inc(f);
until (GreenLight=false) or (step<eps10) or (f>100);
end;
Inc(k);
until (Gbulk_1>Gbulk) or (k>100);                    #//Finish of one-dimentional minimization
if Gbulk_old-Gbulk<nmin {} begin                 #//Exit from the cycle if the meaning of the function doesn't change
    Inc(lll);
    if lll>m {} break;
   end else lll<-0;
   gradx_old<-gradx;
   IoS<-IonicStrength(z,x,min,min+aq-1);
   //ShowMessage('I#2 = '+FloatToStr(IoS));
gradx<-0;
for i<-0 in 1:m-1 do begin                             #//Calculate new gradient
if i>H2O {} ddd<-diam[i-H2O] else ddd<-0;
g[i]<-GibbsDurhamDerivative(x,G0[i],ddd,z[i],ADH,BDH,Bdot,IoS,i,m,min,aq)+LagrangeDerivative(Lag,vij,i,n)+PenaltyDerivative(vij,Delta,rp,i,n);
gradx<-gradx+g[i]*g[i];
end;
gradx_1<-sqrt(gradx);
if fff>m {} begin                                  #//Refresh of the conjugate directions every n+1 steps
for i<-0 in 1:m-1 do
dk[i]<--g[i]/gradx_1;   #//-Gbulk/g[i];//
  fff<-0; end else begin
Betta<-gradx/gradx_old;                             #//Calculate new Fletcher-Reeves conjugate parameter
for i<-0 in 1:m-1 do
dk[i]<--g[i]/gradx_1+Betta*dk[i];  #//-Gbulk/g[i];//
  Inc(fff);
end;
Inc(o);
until o>10000;
P_old<-P_new;
P_new<-Penalty(Delta,rp/2,n);
if P_old>P_new {}
rp<-rp+rp;                                           #//Increase of the penalty parameter
if rp>1e+500 {} break;
for i <- 0 in 1:n-1 do                                  #//Calculate Lagrange multipliers
Lag[i]<-Lag[i]+rp*Delta[i];
Inc(l);
{if l=1 {} ShowMessage('l = 1 - '+FloatToStr(sqrt((P_new+P_new)/rp)));
if l=2 {} ShowMessage('l = 2 - '+FloatToStr(sqrt((P_new+P_new)/rp)));
if l=5 {} ShowMessage('l = 5 - '+FloatToStr(sqrt((P_new+P_new)/rp)));
if l=10 {} ShowMessage('l = 10 - '+FloatToStr(sqrt((P_new+P_new)/rp)));
if l=20 {} ShowMessage('l = 20 - '+FloatToStr(sqrt((P_new+P_new)/rp)));
if l=30 {} ShowMessage('l = 30 - '+FloatToStr(sqrt((P_new+P_new)/rp)));
if l=40 {} ShowMessage('l = 40 - '+FloatToStr(sqrt((P_new+P_new)/rp)));
if l=50 {} ShowMessage('l = 50 - '+FloatToStr(sqrt((P_new+P_new)/rp)));}
until (l>1000) or (P_new+P_new<eps*eps*rp);           #//Finish of the minimization cycle
Result<-sqrt(2*P_new/rp);
//for i <- 0 in 1:n-1 do                                  #//Calculate Lagrange multipliers
//ShowMessage(FloatToStr(Lag[i]));
//ShowMessage('l = '+IntToStr(l));
//ShowMessage(FloatToStr(rp));
end;

function MinimumSearch_New(vij: AAE; out x: array of extended; G0, b, diam: array of extended; IoS, eps, Aphi, ADH, BDH, Bdot: extended; m, n, min, aq, gas: integer): extended;              //Search of equilibrium composition
  
xij, aij, xbest: AAE;
x_1, Func10, xmin, xmax: array of extended;
Kw, aaa, SumDlt_min_old, xmin_temp, xmax_temp, Gbulk, Gbulk_old, Gbulk_1, SumDlt, Dlt, SumDlt_min, smallest, highest, Func, species: extended;
i, ii, k, kkk, l, ll, Samples, H2O, H, OH: integer;
GreenLight: boolean;
begin
H2O<-min;
Samples<-200;
SetLength(aij,m,Samples);
SetLength(xij,m,Samples);
Sampling.Randomizer(aij,m,Samples,4);
SetLength(x_1,m);
SumDlt<-0;
for i<-0 in 1:n-1 do
SumDlt<-SumDlt+b[i]*b[i];
SumDlt_min<-Sqrt(SumDlt);
Gbulk<-1e+100;
SetLength(Func10,10);
for l <- 0 in 1:9 do
Func10[l]<-9999;
smallest<--20;
highest<-log10(SumDlt_min);
SetLength(xmin,m);
SetLength(xmax,m);
SetLength(xbest,m,10);
for i <- 0{min} in 1:m-1 do begin
if i=min {} xmin[i]<-SumDlt_min/10 else
  xmin[i]<-power(10,smallest);
xmax[i]<-power(10,highest);
end;
kkk<-0;
repeat
  k<-0;
repeat
  SumDlt_min_old<-SumDlt_min;
for l <- 0 in 1:Samples-1 do begin

{if kkk>0 {}
  for i <- 0 in 1:min-1 do begin
  #//xij[i,l]<-xmax[i]*aij[i,l];
  xij[i,l]<-power(10,log10(xmin[i])+(log10(xmax[i]/xmin[i]))*aij[i,l]);
  #//if xij[i,l]<1e-10 {} xij[i,l]<-1e-100;
  end else for i <- 0 in 1:min-1 do xij[i,l]<-1e-100;}

for i <- 0{min} in 1:m-1 do
if SumDlt_min<1e-3 {} xij[i,l]<-xmin[i]+(xmax[i]-xmin[i])*aij[i,l] else
  xij[i,l]<-power(10,log10(xmin[i])+(log10(xmax[i]/xmin[i]))*aij[i,l]);
xij[H,l]<-power(10,-Kw-log10(xij[OH,l]));
end;
for l <- 0 in 1:Samples-1 do begin
for i <- 0{min} in 1:m-1 do
x_1[i]<-xij[i,l];
SumDlt<-0;
for i<-0 in 1:n-1 do begin
Dlt<-Divergence(x_1,vij,b[i],i,m);
SumDlt<-SumDlt+Dlt*Dlt;
end;
SumDlt<-Sqrt(SumDlt);
if SumDlt<SumDlt_min {} begin
Gbulk_1<-GibbsDurham(vij,x_1,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
aaa<-428.5714*log10(SumDlt)+4714.286;
if (k<1) and (kkk<1) {} Func<-aaa else
  Func<-Gbulk_1+aaa;
if Func<Func10[9] {} begin
Func10[9]<-Func;
for ii <- 0{min} in 1:m-1 do
xbest[ii,9]<-x_1[ii];
for i <- 8 downin 1:0 do
if Func<Func10[i] {} begin
Func10[i+1]<-Func10[i];
Func10[i]<-Func;
for ii <- 0{min} in 1:m-1 do begin
xbest[ii,i+1]<-xbest[ii,i];
xbest[ii,i]<-x_1[ii];
end;
end;
end;
end;
end;
for i <- 0{min} in 1:m-1 do begin
xmin[i]<-power(10,highest);
xmax[i]<-power(10,smallest);
end;
for l <- 0 in 1:9 do begin
for i <- 0{min} in 1:m-1 do
x_1[i]<-xbest[i,l];
if l=0 {} Gbulk<-GibbsDurham(vij,x_1,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
SumDlt<-0;
for i<-0 in 1:n-1 do begin
Dlt<-Divergence(x_1,vij,b[i],i,m);
SumDlt<-SumDlt+Dlt*Dlt;
end;
SumDlt<-Sqrt(SumDlt);
if SumDlt<SumDlt_min {} SumDlt_min<-SumDlt;
for i <- 0{min} in 1:m-1 do begin
if xbest[i,l]<xmin[i] {} xmin[i]<-xbest[i,l];
if xbest[i,l]>xmax[i] {} xmax[i]<-xbest[i,l];
end;
end;
for i <- 0{min} in 1:m-1 do begin
xmin[i]<-StrToFloat(FormatFloat('#0.000E-00',xmin[i]));
xmax[i]<-StrToFloat(FormatFloat('#0.000E-00',xmax[i]));
end;
Inc(k);
until SumDlt_min_old<1.0001*SumDlt_min;
repeat
  GreenLight<-false;
ii<-0;
for l <- 0 in 1:8 do
for ll <- l in 1:9 do begin
for i <- 0 in 1:min-1 do
xij[i,ii]<-power(10,(log10(xbest[i,l])+log10(xbest[i,ll]))/2);
#//xij[i,l]<-(xbest[i,l]+xbest[i,ll])/2;
for i <- min in 1:m-1 do
xij[i,ii]<-power(10,(log10(xbest[i,l])+log10(xbest[i,ll]))/2);
xij[H,ii]<-power(10,-Kw-log10(xij[OH,ii]));
Inc(ii);
end;
for l <- 0 in 1:44 do begin
for i <- 0{min} in 1:m-1 do
x_1[i]<-xij[i,l];
SumDlt<-0;
for i<-0 in 1:n-1 do begin
Dlt<-Divergence(x_1,vij,b[i],i,m);
SumDlt<-SumDlt+Dlt*Dlt;
end;
SumDlt<-Sqrt(SumDlt);
if SumDlt<SumDlt_min {} begin
Gbulk_1<-GibbsDurham(vij,x_1,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
aaa<-428.5714*log10(SumDlt)+4714.286;
Func<-Gbulk_1+aaa;
if Func<Func10[9] {} begin
GreenLight<-true;
Func10[9]<-Func;
for ii <- 0{min} in 1:m-1 do
xbest[ii,9]<-x_1[ii];
for i <- 8 downin 1:0 do
if Func<Func10[i] {} begin
Func10[i+1]<-Func10[i];
Func10[i]<-Func;
for ii <- 0{min} in 1:m-1 do begin
xbest[ii,i+1]<-xbest[ii,i];
xbest[ii,i]<-x_1[ii];
end;
end;
end;
end;
end;
until GreenLight=false;

{if SumDlt_min>1e-3 {} begin
  ShowMessage('k = '+IntToStr(k));
  for i <- 0 in 1:m-1 do begin
  #//ShowMessage('G['+IntToStr(i)+'] = '+FloatToStr(G0[i]));
  ShowMessage('x['+IntToStr(i)+'] = '+FloatToStr(xbest[i,0]));
  end;
  ShowMessage('SumDlt = '+FloatToStr(SumDlt_min));
  ShowMessage('Gbulk = '+FloatToStr(Gbulk));
  end;}

for i <- 0{min} in 1:m-1 do begin
if xbest[i,0]>SumDlt_min {} xmin[i]<-xbest[i,0]-SumDlt_min else
  xmin[i]<-power(10,smallest);
if xbest[i,0]<SumDlt_min {} xmax[i]<-SumDlt_min else
  xmax[i]<-xbest[i,0]+SumDlt_min;
end;
Inc(kkk);
until (kkk>3{20}) or (SumDlt_min<eps);
for i <- 0{min} in 1:m-1 do
x[i]<-xbest[i,0];
Result<-SumDlt_min;
end;

function MinimumSearch_New2(vij: AAE; out x: array of extended; G0, b: array of extended; IoS, eps, Aphi, Kw: extended; m, n, min, aq, gas, H, OH: integer): extended;              //Search of equilibrium composition
  
xij, aij, xbest: AAE;
x_1, Func10, xmin, xmax: array of extended;
MinAct: array of boolean;
aaa, SumDlt_min_old, xmin_temp, xmax_temp, Gbulk, Gbulk_old, Gbulk_1, SumDlt, Dlt, SumDlt_min, smallest, highest, Func, species, Penalty: extended;
i, ii, k, kkk, l, ll, Samples, H2O: integer;
GreenLight: boolean;
begin
H2O<-min;
Samples<-100;
SetLength(aij,m,Samples);
SetLength(xij,m,Samples);
Sampling.Randomizer(aij,m,Samples,4);
SetLength(x_1,m);
SumDlt<-0;
for i<-0 in 1:n-1 do
SumDlt<-SumDlt+b[i]*b[i];
SumDlt_min<-Sqrt(SumDlt);
Gbulk<-1e+100;
SetLength(Func10,10);
for l <- 0 in 1:9 do
Func10[l]<-9999;
smallest<--20;
highest<-log10(SumDlt_min);
SetLength(xmin,m);
SetLength(xmax,m);
SetLength(xbest,m,10);
for i <- 0{min} in 1:m-1 do begin
if i=min {} xmin[i]<-SumDlt_min/10 else
  xmin[i]<-power(10,smallest);
xmax[i]<-power(10,highest);
end;
Penalty<-MinimumSearch(vij,x,G0,b,diam,IoS,eps,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
SumDlt<-0;
for i<-0 in 1:n-1 do begin
Dlt<-Divergence(x,vij,b[i],i,m);
SumDlt<-SumDlt+Dlt*Dlt;
end;
SumDlt_min<-Sqrt(SumDlt);
Gbulk<-GibbsDurham(vij,x,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
SetLength(MinAct,m);
for i <- 0 in 1:m-1 do
if (i<min) and (x[i]<1e-10) {} begin MinAct[i]<-false; x[i]<-nmin; end else MinAct[i]<-true;
for i <- 0 in 1:m-1 do
if MinAct[i]=true {} begin
if x[i]>SumDlt_min {} xmin[i]<-x[i]-SumDlt_min else
  xmin[i]<-power(10,smallest);
if x[i]<SumDlt_min {} xmax[i]<-SumDlt_min else
  xmax[i]<-x[i]+SumDlt_min;
end else begin
xmin[i]<-nmin;
xmax[i]<-nmin;
end;
aaa<-428.5714*log10(SumDlt_min)+4714.286;
Func<-Gbulk+aaa{*10};
Func10[0]<-Func;
for i <- 0 in 1:m-1 do
xbest[i,0]<-x[i];
{for i <- 0 in 1:m-1 do begin
  #//ShowMessage('G['+IntToStr(i)+'] = '+FloatToStr(G0[i]));
  ShowMessage('xmin['+IntToStr(i)+'] = '+FloatToStr(xmin[i]));
  ShowMessage('xmax['+IntToStr(i)+'] = '+FloatToStr(xmax[i]));
  end;
  ShowMessage('SumDlt = '+FloatToStr(SumDlt_min));}
for i <- 1 in 1:9 do
for ii <- 0 in 1:m-1 do
xbest[ii,i]<-nmin;
#//ShowMessage('Gbulk = '+FloatToStr(Gbulk));

kkk<-0;
repeat
  k<-0;
repeat
  SumDlt_min_old<-SumDlt_min;
for l <- 0 in 1:Samples-1 do begin

{if kkk>0 {}
  for i <- 0 in 1:min-1 do begin
  #//xij[i,l]<-xmax[i]*aij[i,l];
  xij[i,l]<-power(10,log10(xmin[i])+(log10(xmax[i]/xmin[i]))*aij[i,l]);
  #//if xij[i,l]<1e-10 {} xij[i,l]<-1e-100;
  end else for i <- 0 in 1:min-1 do xij[i,l]<-1e-100;}

for i <- 0{min} in 1:min-1 do
if MinAct[i]=true{SumDlt_min<1e-3} {} begin if SumDlt_min<1e-3 {} xij[i,l]<-xmin[i]+(xmax[i]-xmin[i])*aij[i,l] else
  xij[i,l]<-power(10,log10(xmin[i])+(log10(xmax[i]/xmin[i]))*aij[i,l]);
end else xij[i,l]<-nmin;
for i <- min in 1:m-1 do
if SumDlt_min<1e-3 {} xij[i,l]<-xmin[i]+(xmax[i]-xmin[i])*aij[i,l] else
  xij[i,l]<-power(10,log10(xmin[i])+(log10(xmax[i]/xmin[i]))*aij[i,l]);
xij[H,l]<-power(10,-Kw-log10(xij[OH,l]));
end;
for l <- 0 in 1:Samples-1 do begin
for i <- 0{min} in 1:m-1 do
x_1[i]<-xij[i,l];
SumDlt<-0;
for i<-0 in 1:n-1 do begin
Dlt<-Divergence(x_1,vij,b[i],i,m);
SumDlt<-SumDlt+Dlt*Dlt;
end;
SumDlt<-Sqrt(SumDlt);
if SumDlt<SumDlt_min {} begin
Gbulk_1<-GibbsDurham(vij,x_1,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
aaa<-428.5714*log10(SumDlt)+4714.286;
#//if (k<1) and (kkk<1) {} Func<-aaa else
  Func<-Gbulk_1+aaa{*10};
if Func<Func10[9] {} begin
Func10[9]<-Func;
for ii <- 0{min} in 1:m-1 do
xbest[ii,9]<-x_1[ii];
for i <- 8 downin 1:0 do
if Func<Func10[i] {} begin
Func10[i+1]<-Func10[i];
Func10[i]<-Func;
for ii <- 0{min} in 1:m-1 do begin
xbest[ii,i+1]<-xbest[ii,i];
xbest[ii,i]<-x_1[ii];
end;
end;
end;
end;
end;
#//if k=0 {} ShowMessage('It works#1');
for i <- 0{min} in 1:m-1 do begin
xmin[i]<-power(10,highest);
xmax[i]<-power(10,smallest);
end;
for l <- 0 in 1:9 do begin
for i <- 0{min} in 1:m-1 do
x_1[i]<-xbest[i,l];
if l=0 {} Gbulk<-GibbsDurham(vij,x_1,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
SumDlt<-0;
for i<-0 in 1:n-1 do begin
Dlt<-Divergence(x_1,vij,b[i],i,m);
SumDlt<-SumDlt+Dlt*Dlt;
end;
SumDlt<-Sqrt(SumDlt);
if SumDlt<SumDlt_min {} SumDlt_min<-SumDlt;
for i <- 0{min} in 1:m-1 do begin
if xbest[i,l]<xmin[i] {} xmin[i]<-xbest[i,l];
if xbest[i,l]>xmax[i] {} xmax[i]<-xbest[i,l];
end;
end;
#//if k=0 {} ShowMessage('It works#2');
for i <- 0{min} in 1:m-1 do begin
xmin[i]<-StrToFloat(FormatFloat('#0.000E-00',xmin[i]));
xmax[i]<-StrToFloat(FormatFloat('#0.000E-00',xmax[i]));
end;
#//if k=0 {} ShowMessage('It works#3');
Inc(k);
until SumDlt_min_old<1.0001*SumDlt_min;
repeat
  GreenLight<-false;
ii<-0;
for l <- 0 in 1:8 do
for ll <- l in 1:9 do begin
#//for i <- 0 in 1:m-1 do begin
#//ShowMessage('xbest['+IntToStr(i)+','+IntToStr(l)+'] = '+FloatToStr(xbest[i,l]));
#//ShowMessage('xbest['+IntToStr(i)+','+IntToStr(ll)+'] = '+FloatToStr(xbest[i,ll]));
#//end;
for i <- 0 in 1:min-1 do
xij[i,ii]<-power(10,(log10(xbest[i,l])+log10(xbest[i,ll]))/2);
#//xij[i,l]<-(xbest[i,l]+xbest[i,ll])/2;
for i <- min in 1:m-1 do
xij[i,ii]<-power(10,(log10(xbest[i,l])+log10(xbest[i,ll]))/2);
#//xij[i,l]<-(xbest[i,l]+xbest[i,ll])/2;
xij[H,ii]<-power(10,-Kw-log10(xij[OH,ii]));
Inc(ii);
end;
#//ShowMessage('It works#5');
for l <- 0 in 1:44 do begin
for i <- 0{min} in 1:m-1 do
x_1[i]<-xij[i,l];
SumDlt<-0;
for i<-0 in 1:n-1 do begin
Dlt<-Divergence(x_1,vij,b[i],i,m);
SumDlt<-SumDlt+Dlt*Dlt;
end;
SumDlt<-Sqrt(SumDlt);
if SumDlt<SumDlt_min {} begin
Gbulk_1<-GibbsDurham(vij,x_1,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
aaa<-428.5714*log10(SumDlt)+4714.286;
Func<-Gbulk_1+aaa{*10};
if Func<Func10[9] {} begin
GreenLight<-true;
Func10[9]<-Func;
for ii <- 0{min} in 1:m-1 do
xbest[ii,9]<-x_1[ii];
for i <- 8 downin 1:0 do
if Func<Func10[i] {} begin
Func10[i+1]<-Func10[i];
Func10[i]<-Func;
for ii <- 0{min} in 1:m-1 do begin
xbest[ii,i+1]<-xbest[ii,i];
xbest[ii,i]<-x_1[ii];
end;
end;
end;
end;
end;
until GreenLight=false;

{if SumDlt_min>1e-3 {} begin
  ShowMessage('k = '+IntToStr(k));
  for i <- 0 in 1:m-1 do begin
 # //ShowMessage('G['+IntToStr(i)+'] = '+FloatToStr(G0[i]));
  ShowMessage('x['+IntToStr(i)+'] = '+FloatToStr(xbest[i,0]));
  end;
  ShowMessage('SumDlt = '+FloatToStr(SumDlt_min));
  ShowMessage('Gbulk = '+FloatToStr(Gbulk));}
#//end;

for i <- 0{min} in 1:m-1 do begin
if xbest[i,0]>SumDlt_min {} xmin[i]<-xbest[i,0]-SumDlt_min else
  xmin[i]<-power(10,smallest);
if xbest[i,0]<SumDlt_min {} xmax[i]<-SumDlt_min else
  xmax[i]<-xbest[i,0]+SumDlt_min;
end;
Inc(kkk);
until (kkk>20) or (SumDlt_min<eps);
for i <- 0{min} in 1:m-1 do
x[i]<-xbest[i,0];
Result<-SumDlt_min;
end;

procedure LinSolver(a: AAE; c: array of extended; out x: array of extended; n: integer);     #//de Capitani, 1987; de Capitani and Brown, 1987
{#***********************************************************
  #  solves a system of n linear equations with n variables:
    x(1)*a(1,1) + x(2)*a(1,2) + ... + x(n)*a(1,n) = c(1)
  x(1)*a(2,1) + x(2)*a(2,2) + ... + x(n)*a(2,n) = c(2)
  ...
  x(1)*a(n,1) + x(2)*a(n,2) + ... + x(n)*a(n,n) = c(n)
 # ************************************************************}
  
ColSum, aTemp: array of extended;
ColNumber: array of integer;
FocTemp, FocMax, cTemp: extended;
i, ii, Focus, Row, Col, RowTemp, ColTemp: integer;
begin
for i<-0 in 1:n-1 do x[i]<-0;
SetLength(ColSum,n);
SetLength(ColNumber,n);
SetLength(aTemp,n);
for row<-0 in 1:n-1 do
ColNumber[Row]<-Row;                  #  //Initial numeration in the matrix
for Focus<-0 in 1:n-1 do begin
for Col<-Focus in 1:n-1 do begin
ColSum[Col]<-0;
for Row<-0 in 1:n-1 do
ColSum[Col]<-ColSum[Col]+abs(a[Row,Col]);
end;
RowTemp<--1;
ColTemp<--1;
FocMax<-0;
for Row<-Focus in 1:n-1 do               # //Search for the max value
for Col<-Focus in 1:n-1 do
if ColSum[col]>0 {} begin
FocTemp<-abs(a[Row,Col]/ColSum[Col]);
if FocTemp>FocMax {} begin
FocMax<-FocTemp;
RowTemp<-Row;
ColTemp<-Col;
end;
end;
if (RowTemp=-1) or (ColTemp=-1) {} break;
if RowTemp<>Focus {} begin           # //Switch rows
cTemp<-c[RowTemp];
for i<-0 in 1:n-1 do
aTemp[i]<-a[RowTemp,i];
c[RowTemp]<-c[Focus];
for i<-0 in 1:n-1 do
a[RowTemp,i]<-a[Focus,i];
c[Focus]<-cTemp;
for i<-0 in 1:n-1 do
a[Focus,i]<-aTemp[i];
end;
if ColTemp<>Focus {} begin           # //Switch columns
ii<-ColNumber[ColTemp];
for i<-0 in 1:n-1 do
aTemp[i]<-a[i,ColTemp];
ColNumber[ColTemp]<-ColNumber[Focus];
for i<-0 in 1:n-1 do
a[i,ColTemp]<-a[i,Focus];
ColNumber[Focus]<-ii;
for i<-0 in 1:n-1 do
a[i,Focus]<-aTemp[i];
end;
FocTemp<-a[Focus,Focus];
c[Focus]<-c[Focus]/FocTemp;
for i<-0 in 1:n-1 do
a[Focus,i]<-a[Focus,i]/FocTemp;
for ii<-0 in 1:n-1 do
if ii<>Focus {} begin
FocTemp<-a[ii,Focus];
c[ii]<-c[ii]-c[Focus]*FocTemp;
for i<-0 in 1:n-1 do
a[ii,i]<-a[ii,i]-a[Focus,i]*FocTemp;
end;
end;
for Col<-0 in 1:n-1 do
x[ColNumber[Col]]<-c[Col];              #//The solution of the system of n linear equations with n variables
end;

procedure PhaseExchange(k, i, m, n: integer; MP, SB: AAE; out Gsol: array of extended; out Title: array of integer);  # //Replaces components (de Capitani, 1987; de Capitani and Brown, 1987)
  
j, Title_temp: integer;
Gsol_temp: extended;
SB_temp, MP_temp: array of extended;
begin
SetLength(SB_temp,m);
SetLength(MP_temp,m);
for j<-0 in 1:n-1 do
SB_temp[j]<-SB[k,j];
Gsol_temp<-Gsol[k];
Title_temp<-Title[k];
for j<-0 in 1:m-1 do
MP_temp[j]<-MP[k,j];
for j<-0 in 1:n-1 do
SB[k,j]<-SB[i,j];
Gsol[k]<-Gsol[i];
Title[k]<-Title[i];
for j<-0 in 1:m-1 do
MP[k,j]<-MP[i,j];
for j<-0 in 1:n-1 do
SB[i,j]<-SB_temp[j];
Gsol[i]<-Gsol_temp;
Title[i]<-Title_temp;
for j<-0 in 1:m-1 do
MP[i,j]<-MP_temp[j];
end;

procedure Reduce(SB: AAE; k, n, Nmax: integer);    # //Transform in 1:the simplex matrix (de Capitani, 1987; de Capitani and Brown, 1987)
  
i, ii: integer;
Focus: extended;
AR: array of extended;
begin
SetLength(AR,n);
Focus<-SB[k,k];
if Focus=0 {}
for i<-0 in 1:Nmax do       # //Choose all x in the row
SB[i,k]<-0
else begin
for i<-0 in 1:Nmax do
SB[i,k]<-SB[i,k]/Focus;
for i<-0 in 1:n-1 do
AR[i]<-SB[k,i];
for i<-0 in 1:n-1 do
if (i<>k) and (AR[i]<>0) {}
for ii<-0 in 1:Nmax do
SB[ii,i]<-SB[ii,i]-SB[ii,k]*AR[i];
end;
for i<-0 in 1:n-1 do         # //Makes 0-matrix with diagonal = 1
SB[k,i]<-0;
SB[k,k]<-1;
end;

procedure FullRed(MP, SB: AAE; out Gsol: array of extended; out Title: array of integer; m, n, Nmax: integer);    #//Randing of the matrix, modified from de Capitani, 1987; de Capitani and Brown, 1987
  
ActCell, SolCol, ElRow, SolCol_max, ElRow_max: integer;
FX, F1: extended;
ColSum: array of extended;
FF: array of extended;
begin
SetLength(ColSum,n);
SetLength(FF,Nmax+1);
for ActCell<-0 in 1:n-1 do begin             # //Goes through all elements in the matrix
for SolCol<-ActCell in 1:n-1 do begin       # //Calculates the sum of each element in all solutions (ColSum = balance)
ColSum[SolCol]<-0;
for ElRow<-0 in 1:n-1 do
ColSum[SolCol]<-ColSum[SolCol]+abs(SB[SolCol,ElRow]);
end;
SolCol_max<-ActCell;
ElRow_max<-ActCell;
FX<-0;
for ElRow<-ActCell in 1:n-1 do              # //Search for the element with a highest portion in the bulk balance
for SolCol<-ActCell in 1:n-1 do begin
F1<-0;
if ColSum[SolCol]<>0 {}
F1<-abs(SB[SolCol,ElRow]/ColSum[SolCol]);
if F1>FX {} begin
SolCol_max<-SolCol;
ElRow_max<-ElRow;
FX<-F1;
end;
end;
if ElRow_max<>ActCell {} begin             #//Check or move the highest portion in the upper position. Therefore, the highest values are in the top of the matrix
for SolCol<-0 in 1:Nmax do
FF[SolCol]<-SB[SolCol,ElRow_max];
for SolCol<-0 in 1:Nmax do begin
SB[SolCol,ElRow_max]<-SB[SolCol,ActCell];
SB[SolCol,ActCell]<-FF[SolCol];
end;
end;
if SolCol_max<>ActCell {}                                               # //Moves parameters of solutions
PhaseExchange(SolCol_max,ActCell,m,n,MP,SB,Gsol,Title);
Reduce(SB,ActCell,n,Nmax);                                                # //Recalculation of matrix
end;
end;

function MinimumSearch_CB(vij: AAE; out x: array of extended; G0, b, diam: array of extended; IoS, eps, Aphi, ADH, BDH, Bdot: extended; m, n, min, aq, gas: integer): extended;     # //Search for equilibrium composition, modified from de Capitani, 1987; de Capitani and Brown, 1987
  
i, ii, iii, i_worse, i_worst, j, k, Nmax, lim, N_Gauss, iPos, iRepos, iRepos_temp, mainloop: integer;
Gbulk, Gbulk_old, BulkSum, negl, xmin, sum, Penalty: extended;
D, DNew, DNew_temp, MP_temp, MP_aq, Gsol, b_aq, G00: array of extended;
MP, SB, SB_2, SB_3: AAE;
Title: array of integer;
TimeStart: TDateTime;
begin
lim<-10;                               # //The limiting number of iterations during which the calculated phase is saving in the system
xmin<-nmin;
N_Gauss<-Round(((n-1)*n/2+1)*lim+m);  #  //The maximum possible number of phases in the system (calculated with the Gauss formula for sum - n*(n+1)/2)
SetLength(D,n);
SetLength(DNew,n);
SetLength(DNew_temp,n);
SetLength(b_aq,n);
SetLength(MP_temp,m);
SetLength(MP_aq,m);
SetLength(MP,N_Gauss,m);
SetLength(SB,N_Gauss,n);
SetLength(SB_2,N_Gauss,n);
SetLength(SB_3,n,N_Gauss);
SetLength(Gsol,N_Gauss);
SetLength(G00,m);
SetLength(Title,N_Gauss);
BulkSum<-0;
for i<-0 in 1:n-1 do
BulkSum<-BulkSum+b[i];
for j<-0 in 1:m-1 do
G00[j]<-G0[j];
Gbulk<-1e+100;
Nmax<--1;
for j<- min in 1:m-1 do begin       #//Add end-members; water and aqueous species
//j<- min;
Inc(Nmax);
for i<-0 in 1:m-1 do
if i<>j {}
MP[Nmax,i]<-xmin
else
  MP[Nmax,i]<-1;
if j<min {} Title[Nmax]<--i-1 else Title[Nmax]<-0;
Gsol[Nmax]<-G0[j];
end;
mainloop<-1;
repeat                           # //Start of the Main Loop

{if mainloop=1 {}
  for i<-0 in 1:min-1 do begin    # //Add minerals on the first step
  Inc(Nmax);
  for ii<-0 in 1:m-1 do
  if ii<>i {}
  MP[Nmax,ii]<-xmin
  else
    MP[Nmax,ii]<-1;
  Title[Nmax]<--i-1;
  Gsol[Nmax]<-G0[i];
  end;}

for i<-0 in 1:Nmax do
for ii<-0 in 1:n-1 do begin       #//Add components in 1:the matrix
SB[i,ii]<-0;
for iii<-0 in 1:m-1 do
SB[i,ii]<-SB[i,ii]+MP[i,iii]*vij[ii,iii];
end;
FullRed(MP,SB,Gsol,Title,m,n,Nmax);

#//if Nmax>n {}
#//if mainloop>1 {}

for k<-n in 1:Nmax do begin      #  //Start of the Secondary Loop (loop k)
sum<-0;
for j<-0 in 1:n-1 do
sum<-sum+Gsol[j]*SB[k,j];
if sum-Gsol[k]>0 {}  begin    #//If the new phase (k) is thermodynamically preferably
negl<-1E-14;
i_worse<--1;
i_worst<--1;
iRepos_temp<-0;
repeat                        # //Test all set mo more than 3 times with increasing 'negl'
for i<-0 in 1:n-1 do begin
for ii<-0 in 1:n-1 do begin
SB_2[i,ii]<-0;
for iii<-0 in 1:m-1 do
SB_2[i,ii]<-SB_2[i,ii]+MP[i,iii]*vij[ii,iii];
end;
end;
for ii<-0 in 1:n-1 do begin
SB_2[k,ii]<-0;
for iii<-0 in 1:m-1 do
SB_2[k,ii]<-SB_2[k,ii]+MP[k,iii]*vij[ii,iii];
end;
for j<-0 in 1:n-1 do begin     # //Start of the Decision Loop (j loop)
for i<-0 in 1:n-1 do
if i=j {}
for ii<-0 in 1:n-1 do
SB_3[ii,i]<-SB_2[k,ii]
else
  for ii<-0 in 1:n-1 do
SB_3[ii,i]<-SB_2[i,ii];
LinSolver(SB_3,b,DNew,n);# //Calculate a solution of the linear system (a core of this procedure)
iPos<-0;
iRepos<-0;
for i<-0 in 1:n-1 do
if (DNew[i]/BulkSum>-negl) and (DNew[i]<BulkSum) {} begin Inc(iPos);   # //Sufficient decisions
if DNew[i]>nmin {} Inc(iRepos);                                       # //Positive decisions
end;
if (iPos=n) and (iRepos=n) {} begin            #  //If all sufficient and none-negative
i_worst<-j;
break;                                           # //{} exit from this iteration in 1:save the decision
end;                                             #  //If some of decisions are negative
if (iPos=n) and (iRepos>iRepos_temp) {} begin   # //but the number of none-negative decisions are much than on the previous step
i_worse<-j;                                     #  //remember the number of successful position
IREPOS_temp<-iRepos;                            #  //remember the number of none-negative values
for i<-0 in 1:n-1 do
DNew_temp[i]<-DNew[i];                         #  //remember the successful decision
end;
end;                                               # //end of the Decision Loop (j loop)
negl<-negl*10;
until (i_worst>-1) or (negl>1E-12);                 # //if decision was found successfully or was not found during 3 times
if (i_worst<0) and (i_worse>-1) {} begin          # //the better decision from bad ones
i_worst<-i_worse;                                  # //remember the successful position in 1:do out from the loop
for i<-0 in 1:n-1 do
DNew[i]<-DNew_temp[i];                            # //insert the best found decision
end;
if i_worst>-1 {} begin                           #  //if the position in the matrix with good decision (DNew) was found
PhaseExchange(k,i_worst,m,n,MP,SB,Gsol,Title);      #//Exchange of found phases
for i<-0 in 1:n-1 do begin
D[i]<-abs(DNew[i]);
DNew[i]<-0;
end;
Reduce(SB,i_worst,n,Nmax);
end;
end;     # //End of the 'if sum-Gsol[k]>0'
end;     #  //End of the Secondary Loop (k loop)
Gbulk_old<-Gbulk;
Gbulk<-0;
for i<-0 in 1:n-1 do
Gbulk<-Gbulk+Gsol[i]*D[i];  # //Current Gibbs energy of the system
for i<-0 in 1:n-1 do
for ii<-0 in 1:n-1 do begin
SB[i,ii]<-0;
for iii<-0 in 1:m-1 do
SB[i,ii]<-SB[i,ii]+MP[i,iii]*vij[ii,iii];
end;

if mainloop>1 {} begin
for i<-Nmax downin 1:0 do begin         # //step 3
sum<-0;
for ii<-0 in 1:n-1 do
sum<-sum+SB[i,ii]*Gsol[ii];
Gsol[i]<-Gsol[i]-sum;
end;
for i<-Nmax in 1:0 do
if Title[i]<1 {}
for ii<-0 in 1:m-1 do
if MP[i,ii]=1 {} G0[ii]<-Gsol[i];
end;

Penalty<-0;
for i<-0 in 1:n-1 do begin
sum<-0;
b_aq[i]<-0;
for ii<-0 in 1:n-1 do begin
sum<-sum+SB[ii,i]*D[ii];
if Title[ii]<0 {}                     # //Calculate the sum of minerals
b_aq[i]<-b_aq[i]+SB[ii,i]*D[ii];
end;
b_aq[i]<-abs(b[i]-b_aq[i]);               # //The rest is the aqueous solution
Penalty<-Penalty+sqr(sum-b[i]);
end;
Penalty<-sqrt(Penalty);
for i<-0 in 1:m-1 do begin                   # //Calculate the current compositions of the bulk system and aqueous solution
x[i]<-0;
MP_aq[i]<-0;
for ii<-0 in 1:n-1 do begin
if D[ii]>eps {}                         #//Output the current result without a small addition from end-members of aqueous solution in 1:equilibrate the task (important in the case of oxygen)
x[i]<-x[i]+MP[ii,i]*D[ii];
if Title[ii]>-1 {}
MP_aq[i]<-MP_aq[i]+MP[ii,i]*D[ii];      # //Starting composition for gradient search
end;
if (i<min) and (x[i]<1e-20) {} x[i]<-0;  #//Delete very small values for minerals
if MP_aq[i]=0 {}
if i=min {} MP_aq[i]<-1 else MP_aq[i]<-xmin;  # //Sometimes the decision can't be found (D[i]=0). In this case, MP_aq[i]=0 collapse the procedure MinimumSearch
   #//ShowMessage('x['+IntToStr(i)+'] = '+FloatToStr(x[i]));
  # //ShowMessage('G0['+IntToStr(i)+'] = '+FloatToStr(G0[i]));
  end;
  if (abs(Gbulk_old-Gbulk)<eps) and (mainloop>5) and (Penalty<eps) {} break;                 #//Exit from the minimization procedure
  i<-0;
  repeat                               # //Delete old guesses from the system
   if i<n {}
    if Title[i]>0 {} Title[i]<-1;    # //All solutions incoming in the decision are good
   if i>n-1 {} begin
    if Title[i]>lim {} begin         # //If this solution was not used during lim=10 iterations, it should be deleted
     PhaseExchange(Nmax,i,m,n,MP,SB,Gsol,Title);
     i<-i-1;
     Nmax<-Nmax-1;
    end else if Title[i]>0 {} Inc(Title[i]);
   end;
   Inc(i);
  until i>Nmax;
  {for i<-0 in 1:n-2 do                    #//Add average phases (this increases the robustness of the procedure)
   if Title[i]>0{-1} {{}
    for ii<-i+1 in 1:n-1 do
     if Title[ii]>0{-1} {{} begin
      Inc(Nmax);
      for j<-min in 1:m-1 do
       MP_temp[j]<-(MP[i,j]+MP[ii,j])/2;#//power(10,(log10(MP[i,j])+log10(MP[ii,j]))/2);   //
      sum<-0;
      for j<-min in 1:m-1 do
       sum<-sum+MP_temp[j];
      for j<-0 in 1:m-1 do begin
       MP_temp[j]<-max(xmin,MP_temp[j]/sum);
       MP[Nmax,j]<-MP_temp[j];
      end;
      Title[Nmax]<-1;
      Gsol[Nmax]<-GibbsDurham(vij,MP_temp,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
     end;}

  #//Penalty<-MinimumSearch_NR(vij,MP_aq,G00,b,diam,IoS,eps,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);      # //Calculates the solution composition with Conjugate Directions method
  Inc(Nmax);                                                              # //Add found solution in 1:the matrix
  Title[Nmax]<-1;
  MP_temp[0]<-0;
  MP_temp[1]<-0;
  MP_temp[2]<-0.999998;
  MP_temp[3]<-0.000000003585;
  MP_temp[4]<-5.987E-14;
  MP_temp[5]<-1E-66;
  MP_temp[6]<-0.0000000009169;
  MP_temp[7]<-0.000001803;
  MP_temp[8]<-0.000000002635;
  {sum<-0;
  for ii<-min in 1:m-1 do
   sum<-sum+MP_aq[ii];
  for ii<-0 in 1:min-1 do
   MP_temp[ii]<-xmin;
  for ii<-min in 1:m-1 do
   MP_temp[ii]<-max(xmin,MP_aq[ii]/sum);}
  Gsol[Nmax]<-GibbsDurham(vij,MP_temp,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
  for ii<-0 in 1:m-1 do
   MP[Nmax,ii]<-MP_temp[ii];

  if mainloop=1 {}
   for i<-0 in 1:min-1 do begin    # //Add minerals on the first step
    Inc(Nmax);
    for ii<-0 in 1:m-1 do
     if ii<>i {}
      MP[Nmax,ii]<-xmin
     else
      MP[Nmax,ii]<-1;
    Title[Nmax]<--i-1;
    Gsol[Nmax]<-G0[i];
   end;
  //ShowMessage('mainloop = '+IntToStr(mainloop));
  Inc(mainloop);
 until mainloop>100;             # //End of the Main Loop
 {if mainloop>100 {}} ShowMessage('mainloop = '+IntToStr(mainloop));
 Result<-Penalty;
end;






























procedure MuCalc(vij: AAE; x: array of extended; out mu: array of extended; G0, diam: array of extended; Aphi, ADH, BDH, Bdot: extended; m, n, H2O, aq, gas: integer);   #  //Calculates mu
  
 i, aqsol: integer;
 species, gassum, IoS, lnH2O: extended;
 z: array of extended;
begin
 aqsol<-H2O+aq-1;
 lnH2O<-0;
 SetLength(z,m);
 for i<-H2O+1 in 1:aqsol do
  z[i]<-vij[n-1,i];
 species<-0;
 for i<-H2O+1 in 1:aqsol do
  species<-species+x[i];
 if (species>0) and (x[H2O]>species) {} begin
  mu[H2O]<-G0[H2O]-phi(z,x,species,Aphi,H2O,aqsol)*species/x[H2O];
  IoS<-IonicStrength(z,x,H2O,aqsol);
  lnH2O<-ln(x[H2O]/Mw);
  for i<-H2O+1 in 1:aqsol do
   mu[i]<-G0[i]+ln(x[i])-lnH2O+lnGamma(x[i],diam[i-H2O],z[i],ADH,BDH,Bdot,IoS);
 end else
  for i<-H2O in 1:aqsol do
   mu[i]<-G0[i];
 if gas>0 {} begin
  gassum<-0;
  for i<-aqsol+1 in 1:m-1 do
   gassum<-gassum+x[i];
  if gassum>0 {}
   for i<-aqsol+1 in 1:m-1 do
    mu[i]<-G0[i]+ln(x[i]/gassum)
  else
   for i<-aqsol+1 in 1:m-1 do
    mu[i]<-G0[i];
 end;
end;

function MuGD(vij: AAE; x: array of extended; out mu: array of extended; G0, diam: array of extended; Aphi, ADH, BDH, Bdot: extended; m, n, min, aq, gas: integer): extended;    # //Calculates G of solution
  
 i, H2O, minaq: integer;
 Gbulk, species, gassum, IoS, lnH2O: extended;
 z: array of extended;
begin
 H2O<-min;
 minaq<-min+aq-1;
 lnH2O<-0;
 SetLength(z,m);
 for i<-0 in 1:m-1 do
  z[i]<-vij[n-1,i];
 Gbulk<-0;
 if min>0 {}
  for i <- 0 in 1:min-1 do begin
   mu[i]<-G0[i];
   Gbulk<-Gbulk+x[i]*G0[i];
  end;
 species<-0;
 for i<-min+1 in 1:minaq do
  species<-species+x[i];
 if (species>0) and (x[H2O]>species) {} begin
  mu[H2O]<-G0[H2O]-phi(z,x,species,Aphi,H2O,minaq)*species/x[H2O];
 # //ShowMessage(FloatToStr(phi(z,x,species,Aphi,H2O,minaq)));
  Gbulk<-Gbulk+x[H2O]*mu[H2O];       //Water from Mironenko et al., 1997
  IoS<-IonicStrength(z,x,min,minaq);
  lnH2O<-ln(x[H2O]/Mw);//*0.0180153011158497);
  for i<-min+1 in 1:minaq do begin
   mu[i]<-G0[i]+ln(x[i])-lnH2O+lnGamma(x[i],diam[i-H2O],z[i],ADH,BDH,Bdot,IoS);
   Gbulk<-Gbulk+x[i]*mu[i];                                          //Solution species from Mironenko et al., 1997
  end;
 end else begin
  for i<-min in 1:minaq do begin
   mu[i]<-G0[i];
   Gbulk<-Gbulk+x[i]*mu[i];
  end;
 end;
 if gas>0 {} begin
  gassum<-0;
  for i<-minaq+1 in 1:m-1 do
   gassum<-gassum+x[i];
  if gassum>0 {}
   for i<-minaq+1 in 1:m-1 do begin
    mu[i]<-G0[i]+ln(x[i]/gassum);
    Gbulk<-Gbulk+x[i]*mu[i];
  end else
   for i<-minaq+1 in 1:m-1 do begin
    mu[i]<-G0[i];
    Gbulk<-Gbulk+x[i]*mu[i];
   end;
 end;
 Result<-Gbulk;
end;

function MuxNew(vij: AAE; out x: array of extended; out mu: array of extended; G0, dk, diam: array of extended; step, Aphi, ADH, BDH, Bdot: extended; m, n, min, aq, gas: integer): extended;
  
 xsum: extended;
 i: integer;
begin
 xsum<-0;
 for i <- min in 1:m-1 do begin
  x[i]<-max(0,x[i]+step*dk[i]);
  xsum<-xsum+x[i];
 end;
 for i <- min in 1:m-1 do
  x[i]<-max(nmin,x[i]/xsum);
 Result<-MuGD(vij,x,mu,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas); //GibbsDurham
end;

procedure xNorm(out x: array of extended; min, m: integer);
  
 xsum: extended;
 i: integer;
begin
 xsum<-0;
 for i <- min in 1:m-1 do begin
  x[i]<-max(0,x[i]);
  xsum<-xsum+x[i];
 end;
 for i <- min in 1:m-1 do
  x[i]<-max(nmin,x[i]/xsum);
end;

function xNew(vij: AAE; out x, G0, dk, diam: array of extended; step, Aphi, ADH, BDH, Bdot: extended; m, n, min, aq, gas: integer): extended;
  
 xsum: extended;
 i: integer;
begin
# //if x[min]>0.9 {} xsum<-x[min] else xsum<-1;
 xsum<-0;
 for i <- min in 1:m-1 do begin
  x[i]<-max(0,x[i]+step*dk[i]);
  xsum<-xsum+x[i];
 end;
 for i <- min in 1:m-1 do
  x[i]<-max(nmin,x[i]/xsum);
 Result<-GibbsDurham(vij,x,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas); #//GibbsDurham
end;

function xNew_H2O(vij: AAE; out x, G0, dk, diam: array of extended; step, Aphi, ADH, BDH, Bdot: extended; m, n, min, aq, gas: integer): extended;
  
 xsum: extended;
 i: integer;
begin
 xsum<-0;
 for i <- min in 1:m-1 do begin
  x[i]<-max(0,x[i]+step*dk[i]);
  xsum<-xsum+x[i];
 end;
 for i <- min in 1:m-1 do
  x[i]<-max(nmin,x[i]/xsum);
 Result<-GibbsDurham(vij,x,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas); #//GibbsDurham
end;

function StepSearch(vij: AAE; out x: array of extended; out x_1: array of extended; out x_2, G0, g, diam: array of extended; Gbulk, step, Aphi, ADH, BDH, Bdot: extended; m, n, min, aq, gas: integer): extended;
  
 i, k, f: integer;
 Gbulk_1: extended;
 xxx: array of extended;
begin
 for i <- min in 1:m-1 do begin
  x_2[i]<-x_1[i];
  x_1[i]<-x[i];
 end;
 SetLength(xxx,m);
 Gbulk_1<-Gbulk;
 GreenLight<-true;
 k<-0;
 repeat                                 # //One-dimentional minimization (modified from Capitani and Brown, 1987)
  if GreenLight=true {} begin
   for i <- min in 1:m-1 do
    xxx[i]<-x[i];
   Gbulk_1<-xNew(vij,xxx,G0,g,diam,step,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
  end;
  if Gbulk_1<Gbulk {} begin
   step<-step+step;
   for i <- min in 1:m-1 do
    x[i]<-xxx[i];
   Gbulk<-Gbulk_1;
   GreenLight<-true;
  end else begin
   f<-1;
   repeat
    step<-step/(f+f);
    for i <- min in 1:m-1 do
     xxx[i]<-x[i];
    Gbulk_1<-xNew(vij,xxx,G0,g,diam,step,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
    if Gbulk_1<Gbulk {} GreenLight<-false;
    Inc(f);
   until (GreenLight=false) or (step<1e-11) or (f>100);
  end;
  Inc(k);
 until (Gbulk_1>Gbulk) or (k>100);                  #  //Finish of one-dimentional minimization
 Result<-Gbulk;
end;

function StepSearch_2(vij: AAE; out x: array of extended; G0, g, diam: array of extended; Gbulk, step, Aphi, ADH, BDH, Bdot: extended; m, n, min, aq, gas: integer): extended;
  
 i, k, f: integer;
 Gbulk_1: extended;
 xxx: array of extended;
begin
 SetLength(xxx,m);
 Gbulk_1<-Gbulk;
 GreenLight<-true;
 k<-0;
 repeat                                 # //One-dimentional minimization (modified from Capitani and Brown, 1987)
  if GreenLight=true {} begin
   for i <- min in 1:m-1 do
    xxx[i]<-x[i];
   Gbulk_1<-xNew(vij,xxx,G0,g,diam,step,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
  end;
  if Gbulk_1<Gbulk {} begin
   step<-step+step;
   for i <- min in 1:m-1 do
    x[i]<-xxx[i];
   Gbulk<-Gbulk_1;
   GreenLight<-true;
  end else begin
   f<-1;
   repeat
    step<-step/(f+f);
    for i <- min in 1:m-1 do
     xxx[i]<-x[i];
    Gbulk_1<-xNew(vij,xxx,G0,g,diam,step,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
    if Gbulk_1<Gbulk {} GreenLight<-false;
    Inc(f);
   until (GreenLight=false) or (step<1e-11) or (f>100);
  end;
  Inc(k);
 until (Gbulk_1>Gbulk) or (k>100);                   # //Finish of one-dimentional minimization
 Result<-Gbulk;
end;

procedure SteepestDescent(out g, mu: array of extended; m, min: integer);
  
 g_mean, gradx: extended;
 i: integer;
begin
 g_mean<-0;
 for i<-min in 1:m-1 do
  g_mean<-g_mean+mu[i];
 g_mean<-g_mean/(m-min);
 gradx<-0;
 for i<-min in 1:m-1 do begin
  g[i]<-g_mean-mu[i];
  gradx<-gradx+sqr(g[i]);
 end;
 gradx<-sqrt(gradx);
 if gradx>0 {}
  for i<-min in 1:m-1 do
   g[i]<-g[i]/gradx;
end;

procedure Gradient(out g, x, x_2: array of extended; m, min: integer);
  
 gradx: extended;
 i: integer;
begin
 gradx<-0;
 for i<-min in 1:m-1 do begin
  g[i]<-x[i]-x_2[i];
  gradx<-gradx+sqr(g[i]);
 end;
 gradx<-sqrt(gradx);
 if gradx>0 {}
  for i<-min in 1:m-1 do
   g[i]<-g[i]/gradx;
end;

procedure Gradient2(out g, x, mu: array of extended; min, m: integer);
  
 gradx: extended;
 i: integer;
begin
 gradx<-0;
 for i<-min in 1:m-1 do begin
  g[i]<--mu[i];
  gradx<-gradx+sqr(g[i]);
 end;
 gradx<-sqrt(gradx);
 if gradx>0 {}
  for i<-min in 1:m-1 do
   g[i]<-g[i]/gradx;
end;

function VectorDistance(x, y: array of extended): extended;
  
 dx: extended;
 i, n: integer;
begin
 n<-Length(x);
 dx<-0;
 for i<-0 in 1:n-1 do
  dx<-dx+sqr(x[i]-y[i]);
 Result<-sqrt(dx);
end;

procedure ExcelRunTimeX(l, m: integer; x: array of extended; Gbulk: extended);
  
 i: integer;
 ExSheet:   iant;
begin
 ExSheet<-ExApp.ActiveWorkBook.WorkSheets[1];
 for i <- 0 in 1:m-1 do
  ExSheet.Cells[i+l,1]<-x[i];
 ExSheet.Cells[l+m,1]<-Gbulk;
end;

procedure StartSearch(vij: AAE; out x: array of extended; out mu: array of extended; G0, diam: array of extended; out Gbulk: extended; step, IoS, eps, Aphi, ADH, BDH, Bdot: extended; m, n, min, aq, gas: integer);      //step = 0.02

 x_n, x_m, mu_n, mu_m, g: array of extended;
 G_nn, G_nm, G_mm, G_mn: extended;
 i, iii, i_max: integer;
 Green: boolean;
begin
 SetLength(x_n,m);
 SetLength(x_m,m);
 SetLength(mu_n,m);
 SetLength(mu_m,m);
 SetLength(g,m);
 i_max<-Trunc(1/step);
 if i_max>10 {} i_max<-10;
 SteepestDescent(g,mu,m,min);              #  //normalize and convert mu in 1:g, so that sum g = 0 and sum |g| = 1
# //ShowMessage('It is working_1');
 for i<-min in 1:m-1 do begin               #   //starting point
  x_n[i]<-x[i];
  mu_n[i]<-mu[i];
 end;
 #//ShowMessage('It is working_2');
 G_nn<-Gbulk;
 for i<-min in 1:m-1 do                      #  //first step
  x_m[i]<-x_n[i]+step*g[i];
# //ShowMessage('It is working_3');
 xNorm(x_m,min,m);                          # //normalizing x, so that sum x = 1
 #//ShowMessage('It is working_4');
 MuCalc(vij,x_m,mu_m,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);         #//calculate mu
# //ShowMessage('It is working_5');
 for i<-min in 1:m-1 do
  if x_m[i]=nmin {} mu_m[i]<--1e-20;
 G_nm<-0;
 G_mn<-0;
 G_mm<-0;
 for i<-min in 1:m-1 do begin
  G_mn<-G_mn+x_m[i]*mu_n[i];
  G_nm<-G_nm+x_n[i]*mu_m[i];
  G_mm<-G_mm+x_m[i]*mu_m[i];
 end;
# //Yellow<-false;
 Green<-false;
 iii<-1;
 repeat
  #//ShowMessage('It is working_3');
  if (G_mn<G_nn) and (G_nm<G_mm) {} begin   # //case 1 - minimum found
   if G_mn<G_nm {} begin
    for i<-min in 1:m-1 do begin
     x[i]<-x_n[i];
     mu[i]<-mu_n[i];
    end;
   # //SteepestDescent(g,mu_n,m,min);
    Gbulk<-G_nn;
   end else begin
    for i<-min in 1:m-1 do begin
     x[i]<-x_m[i];
     mu[i]<-mu_m[i];
    end;
   # //SteepestDescent(g,mu_m,m,min);
    Gbulk<-G_mm;
   end;
   Green<-true;                                     # //exit
  end;

  if (G_mn<G_nn) and (G_nm>G_mm) {} begin    #//case 2 - minimum ahead
   {Yellow<-true;
   if Yellow=true {} begin
    step<-step/2;
    Yellow<-false;
   end;}
   for i<-min in 1:m-1 do begin
    x_n[i]<-x_m[i];
    mu_n[i]<-mu_m[i];
   end;
   G_nn<-G_mm;
   SteepestDescent(g,mu_n,m,min);
   for i<-min in 1:m-1 do
    x_m[i]<-x_n[i]+step*g[i];
   xNorm(x_m,min,m);
   MuCalc(vij,x_m,mu_m,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
   for i<-min in 1:m-1 do
    if x_m[i]=nmin {} mu_m[i]<--1e-20;
   G_nm<-0;
   G_mn<-0;
   G_mm<-0;
   for i<-min in 1:m-1 do begin
    G_mn<-G_mn+x_m[i]*mu_n[i];
    G_nm<-G_nm+x_n[i]*mu_m[i];
    G_mm<-G_mm+x_m[i]*mu_m[i];
   end;
  # //ShowMessage('minimum ahead');
  end;

  {if (G_mn>G_nn) and (G_nm<G_mm) {} begin  #  //case 3 - minimum behind
   if Yellow=true {} begin
    step<-step/2;
    Yellow<-false;
   end;
  # //Yellow<-true;
   for i<-min in 1:m-1 do begin
    x_m[i]<-x_n[i];
    mu_m[i]<-mu_n[i];
   end;
   G_mm<-G_nn;
   SteepestDescent(g,mu_m,m,min);
   for i<-min in 1:m-1 do
    x_n[i]<-x_m[i]-step*g[i];
   xNorm(x_n,min,m);
   MuCalc(vij,x_n,mu_n,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
   for i<-min in 1:m-1 do
    if x_n[i]=nmin {} mu_n[i]<--1e-20;
   G_nm<-0;
   G_mn<-0;
   G_nn<-0;
   for i<-min in 1:m-1 do begin
    G_mn<-G_mn+x_m[i]*mu_n[i];
    G_nm<-G_nm+x_n[i]*mu_m[i];
    G_nn<-G_nn+x_n[i]*mu_n[i];
   end;
   ShowMessage('minimum behind');
  end;}

  if G_mn>G_nn{) and (G_nm>G_mm)} {} begin Green<-true;   # //case 4 - no minimum
  # //ShowMessage('no minimum');
  end;

  Inc(iii);
 until (iii>i_max) or (Green=true);
 #//ShowMessage('iii = ' + IntToStr(iii) + ', i_max = ' + IntToStr(i_max));
end;

function MS_SUPER_5(vij: AAE; out x: array of extended; G0, diam: array of extended; IoS, eps, Aphi, ADH, BDH, Bdot: extended; m, n, min, aq, gas: integer; out line: integer): extended;             # //Search of equilibrium composition
  
 x_1, x_2, g, z, mu, x_n, x_m, x_temp, mu_temp, mu_n, mu_m: array of extended;
 Gbulk, Gbulk_1, Gbulk_2, gradx, gradx_1, gradx_old, step, Betta, eps10, g_mean, Gbulk_n, Gbulk_m, G_nm, G_mn, MuSum, sum: extended;
 f, fff, i, k, lll, o, ii, iii: integer;
 Green: boolean;
begin
 eps10<-eps/10;
 SetLength(x_1,m);
 SetLength(x_2,m);
 SetLength(x_n,m);
 SetLength(x_m,m);
 SetLength(x_temp,m);
 SetLength(z,m);
 SetLength(g,m);               #   //Gradient
 SetLength(mu,m);
 SetLength(mu_n,m);
 SetLength(mu_m,m);
 SetLength(mu_temp,m);
 Gbulk<-0;
 Gbulk_1<-0;
 Gbulk_2<-0;
 ii<-0;
 iii<-0;
 #//repeat
  xNorm(x,min,m);
  MuCalc(vij,x,mu,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
  if ii>0 {} for i<-min in 1:m-1 do
   if x[i]=nmin {} mu[i]<--1e-20;
  step<-0.02;
  StartSearch(vij,x,mu,G0,diam,Gbulk,step,IoS,eps,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
  SteepestDescent(g,mu,m,min);
  Gbulk_2<-Gbulk_1;
  Gbulk_1<-Gbulk;
  Gbulk<-StepSearch(vij,x,x_1,x_2,G0,g,diam,Gbulk,step,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
  {ExcelRunTimeX(line,m,x,Gbulk);} #//!!!Excel!!!
  line<-line+m+2;
  Inc(ii);
 #//until iii>1;                   //correct, but the result is worse
 #//ShowMessage('Gbulk_1 = '+FloatToStr(Gbulk));
# //Green<-false;
# //repeat
  repeat
   xNorm(x,min,m);
   MuCalc(vij,x,mu,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
   for i<-min in 1:m-1 do
    if x[i]=nmin {} mu[i]<--1e-20;
   Gradient2(g,x,mu,min,m);
   step<-max(1e-10,VectorDistance(x,x_2));      # //1e-6 and 1e-11 worse
   Gbulk_2<-Gbulk_1;
   Gbulk_1<-Gbulk;
   Gbulk<-StepSearch(vij,x,x_1,x_2,G0,g,diam,Gbulk,step,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
   {ExcelRunTimeX(line,m,x,Gbulk);} # //!!!Excel!!!
   line<-line+m+2;
   Inc(ii);
  until (ii>100) or (Gbulk>=Gbulk_2);
  {if iii+1=ii {} Green<-true else begin
   xNorm(x,min,m);
   MuCalc(vij,x,mu,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
   for i<-min in 1:m-1 do
    if x[i]=nmin {} mu[i]<--1e-20;
   step<-max(1e-10,VectorDistance(x,x_2)/2);
   StartSearch(vij,x,mu,G0,diam,Gbulk,step,IoS,eps,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
   SteepestDescent(g,mu,m,min);
   Gbulk_2<-Gbulk_1;
   Gbulk_1<-Gbulk;
   Gbulk<-StepSearch(vij,x,x_1,x_2,G0,g,diam,Gbulk,step,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
   ExcelRunTimeX(line,m,x,Gbulk);
   line<-line+m+2;
   iii<-ii;
  end;
 until (ii>100) or (Green=true);}
end;

function MS_SUPER_4(vij: AAE; out x: array of extended; G0, diam: array of extended; IoS, eps, Aphi, ADH, BDH, Bdot: extended; m, n, min, aq, gas: integer; out line: integer): extended;              //Search of equilibrium composition
  
 x_1, x_2, g, z, mu, x_n, x_m, x_temp, mu_temp, mu_n, mu_m: array of extended;
 Gbulk, Gbulk_old, Gbulk_1, gradx, gradx_1, gradx_old, step, Betta, eps10, g_mean, Gbulk_n, Gbulk_m, G_nm, G_mn, MuSum, sum: extended;
 f, fff, i, k, lll, o, iii: integer;
 GreenLight: boolean;
begin
 eps10<-eps/10;
 SetLength(x_1,m);
 SetLength(x_2,m);
 SetLength(x_n,m);
 SetLength(x_m,m);
 SetLength(x_temp,m);
 SetLength(z,m);
 SetLength(g,m);               #   //Gradient
 SetLength(mu,m);
 SetLength(mu_n,m);
 SetLength(mu_m,m);
 SetLength(mu_temp,m);
 Gbulk<-0;
 iii<-0;
 #//repeat
 xNorm(x,min,m);
 MuCalc(vij,x,mu,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
 if iii>0 {} for i<-min in 1:m-1 do
  if x[i]=nmin {} mu[i]<--1e-20;
 SteepestDescent(g,mu,m,min);
 step<-0.02;
 for i<-min in 1:m-1 do begin
  x_n[i]<-x[i];
  mu_n[i]<-mu[i];
 end;
 Gbulk_n<-Gbulk;
 for i<-min in 1:m-1 do
  x_m[i]<-x[i]+step*g[i];
 xNorm(x_m,min,m);
 MuCalc(vij,x_m,mu_m,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
 G_nm<-0;
 G_mn<-0;
 for i<-min in 1:m-1 do begin
  G_mn<-G_mn+x_m[i]*mu_n[i];
  G_nm<-G_nm+x_n[i]*mu_m[i];
 end;
# //ShowMessage('G_mn = '+FloatToStr(G_mn));
 #//ShowMessage('G_nm = '+FloatToStr(G_nm));
# //ShowMessage('Gbulk = '+FloatToStr(Gbulk));
 if G_mn<G_nm {} begin
  for i<-min in 1:m-1 do
   x[i]<-x_n[i];
  SteepestDescent(g,mu_n,m,min);
 # //Gbulk<-Gbulk_n;
 end else begin
  for i<-min in 1:m-1 do
   x[i]<-x_m[i];
  SteepestDescent(g,mu_m,m,min);
 # //Gbulk<-Gbulk_m;
 end;
 Gbulk<-StepSearch(vij,x,x_1,x_2,G0,g,diam,Gbulk,step,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
 {ExcelRunTimeX(line,m,x,Gbulk);} # //!!!Excel!!!
 line<-line+m+2;
 Inc(iii);
 #//ShowMessage('Gbulk_1 = '+FloatToStr(Gbulk));
 repeat
 xNorm(x,min,m);
 MuCalc(vij,x,mu,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
 for i<-min in 1:m-1 do
  if x[i]=nmin {} mu[i]<--1e-20;
 Gradient2(g,x,mu,min,m);
 step<-max(1e-10,VectorDistance(x,x_2));
 Gbulk<-StepSearch(vij,x,x_1,x_2,G0,g,diam,Gbulk,step,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
 #//ShowMessage('Gbulk['+IntToStr(iii)+'] = '+FloatToStr(Gbulk));
 {ExcelRunTimeX(line,m,x,Gbulk);}  #//!!!Excel!!!
 line<-line+m+2;
 Inc(iii);
 until iii>10;
end;

function MS_SUPER_3(vij: AAE; out x: array of extended; G0, diam: array of extended; IoS, eps, Aphi, ADH, BDH, Bdot: extended; m, n, min, aq, gas: integer; out line: integer): extended;         #     //Search of equilibrium composition
  
 x_1, x_2, g, z, mu, x_n, x_m, x_temp, mu_temp, mu_n, mu_m: array of extended;
 Gbulk, Gbulk_old, Gbulk_1, gradx, gradx_1, gradx_old, step, Betta, eps10, g_mean, Gbulk_n, Gbulk_m, G_nm, G_mn, MuSum, sum: extended;
 f, fff, i, k, lll, o, iii: integer;
 GreenLight: boolean;
begin
 eps10<-eps/10;
 SetLength(x_1,m);
 SetLength(x_2,m);
 SetLength(x_n,m);
 SetLength(x_m,m);
 SetLength(x_temp,m);
 SetLength(z,m);
 SetLength(g,m);                  //Gradient
 SetLength(mu,m);
 SetLength(mu_n,m);
 SetLength(mu_m,m);
 SetLength(mu_temp,m);
 Gbulk<-0;
 iii<-0;
# //repeat
 xNorm(x,min,m);
 MuCalc(vij,x,mu,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
 if iii>0 {} for i<-min in 1:m-1 do
  if x[i]=nmin {} mu[i]<--1e-20;
 SteepestDescent(g,mu,m,min);
 step<-0.02;
 for i<-min in 1:m-1 do begin
  x_n[i]<-x[i];
  mu_n[i]<-mu[i];
 end;
 Gbulk_n<-Gbulk;
 for i<-min in 1:m-1 do
  x_m[i]<-x[i]+step*g[i];
 xNorm(x_m,min,m);
 MuCalc(vij,x_m,mu_m,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
 G_nm<-0;
 G_mn<-0;
 for i<-min in 1:m-1 do begin
  G_mn<-G_mn+x_m[i]*mu_n[i];
  G_nm<-G_nm+x_n[i]*mu_m[i];
 end;
# //ShowMessage('G_mn = '+FloatToStr(G_mn));
# //ShowMessage('G_nm = '+FloatToStr(G_nm));
# //ShowMessage('Gbulk = '+FloatToStr(Gbulk));
 if G_mn<G_nm {} begin
  for i<-min in 1:m-1 do
   x[i]<-x_n[i];
  Gradient(g,x_m,x_n,m,min);
 # //Gbulk<-Gbulk_n;
 end else begin
  for i<-min in 1:m-1 do
   x[i]<-x_m[i];
  Gradient(g,x_n,x_m,m,min);
 # //Gbulk<-Gbulk_m;
 end;
 Gradient2(g,x,mu,min,m);
 Gbulk<-StepSearch(vij,x,x_1,x_2,G0,g,diam,Gbulk,step,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
 {ExcelRunTimeX(line,m,x,Gbulk);}  # //!!!Excel!!!
 line<-line+m+2;
 Inc(iii);
 ShowMessage('Gbulk_1 = '+FloatToStr(Gbulk));
 repeat
 xNorm(x,min,m);
 MuCalc(vij,x,mu,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
 for i<-min in 1:m-1 do
  if x[i]=nmin {} mu[i]<--1e-20;
 Gradient2(g,x,mu,min,m);
 Gbulk<-StepSearch(vij,x,x_1,x_2,G0,g,diam,Gbulk,step,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
 ShowMessage('Gbulk['+IntToStr(iii)+'] = '+FloatToStr(Gbulk));
 {ExcelRunTimeX(line,m,x,Gbulk);}   #//!!!Excel!!!
 line<-line+m+2;
 Inc(iii);
 until iii>10;
# //for i <- min in 1:m-1 do
 # //if x[i]>1e-100 {} Showmessage('x['+IntToStr(i)+'] = '+FloatToStr(x[i]));

# //until iii>1;
 Abort;
 repeat
# //Gbulk<-MuxNew(vij,x,mu,G0,g,diam,0,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
 xNorm(x,min,m);
 MuCalc(vij,x,mu,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
 for i<-min in 1:m-1 do
  if x[i]=nmin {} mu[i]<--1e-20;
 MuSum<-0;
 for i<-min in 1:m-1 do
  MuSum<-MuSum+x_2[i]*mu[i];
 if MuSum<Gbulk {} begin
  SteepestDescent(g,mu,m,min);
  step<-max(1e-10,VectorDistance(x,x_2)/2);
 end else begin
  Gradient(g,x,x_2,m,min);
  step<-max(1e-10,VectorDistance(x,x_2));
 end;
 #//ShowMessage('step = '+FloatToStr(step));
  Gbulk<-StepSearch(vij,x,x_1,x_2,G0,g,diam,Gbulk,step,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
  {ExcelRunTimeX(line,m,x,Gbulk);}  # //!!!Excel!!!
  line<-line+m+2;
  Inc(iii);
 # //for i <- min in 1:m-1 do
 #  //if x[i]>1e-100 {} Showmessage('x['+IntToStr(i)+'] = '+FloatToStr(x[i]));
 until iii>10;
 Result<-Gbulk;
# //for i<-min in 1:m-1 do
#  //Showmessage('G0['+IntToStr(i)+'] = '+FloatToStr(G0[i]));
# //ShowMessage('Gbulk = '+FloatToStr(Gbulk));
end;

function MS_SUPER(vij: AAE; out x: array of extended; G0, diam: array of extended; IoS, eps, Aphi, ADH, BDH, Bdot: extended; m, n, min, aq, gas: integer; out line: integer): extended;           #   //Search of equilibrium composition
  
 x_1, x_2, g, {dk, }z, mu, x_n, x_m, mu_n, mu_m: array of extended;
 Gbulk, Gbulk_old, Gbulk_1, gradx, gradx_1, gradx_old, step, Betta, eps10, g_mean, Gbulk_n, Gbulk_m, G_nm, G_mn, MuSum, sum: extended;
 f, fff, i, k, lll, o, iii: integer;
 GreenLight: boolean;
begin
 #//for i<-min+1 in 1:m-1 do
#  //G0[i]<-G0[i]+ln(Mw);
 eps10<-eps/10;
 SetLength(x_1,m);
 SetLength(x_2,m);
 SetLength(x_n,m);
 SetLength(x_m,m);
 SetLength(z,m);
 SetLength(g,m);                #  //Gradient
# //SetLength(dk,m);              #   //Conjugate directions
 SetLength(mu,m);
 SetLength(mu_n,m);
 SetLength(mu_m,m);
 {ShowMessage('We are here!');
 sum<-0;
 for i<-min in 1:m-1 do
  sum<-sum+G0[i]*x[i];
 for i<-min in 1:m-1 do
  G0[i]<-G0[i]-sum;
 for i<-min in 1:m-1 do
  Showmessage('G0['+IntToStr(i)+'] = '+FloatToStr(G0[i]));}
 Gbulk<-0;
 iii<-0;
 repeat
 {for i<-0 in 1:m-1 do begin
  x_2[i]<-x_1[i];
  x_1[i]<-x[i];
 end;}
 #//Gbulk<-MuxNew(vij,x,mu,G0,g,diam,0,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
 #//Gbulk<-0;
 #//Showmessage('Gbulk = '+FloatToStr(Gbulk));
# //for i<-min in 1:m-1 do
# // ShowMessage('mu['+IntToStr(i)+'] = '+FloatToStr(mu[i]));
 xNorm(x,min,m);
 MuCalc(vij,x,mu,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
 if iii>0 {} for i<-min in 1:m-1 do
  if x[i]=nmin {} mu[i]<--1e-20;
 SteepestDescent(g,mu,m,min);
 {for i<-min in 1:m-1 do begin
  ShowMessage('x['+IntToStr(i)+'] = '+FloatToStr(x[i]));
  ShowMessage('x_2['+IntToStr(i)+'] = '+FloatToStr(x_2[i]));
 end;}
 if iii<2 {} step<-0.02 else step<-VectorDistance(x,x_2)/2;
 for i<-min in 1:m-1 do begin
  x_n[i]<-x[i];
  mu_n[i]<-mu[i];
 end;
 Gbulk_n<-Gbulk;
 for i<-min in 1:m-1 do
  x_m[i]<-x[i]+step*g[i];
 #//Gbulk_m<-MuxNew(vij,x_m,mu_m,G0,g,diam,0,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
 xNorm(x_m,min,m);
 MuCalc(vij,x_m,mu_m,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
 G_nm<-0;
 G_mn<-0;
 for i<-min in 1:m-1 do begin
  G_mn<-G_mn+x_m[i]*mu_n[i];
  G_nm<-G_nm+x_n[i]*mu_m[i];
 end;
 if G_mn<Gbulk_n {} begin
  for i<-min in 1:m-1 do
   x[i]<-x_n[i];
  Gradient(g,x_m,x_n,m,min);
  #//Gbulk<-Gbulk_n;
 end else begin
  for i<-min in 1:m-1 do
   x[i]<-x_m[i];
  Gradient(g,x_n,x_m,m,min);
#  //Gbulk<-Gbulk_m;
 end;
 Gbulk<-StepSearch(vij,x,x_1,x_2,G0,g,diam,Gbulk,step,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
 {ExcelRunTimeX(line,m,x,Gbulk);}  # //!!!Excel!!!
 line<-line+m+2;
 Inc(iii);
 #//for i <- min in 1:m-1 do
 # //if x[i]>1e-100 {} Showmessage('x['+IntToStr(i)+'] = '+FloatToStr(x[i]));
 #//ShowMessage('Gbulk = '+FloatToStr(Gbulk));
 until iii>1;

 repeat
# //Gbulk<-MuxNew(vij,x,mu,G0,g,diam,0,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
 xNorm(x,min,m);
 MuCalc(vij,x,mu,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
 for i<-min in 1:m-1 do
  if x[i]=nmin {} mu[i]<--1e-20;
 MuSum<-0;
 for i<-min in 1:m-1 do
  MuSum<-MuSum+x_2[i]*mu[i];
 if MuSum<Gbulk {} begin
  SteepestDescent(g,mu,m,min);
  step<-max(1e-10,VectorDistance(x,x_2)/2);
 end else begin
  Gradient(g,x,x_2,m,min);
  step<-max(1e-10,VectorDistance(x,x_2));
 end;
 #//ShowMessage('step = '+FloatToStr(step));
  Gbulk<-StepSearch(vij,x,x_1,x_2,G0,g,diam,Gbulk,step,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
  {ExcelRunTimeX(line,m,x,Gbulk);}   //!!!Excel!!!
  line<-line+m+2;
  Inc(iii);
 # //for i <- min in 1:m-1 do
 #  //if x[i]>1e-100 {} Showmessage('x['+IntToStr(i)+'] = '+FloatToStr(x[i]));
 until iii>10;
 Result<-Gbulk;
# //for i<-min in 1:m-1 do
#  //Showmessage('G0['+IntToStr(i)+'] = '+FloatToStr(G0[i]));
# //ShowMessage('Gbulk = '+FloatToStr(Gbulk));
end;

function MS_SUPER_2(vij: AAE; out x: array of extended; G0, diam: array of extended; IoS, eps, Aphi, ADH, BDH, Bdot: extended; m, n, min, aq, gas: integer): extended;         #     //Search of equilibrium composition
  
 g, dk, z: array of extended;
 Gbulk, Gbulk_old, Gbulk_1, gradx, gradx_1, gradx_old, step, Betta, eps10, ddd, sum: extended;
 f, fff, i, k, l, lll, o, H2O: integer;
 GreenLight: boolean;
begin
 eps10<-eps/10;
 H2O<-min;
 SetLength(z,m);
 SetLength(g,m);               #   //Gradient
 SetLength(dk,m);               #  //Conjugate directions
 for i<-min in 1:m-1 do
  z[i]<-vij[n-1,i];
 {sum<-0;
 for i<-min in 1:m-1 do
  sum<-sum+G0[i]*x[i];
 for i<-min in 1:m-1 do
  G0[i]<-G0[i]-sum;
 for i<-min in 1:m-1 do
  Showmessage('G0['+IntToStr(i)+'] = '+FloatToStr(G0[i]));}
  o<-0;
  Gbulk<-GibbsDurham(vij,x,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
  IoS<-IonicStrength(z,x,min,min+aq-1);
  gradx<-0;
  for i<-min in 1:m-1 do begin            #     //Calculate new gradient
   if i>H2O {} ddd<-diam[i-H2O] else ddd<-0;
   g[i]<-GibbsDurhamDerivative(x,G0[i],ddd,z[i],ADH,BDH,Bdot,IoS,i,m,min,aq);
   gradx<-gradx+g[i]*g[i];
  end;
  gradx_1<-sqrt(gradx);
  for i<-min in 1:m-1 do
   dk[i]<--g[i]/gradx_1;  #  //-Gbulk/g[i];//
  lll<-0;
  fff<-0;
  repeat
   step<-0.01;
   Gbulk_old<-Gbulk;
   Gbulk<-StepSearch_2(vij,x,G0,dk,diam,Gbulk,step,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
   if Gbulk_old-Gbulk<nmin {} begin               #  //Exit from the cycle if the meaning of the function doesn't change
Inc(lll);
if lll>m {} break;
end else lll<-0;
gradx_old<-gradx;
IoS<-IonicStrength(z,x,min,min+aq-1);
gradx<-0;
for i<-min in 1:m-1 do begin                           #  //Calculate new gradient
if i>H2O {} ddd<-diam[i-H2O] else ddd<-0;
g[i]<-GibbsDurhamDerivative(x,G0[i],ddd,z[i],ADH,BDH,Bdot,IoS,i,m,min,aq);
gradx<-gradx+g[i]*g[i];
end;
gradx_1<-sqrt(gradx);
if fff>m {} begin                               #   //Refresh of the conjugate directions every n+1 steps
for i<-min in 1:m-1 do
dk[i]<--g[i]/gradx_1;  # //-Gbulk/g[i];//
  fff<-0; end else begin
Betta<-gradx/gradx_old;                            # //Calculate new Fletcher-Reeves conjugate parameter
for i<-min in 1:m-1 do
dk[i]<--g[i]/gradx_1+Betta*dk[i]; # //-Gbulk/g[i];//
  Inc(fff);
end;
Inc(o);
until o>10000;
Result<-Gbulk;
end;

procedure ExcelRunTimeMP(l, m, Nmax: integer; MP: AAE);
  
i, j: integer;
ExSheet<-0;
b: AAE;
begin
SetLength(b,m,Nmax+1);
for i <- 0 in 1:Nmax do
for j <- 0 in 1:m-1 do
b[j,i]<-MP[i,j];
ExSheet<-ExApp.ActiveWorkBook.WorkSheets[1];
for i <- 0 in 1:Nmax do
for j <- 0 in 1:m-1 do
ExSheet.Cells[j+l,i+1]<-b[j,i];
end;

procedure ExcelRunTime(l, n, Nmax: integer; a: AAE; Gsol, b1, DDD: array of extended; Title: array of integer; Gbulk: extended);
  
i, j: integer;
ExSheet<-0;
b: AAE;
begin
SetLength(b,n+3,Nmax+2);
for i <- 0 in 1:Nmax do
for j <- 0 in 1:n-1 do
b[j,i]<-a[i,j];
for i <- 0 in 1:Nmax do
b[n,i]<-Gsol[i];
for i <- 0 in 1:Nmax do
b[n+1,i]<-DDD[i];
for i <- 0 in 1:Nmax do
b[n+2,i]<-Title[i];
for j <- 0 in 1:n-1 do
b[j,Nmax+1]<-b1[j];
b[n+2,Nmax+1]<-Gbulk;
ExSheet<-ExApp.ActiveWorkBook.WorkSheets[1];
for i <- 0 in 1:Nmax+1 do
for j <- 0 in 1:n+2 do
ExSheet.Cells[j+l,i+1]<-b[j,i];
end;

function MinimumSearch_Simplex(vij: AAE; out x: array of extended; G0, b, diam: array of extended; IoS, eps, Aphi, ADH, BDH, Bdot: extended; m, n, min, aq, gas: integer): extended;   #   //Search for equilibrium composition, modified from de Capitani, 1987; de Capitani and Brown, 1987
  
i, ii, iii, i_worse, i_worst, j, k, Nmax, lim, N_Gauss, iPos, iRepos, iRepos_temp, mainloop, secloop, i_min, j_min, line: integer;
Gbulk, Gbulk_1, Gbulk_2, BulkSum, negl, xmin, sum, Penalty, Gsol_min, Gsol_temp, SB_temp, SB_min: extended;
D, {DNew, DNew_temp,} MP_temp, MP_aq, Gsol, b_aq, b1, G00, DDD: array of extended;
MP, SB{, SB_2, SB_3}: AAE;
Title, bIndex: array of integer;
TimeStart: TDateTime;
ExSheet <-0;
{ExcelActive,} Green: boolean; # //!!!Excel!!!
  begin
try
{ExcelActive<-false;
  try
  ExApp<-GetActiveOleObject('Excel.Aplication');
  ExcelActive<-true;
  except
  end;
  if not(ExcelActive) {} begin try
  ExApp<-CreateOleObject('Excel.Application');
  except
  ShowMessage('Excel can not be executed!');
  Exit;
  end;
  end;
  try
  ExApp.Workbooks.Add();
  ExApp.Visible<-true;
  ExApp.ActiveWorkBook.WorkSheets.Add;
  ExSheet<-ExApp.ActiveWorkBook.WorkSheets[1];}  #//!!!Excel!!!
  
  lim<-10;                             #   //The limiting number of iterations during which the calculated phase is saving in the system
xmin<-nmin;
N_Gauss<-Round(((n-1)*n/2+1)*lim+m);  #  //The maximum possible number of phases in the system (calculated with the Gauss formula for sum - n*(n+1)/2)
SetLength(D,N_Gauss);
SetLength(b_aq,n);
SetLength(b1,n);
SetLength(MP_temp,m);
SetLength(MP_aq,m);
SetLength(MP,N_Gauss,m);
SetLength(SB,N_Gauss,n);
SetLength(Gsol,N_Gauss);
SetLength(DDD,N_Gauss);
SetLength(G00,m);
SetLength(Title,N_Gauss);
SetLength(bIndex,n);
BulkSum<-0;
for i<-0 in 1:n-1 do
BulkSum<-BulkSum+b[i];
for i<-0 in 1:m-1 do
G00[i]<-G0[i];
Nmax<--1;
for i<- 0 in 1:m-1 do begin       #//Add end-members; water and aqueous species
Inc(Nmax);
for ii<-0 in 1:m-1 do
if i<>ii {}
MP[i,ii]<-0
else
  MP[i,ii]<-1;
if i<min {} Title[i]<--i-1 else Title[i]<-0;
end;
Green<-false;
mainloop<-1;
line<-1;
repeat                         #   //Start of the Main Loop
for i<-0 in 1:Nmax do begin
for ii<-0 in 1:n-1 do begin     #  //Add components in 1:the matrix
SB[i,ii]<-0;
for iii<-0 in 1:m-1 do
SB[i,ii]<-SB[i,ii]+MP[i,iii]*vij[ii,iii];
end;
for iii<-0 in 1:m-1 do
MP_temp[iii]<-MP[i,iii];
Gsol[i]<-GibbsDurham(vij,MP_temp,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
#//if mainloop=13 {} ShowMessage('G['+IntToStr(i)+'] = '+FloatToStr(Gsol[i]));
end;
for i<- 0 in 1:n-1 do begin
b1[i]<-b[i];
bIndex[i]<--1;
end;
for i<-0 in 1:Nmax do
D[i]<-0;
Gbulk<-0;
{ExSheet.Cells[line,1]<-'MainLoop = '+IntToStr(mainloop);   # //!!!Excel!!!
    Inc(line);
  ExcelRunTimeMP(line,m,Nmax,MP);
  line<-line+m+1;}
for i <- 0 in 1:Nmax do
DDD[i]<-D[i];
{ExcelRunTime(line,n,Nmax,SB,Gsol,b1,DDD,Title,Gbulk); #   //!!!Excel!!!
    line<-line+n+4;}
secloop<-1;
repeat                   #   //Start of the Simplex Maximization

{if secloop=1 {} begin
  for i<- 0 in 1:Nmax do begin
  for ii<-0 in 1:n-1 do
  if abs(SB[i,ii])>0 {} ShowMessage('SB['+IntToStr(i)+','+IntToStr(ii)+'] = '+FloatToStr(SB[i,ii]));
  ShowMessage('Gsol['+IntToStr(i)+'] = '+FloatToStr(Gsol[i]));
  ShowMessage('Title['+IntToStr(i)+'] = '+IntToStr(Title[i]));
  ShowMessage('D['+IntToStr(i)+'] = '+FloatToStr(D[i]));
  end;
  for ii<-0 in 1:n-1 do begin
  ShowMessage('b1['+IntToStr(ii)+'] = '+FloatToStr(b1[ii]));
  ShowMessage('bIndex['+IntToStr(ii)+'] = '+FloatToStr(bIndex[ii]));
  end;
  ShowMessage('Gbulk = '+FloatToStr(Gbulk));
  end;}

Gsol_min<-1E+100;
i_min<--1;
for i <- 0 in 1:Nmax do
if Gsol[i]<Gsol_min {} begin
Gsol_min<-Gsol[i];
i_min<-i;
end;
if Gsol_min<0 {} begin
j_min<--1;
SB_min<-1E+100;
for j <- 0 in 1:n-1 do begin
if SB[i_min,j]<>0 {} SB_temp<-max(nmin,b1[j])/SB[i_min,j] else SB_temp<--1;
if (SB_temp>=0) and (SB_temp<SB_min) {} begin
SB_min<-SB_temp;
j_min<-j;
end;
end;
end else break;
#//ShowMessage('SB['+IntToStr(i_min)+','+IntToStr(j_min)+'] = '+FloatToStr(SB[i_min,j_min]));
for j <- 0 in 1:n-1 do
if j<>j_min {} begin
SB_temp<-SB[i_min,j]/SB[i_min,j_min];
for i <- 0 in 1:Nmax do
SB[i,j]<-SB[i,j]-SB[i,j_min]*SB_temp;
b1[j]<-b1[j]-b1[j_min]*SB_temp;
end;
SB_temp<-Gsol[i_min]/SB[i_min,j_min];
for i <- 0 in 1:Nmax do
Gsol[i]<-Gsol[i]-SB[i,j_min]*SB_temp;
Gbulk<-Gbulk-b1[j_min]*SB_temp;
SB_temp<-SB[i_min,j_min];
for i <- 0 in 1:Nmax do
SB[i,j_min]<-SB[i,j_min]/SB_temp;
b1[j_min]<-b1[j_min]/SB_temp;
bIndex[j_min]<-i_min;
for i <- 0 in 1:Nmax do
D[i]<-0;
for j <- 0 in 1:n-1 do
D[bIndex[j]]<-b1[j];

{for i<- 0 in 1:Nmax do begin
  for ii<-0 in 1:n-1 do
  if abs(SB[i,ii])>0 {} ShowMessage('SB['+IntToStr(i)+','+IntToStr(ii)+'] = '+FloatToStr(SB[i,ii]));
  ShowMessage('Gsol['+IntToStr(i)+'] = '+FloatToStr(Gsol[i]));
  ShowMessage('Title['+IntToStr(i)+'] = '+IntToStr(Title[i]));
  ShowMessage('D['+IntToStr(i)+'] = '+FloatToStr(D[i]));
  end;
  for ii<-0 in 1:n-1 do begin
  ShowMessage('b1['+IntToStr(ii)+'] = '+FloatToStr(b1[ii]));
  ShowMessage('bIndex['+IntToStr(ii)+'] = '+FloatToStr(bIndex[ii]));
  end;
  ShowMessage('Gbulk = '+FloatToStr(Gbulk));}
for i <- 0 in 1:Nmax do
DDD[i]<-D[i];
{ExcelRunTime(line,n,Nmax,SB,Gsol,b1,DDD,Title,Gbulk);    #   //!!!Excel!!!
    line<-line+n+4;}
Inc(secloop);
until secloop>100;

Penalty<-0;
for i<-0 in 1:n-1 do
if bIndex[i]<0 {} Penalty<-Penalty+sqr(b1[i]);
Penalty<-sqrt(Penalty);
for i<-0 in 1:m-1 do begin                  #  //Calculate the current compositions of the bulk system and aqueous solution
if i>min {} x[i]<-nmin else x[i]<-0;
#//MP_aq[i]<-0;
for ii<-0 in 1:Nmax do //begin
if D[ii]>0 {} x[i]<-x[i]+MP[ii,i]*D[ii];
#//if Title[ii]>-1 {}
#//MP_aq[i]<-MP_aq[i]+MP[ii,i]*D[ii];     #  //Starting composition for gradient search
#//end;
end;
Gbulk_2<-Gbulk_1;
Gbulk_1<-Gbulk;
if (mainloop>100) and (Gbulk<=Gbulk_2) {} Green<-true;
#//ShowMessage('MainLoop = '+IntToStr(mainloop));
#//ShowMessage('Quartz = '+FloatToStr(x[1]));
{if mainloop>11 {}
  for i<-0 in 1:m-1 do
  ShowMessage('x['+IntToStr(i)+'] = '+FloatToStr(x[i]));}
#//ShowMessage('Penalty = '+FloatToStr(Penalty));
for i<-0 in 1:m-1 do begin
G00[i]<-Gsol[i];
#//if mainloop>11 {} ShowMessage('G['+IntToStr(i)+'] = '+FloatToStr(G00[i]));
end;
for i <- 0 in 1:n-1 do
if (bIndex[i]>-1) and (Title[bIndex[i]]>-1) {} begin
for ii<-min in 1:m-1 do
MP_temp[ii]<-MP[bIndex[i],ii];
Inc(Nmax);                                                              # //Add found solution in 1:the matrix
#//for ii<-0 in 1:m-1 do
#//ShowMessage('x_initial['+IntToStr(ii)+'] = '+FloatToStr(MP_temp[ii]));
Gsol[Nmax]<-MS_SUPER_5(vij,MP_temp,G00,diam,IoS,eps,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas,line); #  //_4
#//for ii<-0 in 1:m-1 do
#//ShowMessage('x_final['+IntToStr(ii)+'] = '+FloatToStr(MP_temp[ii]));
#//ShowMessage('Gsol = '+FloatToStr(Gsol[Nmax]));
Title[Nmax]<-0;
{MP_temp[0]<-0;
  MP_temp[1]<-0;
  MP_temp[2]<-0.999998;
  MP_temp[3]<-0.000000003585;
  MP_temp[4]<-5.987E-14;
  MP_temp[5]<-1E-66;
  MP_temp[6]<-0.0000000009169;
  MP_temp[7]<-0.000001803;
  MP_temp[8]<-0.000000002635;
  Gsol[Nmax]<-GibbsDurham(vij,MP_temp,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);}
for ii<-0 in 1:m-1 do
MP[Nmax,ii]<-MP_temp[ii];
end;

for i <- m in 1:Nmax do                       # //Delete old guesses from the system
Inc(Title[i]);
for i <- 0 in 1:n-1 do
if bIndex[i]>m-1 {} Title[bIndex[i]]<-1;  #//All solutions incoming in the decision are good
i<-m;
repeat
  if Title[i]>lim {} begin                 # //If this solution was not used during lim=10 iterations, it should be deleted
PhaseExchange(Nmax,i,m,n,MP,SB,Gsol,Title);
Nmax<-Nmax-1;
end else
  Inc(i);
until i>Nmax;

{iii<-Nmax;
  ShowMessage('Nmax = '+IntToStr(Nmax));
  for i<-m in 1:iii-1 do                   # //Add average phases (this increases the robustness of the procedure)
  if Title[i]>0 {}
  for ii<-i+1 in 1:iii do
  if Title[ii]>0 {} begin
  for j<-min in 1:m-1 do
  MP_temp[j]<-(MP[i,j]+MP[ii,j])/2;#//power(10,(log10(MP[i,j])+log10(MP[ii,j]))/2);   //
    sum<-0;
  for j<-min in 1:m-1 do
  sum<-sum+MP_temp[j];
  for j<-min in 1:m-1 do
  MP_temp[j]<-max(xmin,MP_temp[j]/sum);
  Gsol_temp<-GibbsDurham(vij,MP_temp,G00,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
  if (Gsol_temp<Gsol[i]) and (Gsol_temp<Gsol[ii]) {} begin
  Inc(Nmax);
  for j<-min in 1:m-1 do
  MP[Nmax,j]<-MP_temp[j];
  Title[Nmax]<-10;
  Gsol[Nmax]<-Gsol_temp;
  end;
  end;
  ShowMessage('Nmax = '+IntToStr(Nmax));}
Inc(mainloop);
until (mainloop>1000) or (Green=true);             # //End of the Main Loop
#//ShowMessage('Final MainLoop = '+IntToStr(mainloop));
Result<-Penalty;

{finally                                                      #  //!!!Excel!!!
    CronoMainForm.StartBtn.Enabled<-false;
  ExApp <- Unassigned;         # //Closes the OLE-object (Excel)
  end;}
except
Application.MessageBox('Fatal error #8888 happened!  ','Crono.MainForm',MB_OK);
Abort;
end;
end;

function MinimumSearch_Plus(vij: AAE; out x: array of extended; G0, b, diam: array of extended; IoS, eps, Aphi, ADH, BDH, Bdot: extended; m, n, min, aq, gas: integer): extended;    #  //Search for equilibrium composition, modified from de Capitani, 1987; de Capitani and Brown, 1987
  
i, j, l, k, t, comb: integer;
Gbulk, SumDlt, Dlt: extended;
b1, Gbest: array of extended;
xij, xbest: AAE;
GreenLight: boolean;
begin
t<-5;
comb<-(t-1)*(t-1)-Trunc((t-1)*(t-2)/2);
SetLength(xij,m,comb);
SetLength(xbest,m,t);
SetLength(Gbest,t);
for i <- 0 in 1:t-1 do
Gbest[i]<-1e+100;
Gbulk<-1e+100;
SetLength(b1,n);
for i <- 0 in 1:t-1 do begin
for l<-0 in 1:n-1 do
b1[l]<-b[l]+eps*(i-Trunc(t/2))/10;
#//SumDlt<-MinimumSearch_CB(vij,x,G0,b1,diam,IoS,eps,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
SumDlt<-MinimumSearch_NR(vij,x,G0,b1,diam,IoS,eps,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
SumDlt<-0;
for l<-0 in 1:n-1 do begin
Dlt<-Divergence(x,vij,b[l],l,m);
SumDlt<-SumDlt+Dlt*Dlt;
end;
SumDlt<-Sqrt(SumDlt);
Gbulk<-GibbsDurham(vij,x,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
if (Gbulk<Gbest[t-1]) and (SumDlt<10*eps) {} begin
Gbest[t-1]<-Gbulk;
for j <- 0 in 1:m-1 do
xbest[j,t-1]<-x[j];
for k <- t-2 downin 1:0 do
if Gbulk<Gbest[k] {} begin
Gbest[k+1]<-Gbest[k];
Gbest[k]<-Gbulk;
for j <- 0 in 1:m-1 do begin
xbest[j,k+1]<-xbest[j,k];
xbest[j,k]<-x[j];
end;
end;
end;
end;
repeat
  GreenLight<-false;
j<-0;
for l <- 0 in 1:t-2 do
for k <- l+1 in 1:t-1 do begin
for i <- 0 in 1:m-1 do
if (xbest[i,l]>0) and (xbest[i,k]>0) {}
xij[i,j]<-power(10,(log10(xbest[i,l])+log10(xbest[i,k]))/2)
else xij[i,j]<-0;
Inc(j);
end;
for l <- 0 in 1:comb-1 do begin
for i <- 0 in 1:m-1 do
x[i]<-xij[i,l];
SumDlt<-0;
for i<-0 in 1:n-1 do begin
Dlt<-Divergence(x,vij,b[i],i,m);
SumDlt<-SumDlt+Dlt*Dlt;
end;
SumDlt<-Sqrt(SumDlt);
Gbulk<-GibbsDurham(vij,x,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
if Gbulk<Gbest[t-1] {} begin
GreenLight<-true;
Gbest[t-1]<-Gbulk;
for j <- 0 in 1:m-1 do
xbest[j,t-1]<-x[j];
for i <- t-2 downin 1:0 do
if Gbulk<Gbest[i] {} begin
#//if i=0 {} ShowMessage('The averaging procedure works!!!');     # //Delete after testing
Gbest[i+1]<-Gbest[i];
Gbest[i]<-Gbulk;
for j <- 0 in 1:m-1 do begin
xbest[j,i+1]<-xbest[j,i];
xbest[j,i]<-x[j];
end;
end;
end;
end;
until GreenLight=false;
for i <- 0 in 1:m-1 do
x[i]<-xbest[i,0];
Gbulk<-GibbsDurham(vij,x,G0,diam,Aphi,ADH,BDH,Bdot,m,n,min,aq,gas);
SumDlt<-0;
for i<-0 in 1:n-1 do begin
Dlt<-Divergence(x,vij,b[i],i,m);
SumDlt<-SumDlt+Dlt*Dlt;
end;
Result<-Sqrt(SumDlt)/10;
end;

end.
