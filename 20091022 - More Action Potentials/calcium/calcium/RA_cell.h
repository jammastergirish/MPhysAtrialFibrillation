#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Atria_Variables.h"

float I_tot_RA ()
{
   float IKf,IKs,Ik1,Ito,IbNa,Ip,INaCa,INa,ICaL,Isus;
   
   float Apa, Bpa;
   float pam, tpa;
   float Api, Bpi;
   float pim, tpi;
   float An, Bn;
   float nm, tn;
   float Ar, Br;
   float rm, tr;
   float s1m, ts1;
   float s2m, ts2;
   float s3m, ts3;
   float Am, Bm;
   float Ah, Bh;
   float hm, th1, th2;
   float Adl, Bdl;
   float dlm, tdl;
   float Afl, Bfl;
   float flm, tfl;


   Apa=9.0*exp(V/25.371);
   Bpa=1.3*exp(-V/13.026);
   pam=1./(1+exp(-(V+12)/5.9));
   tpa=1./(Apa+Bpa);
   Pa+=HT*(pam-Pa)/tpa;
      
   Api=100.*exp(-V/54.645);
   Bpi=656.*exp(V/106.157);
   pim=1./(1+exp((V+24.3921)/18.6603));
   tpi=1.0/(Api+Bpi);
   Pi+=HT*(pim-Pi)/tpi;

   IKf=0.25*3.5*Pa*Pi*(V-Ek);   // 0.35 - LA, 0.25 - RA


   An=1.66*exp(V/69.452);
   Bn=0.3*exp(-V/21.826);
   nm=1.842/(1+exp(-(V-18)/15.8));
   tn=1./(An+Bn)+0.06;
   n+=HT*(nm-n)/tn;

   IKs=2.5*n*(V-Ek);     



   Ik1=10.16*pow(Kc/(Kc+0.59),3)*(V-Ek-10)/(1+exp(1.393*(V-Ek+3.6)/RTONF)); 


   Ar=386.6*exp(V/12.0);
   Br=8.011*exp(-V/7.2);
   rm=1/(1+exp(-(V+12.0)/9.6));   
   tr=1./(Ar+Br)+0.0004;
   r1+=HT*(rm-r1)/tr;

   s1m=1./(1+exp((V+32)/8.4));
   ts1=1.25/(1+exp((V+32.8)/0.5))+0.02;     // 0.65 - LA, 1.25 - RA
   s1+=HT*(s1m-s1)/ts1;

   s2m=1./(1+exp((V+32)/8.4)); 
   ts2=6.65/(1+exp((V+62.8)/0.5))+0.045;   
   s2+=HT*(s2m-s2)/ts2;
    
   s3m=((1./(1+exp((V+50.67)/27.38)))+0.666)/1.666;
   ts3=(7.5/(1+exp((V+23.0)/0.5)))+0.5;      
   s3+=HT*(s3m-s3)/ts3;

   Ito=25.01*r1*(0.8*pow(s1,3)+0.2*pow(s2,3))*(V-Ek)*(0.6*pow(s3,6)+0.4);    
   Isus = 1.25*(V+85); //+0.25*(V-5);  
   Ito+=Isus;

   Ip=64.41*Kc/(Kc+1)*(pow(Nai,1.5)/(pow(Nai,1.5)+pow(11,1.5)))*(1.6/(1.5+exp(-(V+60)/40.)));

   
   INaCa=0.02*(pow(Nai,3)*2.5*exp(0.450*V/RTONF)-pow(140,3)*Cai*exp(V*(0.45-1)/RTONF))/(1+0.0003*(Cai*pow(140,3)+2.5*pow(Nai,3)));  


   if(fabs(V+44.4) < 0.0001)
      Am=460.*12.673;
   else    
     Am=-460.*(V+44.4)/(exp(-(V+44.4)/12.673)-1);
   Bm=18400.0*exp(-(V+44.4)/12.673);
   m = Am/(Am+Bm);

   Ah=44.9*exp(-(V+66.9)/5.57);
   Bh=1491.0/(1+323.3*exp(-(V+94.6)/12.9));

   th1=0.03/(1+exp((V+40)/6.0))+0.00015;  
   th2=0.12/(1+exp((V+60)/2.0))+0.00045;  
   hm=Ah/(Ah+Bh);

   h1=h1+HT*(hm-h1)/th1;
   h2=h2+HT*(hm-h2)/th2;

   if(fabs(V) > 0.0001) 
      INa=0.75*0.0014*pow(m,3)*(0.635*h1+0.365*h2)*140*V*(Faraday/RTONF)*(exp((V-ENa)/RTONF)-1)/(exp(V/RTONF)-1); 
   else INa=0.75*0.0014*pow(m,3)*(0.635*h1+0.365*h2)*140*Faraday*(exp((V-ENa)/RTONF)-1);

   
   IbNa=0.02*(V-ENa)+0.02*(V-ECa);  


   V=V+10;

   Adl=-16.72*(V+35)/(exp(-(V+35)/2.5)-1)-50.0*V/(exp(-V/4.808)-1);
   if(fabs(V+35) < 0.0001)
     Adl=16.72*2.5-50.0*V/(exp(-V/4.808)-1);
   if(fabs(V) < 0.0001) 
     Adl=-16.72*(V+35)/(exp(-(V+35)/2.5)-1)+50.0*4.808;
       
   if(fabs(V-5) < 0.0001)
     Bdl=4.48*2.5;
   else
     Bdl=4.48*(V-5)/(exp((V-5)/2.5)-1.);
   tdl=1/(Adl+Bdl);

   dlm=1./(1+exp(-(V+6)/6.6));
   dL+=HT*(dlm-dL)/tdl;

   V=V-10;   

   if(fabs(V+28) < 0.0001) 
     Afl=8.49*4;
   else
     Afl=8.49*(V+28)/(exp((V+28)/6.5)-1.);
   Bfl=67.922/(1+exp(-(V+28)/6.5));

   flm=Afl/(Afl+Bfl);
   tfl=1.0/(Afl+Bfl);
   fL1+=HT*(flm-fL1)/tfl;

   ICaL=7.2*(dL*fL1+1.0/(1+exp(-(V-18.0)/12.0)))*(V-45);


   return ((IKf+IKs+Ik1+Ito+Ip+INaCa+INa+IbNa+ICaL)/Cm);
}
