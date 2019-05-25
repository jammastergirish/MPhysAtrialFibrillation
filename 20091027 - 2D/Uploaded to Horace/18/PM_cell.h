#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Atria_Variables.h"


float I_tot_PM (int x, int y)
{
   float IKf,IKs,Ik1,Ito,IbNa,Ip,INaCa,IbCa,INa,ICaL,Isus;

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
   float hm, mm;
   float th1, th2;
   float Adl, Bdl;
   float dlm, tdl;
   float Afl, Bfl;
   float flm, tfl;


   Apa=9.0*exp(V[x][y]/25.371);
   Bpa=1.3*exp(-V[x][y]/13.026);
   pam=1./(1+exp(-(V[x][y]+5.1)/7.4));
   tpa=1./(Apa+Bpa);
   Pa[x][y]=Pa[x][y]+HT*(pam-Pa[x][y])/tpa;
      
   Api=100.*exp(-V[x][y]/54.645);
   Bpi=656.*exp(V[x][y]/106.157);
   pim=1./(1+exp((V[x][y]+47.3921)/18.6603));
   tpi=1.0/(Api+Bpi);
   Pi[x][y]=Pi[x][y]+HT*(pim-Pi[x][y])/tpi;

   IKf=3.5*Pa[x][y]*Pi[x][y]*(V[x][y]-Ek);     

   
   An=1.66*exp(V[x][y]/69.452);
   Bn=0.3*exp(-V[x][y]/21.826);
   nm=1.0/(1+exp(-(V[x][y]-0.9)/13.8));
   tn=1./(An+Bn)+0.06;
   n[x][y]=n[x][y]+HT*(nm-n[x][y])/tn;

   IKs=2.5*n[x][y]*(V[x][y]-Ek);     


   Ik1=12.7*pow(Kc[x][y]/(Kc[x][y]+0.59),3)*(V[x][y]-Ek)/(1+exp(1.393*(V[x][y]-Ek+3.6)/RTONF));     


   Ar=386.6*exp(V[x][y]/12.0);
   Br=8.011*exp(-V[x][y]/7.2);
   rm=1/(1+exp(-(V[x][y]+15.0)/5.633));
   tr=1./(Ar+Br)+0.0004;
   r1[x][y]=r1[x][y]+HT*(rm-r1[x][y])/tr;

   s1m=1./(1+exp((V[x][y]+28.29)/7.06));
   ts1=0.5466/(1+exp((V[x][y]+32.8)/0.1))+0.0204;
   s1[x][y]=s1[x][y]+HT*(s1m-s1[x][y])/ts1;

   s2m= 1./(1+exp((V[x][y]+28.29)/7.06)); 
   ts2=5.75/(1+exp((V[x][y]+32.8)/0.1))+0.45/(1+exp(-(V[x][y]-13.54)/13.97));
   s2[x][y]=s2[x][y]+HT*(s2m-s2[x][y])/ts2;
    
   s3m=((1./(1+exp((V[x][y]+50.67)/27.38)))+0.666)/1.666;
   ts3=(7.5/(1+exp((V[x][y]+23.0)/0.5)))+0.5;
   s3[x][y]=s3[x][y]+HT*(s3m-s3[x][y])/ts3;

   Ito=17.5*r1[x][y]*(0.590*pow(s1[x][y],3)+0.410*pow(s2[x][y],3))*(0.600*pow(s3[x][y],6)+0.4)*(V[x][y]-Ek);
   Isus=2.4*(V[x][y]+70);


   Ip=64.41*Kc[x][y]/(Kc[x][y]+1)*(pow(Nai[x][y],1.5)/(pow(Nai[x][y],1.5)+pow(11,1.5)))*(1.6/(1.5+exp(-(V[x][y]+60)/40.)));

   INaCa=0.02*(pow(Nai[x][y],3)*2.5*exp(0.450*V[x][y]/RTONF)-pow(140,3)*Cai[x][y]*exp(V[x][y]*(0.45-1)/RTONF))/(1+0.0003*(Cai[x][y]*pow(140,3)+2.5*pow(Nai[x][y],3)));

   IbNa=0.03*(V[x][y]-ENa);  // 0.09

   IbCa=0.03*(V[x][y]-ECa);  // 0.09

   
   if(fabs(V[x][y]+44.4) < 0.0001)
      Am=460.*12.673;
   else    
     Am=-460.*(V[x][y]+44.4)/(exp(-(V[x][y]+44.4)/12.673)-1);
   Bm=18400.0*exp(-(V[x][y]+44.4)/12.673);
   m[x][y] = Am/(Am+Bm);

   Ah=44.9*exp(-(V[x][y]+66.9)/5.57);
   Bh=1491.0/(1+323.3*exp(-(V[x][y]+94.6)/12.9));

   th1=0.03/(1+exp((V[x][y]+40)/6.0))+0.00015;
   th2=0.12/(1+exp((V[x][y]+60)/2.0))+0.00045;
   hm=Ah/(Ah+Bh);
   h1[x][y]=h1[x][y]+HT*(hm-h1[x][y])/th1;
   h2[x][y]=h2[x][y]+HT*(hm-h2[x][y])/th2;

   if(fabs(V[x][y]) > 0.0001) 
      INa=0.75*0.0014*pow(m[x][y],3)*(0.635*h1[x][y]+0.365*h2[x][y])*140*V[x][y]*(Faraday/RTONF)*(exp((V[x][y]-ENa)/RTONF)-1)/(exp(V[x][y]/RTONF)-1);
   else INa=0.75*0.0014*pow(m[x][y],3)*(0.635*h1[x][y]+0.365*h2[x][y])*140*Faraday*(exp((V[x][y]-ENa)/RTONF)-1);   // 0.8


   V[x][y]=V[x][y]+10;

   Adl=-16.72*(V[x][y]+35)/(exp(-(V[x][y]+35)/2.5)-1)-50.0*V[x][y]/(exp(-V[x][y]/4.808)-1);
   if(fabs(V[x][y]+35) < 0.0001)
     Adl=16.72*2.5-50.0*V[x][y]/(exp(-V[x][y]/4.808)-1);
   if(fabs(V[x][y]) < 0.0001) 
     Adl=-16.72*(V[x][y]+35)/(exp(-(V[x][y]+35)/2.5)-1)+50.0*4.808;
       
   if(fabs(V[x][y]-5.) < 0.0001)
     Bdl=4.48*2.5;
   else
     Bdl=4.48*(V[x][y]-5)/(exp((V[x][y]-5)/2.5)-1.);
   tdl=1/(Adl+Bdl);

   dlm=1./(1+exp(-(V[x][y]+0.95)/6.6));
   dL[x][y]=dL[x][y]+HT*(dlm-dL[x][y])/tdl;

   V[x][y]=V[x][y]-10;
   V[x][y]=V[x][y]-10;

   if(fabs(V[x][y]+28) < 0.0001) 
     Afl=8.49*4;
   else
     Afl=8.49*(V[x][y]+28)/(exp((V[x][y]+28)/4)-1.);
   Bfl=67.922/(1+exp(-(V[x][y]+28)/4));

   flm=Afl/(Afl+Bfl);
   tfl=1.0/(Afl+Bfl);
   fL1[x][y]=fL1[x][y]+HT*(flm-fL1[x][y])/tfl;

   V[x][y]=V[x][y]+10;

   ICaL=8.4*(dL[x][y]*fL1[x][y]+1.0/(1+exp(-(V[x][y]-23.0)/12.0)))*(V[x][y]-50);


   return ((IKf+IKs+Ik1+Ito+Isus+Ip+INaCa+INa+IbNa+IbCa+ICaL)/Cm);
}
