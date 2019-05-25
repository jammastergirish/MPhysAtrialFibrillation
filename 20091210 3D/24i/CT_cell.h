#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Atria_Variables.h"

float I_tot_CT (int x, int y, int z)
{
   float IKf,IKs,IK,Ik1,Ito,IbNa,Ip,INaCa,IbCa,INa,ICap,ICaL,ICaT,Isus;

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
   float hm;
   float th1, th2;
   float Adl, Bdl;
   float dlm, tdl;
   float Afl, Bfl;
   float flm, tfl;
   float Adt, Bdt;
   float dtm, tdt;
   float Aft, Bft;
   float ftm, tft;
   float E0=V[x][y][z]-Ek+3.6;

   IKf=3.5*Pa[x][y][z]*Pi[x][y][z]*(V[x][y][z]-Ek);     

   Apa=9.0*exp(V[x][y][z]/25.371);
   Bpa=1.3*exp(-V[x][y][z]/13.026);
   pam=1./(1+exp(-(V[x][y][z]+5.1)/7.4));
   tpa=1./(Apa+Bpa);
   Pa[x][y][z]+=HT*(pam-Pa[x][y][z])/tpa;
      
   Api=100.*exp(-V[x][y][z]/54.645);
   Bpi=656.*exp(V[x][y][z]/106.157);
   pim=1./(1+exp((V[x][y][z]+47.3921)/18.6603));
   tpi=1.0/(Api+Bpi);
   Pi[x][y][z]+=HT*(pim-Pi[x][y][z])/tpi;

   IKs=2.5*n[x][y][z]*(V[x][y][z]-Ek);      

   An=1.66*exp(V[x][y][z]/69.452);
   Bn=0.3*exp(-V[x][y][z]/21.826);
   nm=1.0/(1+exp(-(V[x][y][z]-0.9)/13.8));
   tn=1./(An+Bn)+0.06;   
   n[x][y][z]+=HT*(nm-n[x][y][z])/tn;

   IK=IKf+IKs;

   Ik1=2.0*5.08*pow(Kc[x][y][z]/(Kc[x][y][z]+0.59),3)*(V[x][y][z]-Ek)/(1+exp(1.393*E0/RTONF));     // CT - 2.0, PM - 2.5

   Ito=0.2*50.02*r1[x][y][z]*(0.590*pow(s1[x][y][z],3)+0.410*pow(s2[x][y][z],3))*(0.600*pow(s3[x][y][z],6)+0.4)*(V[x][y][z]-Ek);     // CT - 0.2, PM - 0.35

   Isus = 1.4*(V[x][y][z]+70);  // CT - 1.4, PM - 2.4
   Ito+=Isus;

   Ar=386.6*exp(V[x][y][z]/12.0);
   Br=8.011*exp(-V[x][y][z]/7.2);
   rm=1/(1+exp(-(V[x][y][z]+15.0)/5.633));   
   tr=1./(Ar+Br)+0.0004;
   r1[x][y][z]+=HT*(rm-r1[x][y][z])/tr;

   s1m=1./(1+exp((V[x][y][z]+28.29)/7.06));
   ts1=0.5466/(1+exp((V[x][y][z]+32.8)/0.1))+0.0204;
   s1[x][y][z]+=HT*(s1m-s1[x][y][z])/ts1;

   s2m= 1./(1+exp((V[x][y][z]+28.29)/7.06)); 
   ts2=5.75/(1+exp((V[x][y][z]+32.8)/0.1))+0.45;
   s2[x][y][z]+=HT*(s2m-s2[x][y][z])/ts2;
    
   s3m=((1./(1+exp((V[x][y][z]+50.67)/27.38)))+0.666)/1.666;
   ts3=(7.5/(1+exp((V[x][y][z]+23.0)/0.5)))+0.5;      
   s3[x][y][z]+=HT*(s3m-s3[x][y][z])/ts3;

   Ip=64.41*Kc[x][y][z]/(Kc[x][y][z]+1)*(pow(Nai[x][y][z],1.5)/(pow(Nai[x][y][z],1.5)+pow(11,1.5)))*(1.6/(1.5+exp(-(V[x][y][z]+60)/40.)));

   INaCa=0.02*(pow(Nai[x][y][z],3)*2.5*exp(0.450*V[x][y][z]/RTONF)-pow(140,3)*Cai[x][y][z]*exp(V[x][y][z]*(0.45-1)/RTONF))/(1+0.0003*(Cai[x][y][z]*pow(140,3)+2.5*pow(Nai[x][y][z],3)));  

   if(fabs(V[x][y][z]+44.4) < 0.0001)
      Am=460.*12.673;
   else    
     Am=-460*(V[x][y][z]+44.4)/(exp(-(V[x][y][z]+44.4)/12.673)-1);
   Bm=18400.0*exp(-(V[x][y][z]+44.4)/12.673);
   m[x][y][z] = Am/(Am+Bm);
//   m[x][y][z]+=HT*(Am*(1-m[x][y][z])-Bm*m[x][y][z]);

   Ah=44.9*exp(-(V[x][y][z]+66.9)/5.57);
   Bh=1491.0/(1+323.3*exp(-(V[x][y][z]+94.6)/12.9));

   th1=0.03/(1+exp((V[x][y][z]+40)/6.0))+0.00015;  //0.00035 - Lindblad   
   th2=0.12/(1+exp((V[x][y][z]+60)/2.0))+0.00045;  //0.00295 - Lindblab    

   hm=Ah/(Ah+Bh);
   h1[x][y][z]+=HT*(hm-h1[x][y][z])/th1; 
   h2[x][y][z]+=HT*(hm-h2[x][y][z])/th2;

   if(fabs(V[x][y][z]) > 0.0001) 
      INa=0.75*0.0014*pow(m[x][y][z],3)*(0.635*h1[x][y][z]+0.365*h2[x][y][z])*140*V[x][y][z]*(Faraday/RTONF)*(exp((V[x][y][z]-ENa)/RTONF)-1)/(exp(V[x][y][z]/RTONF)-1); // multiply by 0.75 for instant activation
   else INa=0.75*0.0014*pow(m[x][y][z],3)*(0.635*h1[x][y][z]+0.365*h2[x][y][z])*140*Faraday*(exp((V[x][y][z]-ENa)/RTONF)-1);   

   ICap=9.509*(Cai[x][y][z]/(Cai[x][y][z]+0.0002));   

   V[x][y][z]=V[x][y][z]+10;    

   Adl=-16.72*(V[x][y][z]+35)/(exp(-(V[x][y][z]+35)/2.5)-1)-50.0*V[x][y][z]/(exp(-V[x][y][z]/4.808)-1);
   if(fabs(V[x][y][z]+35) < 0.0001)
     Adl=16.72*2.5-50.0*V[x][y][z]/(exp(-V[x][y][z]/4.808)-1);
   if(fabs(V[x][y][z]) < 0.0001) 
     Adl=-16.72*(V[x][y][z]+35)/(exp(-(V[x][y][z]+35)/2.5)-1)+50.0*4.808;
       
   if(fabs(V[x][y][z]-5.) < 0.0001)
     Bdl=4.48*2.5;
   else
     Bdl=4.48*(V[x][y][z]-5)/(exp((V[x][y][z]-5)/2.5)-1.);
   tdl=1/(Adl+Bdl);

   dlm=1./(1+exp(-(V[x][y][z]+0.95)/6.6));
   dL[x][y][z]+=HT*(dlm-dL[x][y][z])/tdl;

   V[x][y][z]=V[x][y][z]-10;   

   V[x][y][z]=V[x][y][z]-10;   

   if(fabs(V[x][y][z]+28) < 0.0001) 
     Afl=8.49*4;
   else
     Afl=8.49*(V[x][y][z]+28)/(exp((V[x][y][z]+28)/4)-1.);
   Bfl=67.922/(1+exp(-(V[x][y][z]+28)/4));

   flm=Afl/(Afl+Bfl);
   tfl=1.0/(Afl+Bfl);
   fL1[x][y][z]+=HT*(flm-fL1[x][y][z])/tfl;

   V[x][y][z]=V[x][y][z]+10;   

   V[x][y][z]=V[x][y][z]+10;

   ICaL=1.8*4.0*(dL[x][y][z]*fL1[x][y][z]+1.0/(1+exp(-(V[x][y][z]-33.0)/12.0)))*(V[x][y][z]-60);      // CT - 1.8, PM - 2.1

   V[x][y][z]=V[x][y][z]-10;   

   Adt=674.173*exp((V[x][y][z]+23.3)/30.);
   Bdt=674.173*exp(-(V[x][y][z]+23.3)/30.);
   tdt=1/(Adt+Bdt);
   dtm=1./(1+exp(-(V[x][y][z]+23.0)/6.1));
   dT[x][y][z]+=HT*(dtm-dT[x][y][z])/tdt;

   Aft=9.637*exp(-(V[x][y][z]+75)/83.3);
   Bft=9.637*exp((V[x][y][z]+75)/15.38);
   tft=1.0/(Aft+Bft);
   ftm=Aft/(Aft+Bft);
   fT[x][y][z]+=HT*(ftm-fT[x][y][z])/tft;

   ICaT=6.0*dT[x][y][z]*fT[x][y][z]*(V[x][y][z]-38.0);

   IbNa=0.02*(V[x][y][z]-ENa);  // 0.02 - CT, 0.03 - PM

   IbCa=0.02*(V[x][y][z]-ECa);  // 0.02 - CT, 0.03 - PM

   return ((IKf+IKs+Ik1+Ito+Ip+INaCa+INa+IbNa+ICap+ICaL+ICaT+IbCa)/Cm);

}
