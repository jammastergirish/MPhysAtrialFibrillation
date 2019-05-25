#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Atria_Variables.h"

float I_tot_BB (int x, int y, int z)
{
   float IKf,IKs,IK,Ik1,Ito,IbNa,Ip,INaCa,IbCa,INa, ICap,ICaL,ICaT,Isus /*, Iup,Itr,Irel */;

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
   float flm, tfl1, tfl2;
//   float dfca;
   float Adt, Bdt;
   float dtm, tdt;
   float Aft, Bft;
   float ftm, tft;

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

   Ik1=(5.08*pow(5.0/(5.0+0.59),3)/(1+exp(1.393*(V[x][y][z]-Ek-12)/RTONF))+0.05)*(V[x][y][z]-Ek-12);

   Ito=0.54*50.02*r1[x][y][z]*(0.590*pow(s1[x][y][z],3)+0.410*pow(s2[x][y][z],3))*(0.600*pow(s3[x][y][z],6)+0.4)*(V[x][y][z]-Ek);     
   Isus = 0.3*(V[x][y][z]+85);
   
   Ito+=Isus;
   
   Ar=386.6*exp(V[x][y][z]/12.0);
   Br=8.011*exp(-V[x][y][z]/7.2);
   rm=1/(1+exp(-(V[x][y][z]+15.0)/5.633));   
   tr=1./(Ar+Br)+0.0004;
   r1[x][y][z]+=HT*(rm-r1[x][y][z])/tr;

   s1m=1./(1+exp((V[x][y][z]+28.2)/5.3));
   ts1=0.1199/(1+exp((V[x][y][z]+32.8)/0.1))+0.0157; // t1_r-t1_d / (Exp) + t1_d
   s1[x][y][z]+=HT*(s1m-s1[x][y][z])/ts1;

   s2m= 1./(1+exp((V[x][y][z]+28.2)/5.3)); 
   ts2=3.1616/(1+exp((V[x][y][z]+32.8)/0.1))+0.1103; // t2_r+t2_d / (Exp) + t2_d
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


   Ah=44.9*exp(-(V[x][y][z]+66.9)/5.57);
   Bh=1491.0/(1+323.3*exp(-(V[x][y][z]+94.6)/12.9));

   th1=0.03/(1+exp((V[x][y][z]+40)/6.0))+0.00015;  //0.00035 - Lindblad   
   th2=0.12/(1+exp((V[x][y][z]+60)/2.0))+0.00045;  //0.00295 - Lindblab    

   hm=Ah/(Ah+Bh);
   h1[x][y][z]+=HT*(hm-h1[x][y][z])/th1; 
   h2[x][y][z]+=HT*(hm-h2[x][y][z])/th2;

   if(fabs(V[x][y][z]) > 0.0001) 
        INa=0.0014*pow(m[x][y][z],3)*(0.635*h1[x][y][z]+0.365*h2[x][y][z])*140*V[x][y][z]*(Faraday/RTONF)*(exp((V[x][y][z]-ENa)/RTONF)-1)/(exp(V[x][y][z]/RTONF)-1); // multiply by 0.75 for instant activation
   else INa=0.0014*pow(m[x][y][z],3)*(0.635*h1[x][y][z]+0.365*h2[x][y][z])*140*Faraday*(exp((V[x][y][z]-ENa)/RTONF)-1);   

   ICap=9.509*(Cai[x][y][z]/(Cai[x][y][z]+0.0002));   
   
   Adl=-16.72*(V[x][y][z]+45)/(exp(-(V[x][y][z]+45)/2.5)-1)-50.0*(V[x][y][z]+10)/(exp(-(V[x][y][z]+10)/4.808)-1);
   if(fabs(V[x][y][z]+45) < 0.0001)
     Adl=16.72*2.5-50.0*(V[x][y][z]+10)/(exp(-(V[x][y][z]+10)/4.808)-1);
   if(fabs(V[x][y][z]+10) < 0.0001) 
     Adl=-16.72*(V[x][y][z]+45)/(exp(-(V[x][y][z]+45)/2.5)-1)+50.0*4.808;
       
   if(fabs(V[x][y][z]+5.) < 0.0001)
     Bdl=4.48*2.5;
   else
     Bdl=4.48*(V[x][y][z]+5)/(exp((V[x][y][z]+5)/2.5)-1.);
   tdl=1/(Adl+Bdl);

   dlm=1./(1+exp(-(V[x][y][z]-7.9)/6.3));
   dL[x][y][z]+=HT*(dlm-dL[x][y][z])/tdl;

  flm=1./(1+exp((V[x][y][z]+20.4)/6.3));
  tfl1=(0.0291/(1+exp((V[x][y][z]+18)/4.0)))+0.0132;
  fL1[x][y][z]+=HT*(flm-fL1[x][y][z])/tfl1;
/*
  dfca=1.0/(1+Cai[x][y][z]/0.00015);
  fca[x][y][z]+=HT*(dfca-fca[x][y][z])/0.02;
*/
  ICaL=0.6*16*dL[x][y][z]*fL1[x][y][z]*(V[x][y][z]-61.4); 


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

   ICaT=6*dT[x][y][z]*fT[x][y][z]*(V[x][y][z]-38.0);

   IbNa=0.02*(V[x][y][z]-ENa);  // 0.02 - CT, 0.03 - PM

   IbCa=0.02*(V[x][y][z]-ECa);  // 0.02 - CT, 0.03 - PM

/*
   Nai[x][y][z]+=HT*(-3*Ip-3*INaCa-IbNa-INa)/(Faraday*Vi);

   Iup=2800.0*(((Cai[x][y][z]/0.0003)-(0.4*0.4*Caup[x][y][z]/0.5))/(((Cai[x][y][z]+0.0003)/0.0003)+(0.4*(Caup[x][y][z]+0.5)/0.5)));

   Irel=200000.0*pow(F2[x][y][z]/(F2[x][y][z]+0.25),2)*(Carel[x][y][z]-Cai[x][y][z]);    

   Itr=(Caup[x][y][z]-Carel[x][y][z])*2.*Faraday*Vup/0.01;

   dfac=200000.0*Cai[x][y][z]*(1-fac)-476*fac;
   dfaTc=78400*Cai[x][y][z]*(1-faTc)-392*faTc;
   dfaTmgc=200000*Cai[x][y][z]*(1-faTmgc-faTmgm)-6.6*faTmgc;
   dfaTmgm=2000*2.5*(1-faTmgc-faTmgm)-666*faTmgm;
   dfaCalse=480*Carel[x][y][z]*(1-faCalse)-400*faCalse;       
   dfab=0.08*dfaTc+0.16*dfaTmgc+0.045*dfac;

   Caup[x][y][z]+=HT*(Iup-Itr)/(2.*Faraday*Vup);
   Carel[x][y][z]+=HT*((Itr-Irel)/(2.*Faraday*Vrel)-31.0*(480*Carel[x][y][z]*(1-faCalse[x][y][z])-400*faCalse[x][y][z]));   
   Cai[x][y][z]+=HT*((2.*INaCa-ICaL-ICaT-ICap-IbCa-Iup+Irel)/(2.*VCa*Faraday)-dfab);

   fac+=HT*dfac;
   faTc+=HT*dfaTc;
   faTmgc+=HT*dfaTmgc;
   faTmgm+=HT*dfaTmgm;
   faCalse+=HT*dfaCalse;

   Kc[x][y][z]+=HT*(-2.*Ip+IKf+IKs+Ito+Ik1)/(Faraday*Vc);
   Ki[x][y][z]+=HT*(2.*Ip-IKf-IKs-Ito-Ik1)/(Faraday*Vi);

   ract[x][y][z]=240.0*exp((V[x][y][z]-20.)/12.5)+203.8*pow(Cai[x][y][z]/(Cai[x][y][z]+0.0003),4); 
   rinact=33.96+339.6*pow(Cai[x][y][z]/(Cai[x][y][z]+0.0003),4);

   dF1[x][y][z]=0.815*F3[x][y][z]-ract[x][y][z]*F1[x][y][z];
   dF2[x][y][z]=ract[x][y][z]*F1[x][y][z]-rinact*F2[x][y][z];
   dF3[x][y][z]=rinact*F2[x][y][z]-0.815*F3[x][y][z];

   F1[x][y][z]+=HT*dF1[x][y][z];
   F2[x][y][z]+=HT*dF2[x][y][z];
   F3[x][y][z]+=HT*dF3[x][y][z];
*/
   return ((IKf+IKs+Ik1+Ito+Ip+INaCa+INa+IbNa+ICap+ICaL+ICaT+IbCa)/Cm);

}
