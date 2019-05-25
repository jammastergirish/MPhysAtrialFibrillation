#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Atria_Variables.h"

float I_tot_PV (int x, int y)
{
   float IKf,IKs,IK,Ik1,Ito,IbNa,Ip,INaCa,IbCa,INa, ICap,ICaL,ICaT,Isus,IfNa,IfK,IKH_K, IKH_Na, IKH,Iup,Itr,Irel;
   float ract, rinact,dfac,dfaTc,dfaTmgc,dfaTmgm,dfaCalse,dfab;

   float Ay, By;
   float ym, ty;
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
   float xsusm;
   float Am, Bm;
   float Ah, Bh;
   float hm;
   float th1, th2;
   float Adl, Bdl;
   float dlm, tdl;
   float Afl, Bfl;
   float flm, tfl1, tfl2;
   float dfca;
   float Adt, Bdt;
   float dtm, tdt;
   float Aft, Bft;
   float ftm, tft;
   float E0=V[x][y]-Ek+3.6;
   
   
   Ay = exp(-(2.9 + (0.04*V[x][y])));
   By  = exp((3.6 + (0.11*V[x][y])));

   ym = 1/(1 + exp((V[x][y]+100)/20));
   ty = 0.4/(Ay+By);

   y2[x][y] += HT*(ym-y2[x][y])/ty;

   IKH_K= 0.1*ym*(V[x][y]-Ek);
   IKH_Na = 0.008*ym*(V[x][y]-ENa);
   IKH = (IKH_K + IKH_Na);
   
   IKf=4*0.0035*Pa[x][y]*Pi[x][y]*(V[x][y]-Ek);    

   Apa=9.0*exp(V[x][y]/25.371);
   Bpa=1.3*exp(-V[x][y]/13.026);
   pam=1./(1+exp(-(V[x][y]+5.1)/7.4));
   tpa=1./(Apa+Bpa);
   Pa[x][y]+=HT*(pam-Pa[x][y])/tpa;
      
   Api=100.*exp(-V[x][y]/54.645);
   Bpi=656.*exp(V[x][y]/106.157);
   pim=1./(1+exp((V[x][y]+47.3921)/18.6603));
   tpi=1.0/(Api+Bpi);
   Pi[x][y]+=HT*(pim-Pi[x][y])/tpi;

   IKs=2.5*n[x][y]*(V[x][y]-Ek);      

   An=1.66*exp(V[x][y]/69.452);
   Bn=0.3*exp(-V[x][y]/21.826);
   nm=1.0/(1+exp(-(V[x][y]-0.9)/13.8));
   tn=1./(An+Bn)+0.06;   
   n[x][y]+=HT*(nm-n[x][y])/tn;

   IK=IKf+IKs;

   Ik1=(0.9*5.08*pow(5.0/(5.0+0.59),3)/(1+exp(1.393*(V[x][y]-Ek-12)/RTONF))+0.05)*(V[x][y]-Ek-12);

   Ito=0.45*0.515*50.02*r1[x][y]*(0.8*pow(s1[x][y],3)+0.2*pow(s2[x][y],3))*/*(0.600*pow(s3[x][y],6)+0.4)*/(V[x][y]-Ek);     // CT - 0.2, PM - 0.35
   Isus = 0.46*(V[x][y]+85)+0.0*(V[x][y]-5);
   
   Ito+=Isus;
   
   Ar=386.6*exp(V[x][y]/12.0);
   Br=8.011*exp(-V[x][y]/7.2);
   rm=1/(1+exp(-(V[x][y]+15.0)/5.633));   
   tr=1./(Ar+Br)+0.0004;
   r1[x][y]+=HT*(rm-r1[x][y])/tr;

   s1m=1./(1+exp((V[x][y]+29.4)/5.1));
   ts1=0.1155/(1+exp((V[x][y]+32.8)/0.1))+0.0147; // t1_r-t1_d / (Exp) + t1_d
   s1[x][y]+=HT*(s1m-s1[x][y])/ts1;

   s2m= 1./(1+exp((V[x][y]+29.4)/5.1)); 
   ts2=3.0967/(1+exp((V[x][y]+32.8)/0.1))+0.0926; // t2_r+t2_d / (Exp) + t2_d
   s2[x][y]+=HT*(s2m-s2[x][y])/ts2;
    
   s3m=((1./(1+exp((V[x][y]+50.67)/27.38)))+0.666)/1.666;
   ts3=(7.5/(1+exp((V[x][y]+23.0)/0.5)))+0.5;      
   s3[x][y]+=HT*(s3m-s3[x][y])/ts3;

   Ip=64.41*Kc[x][y]/(Kc[x][y]+1)*(pow(Nai[x][y],1.5)/(pow(Nai[x][y],1.5)+pow(11,1.5)))*(1.6/(1.5+exp(-(V[x][y]+60)/40.)));

//Regular AP
//INaCa=0.02*(pow(Nai[x][y],3)*2.5*exp(0.450*V[x][y]/RTONF)-pow(140,3)*Cai[x][y]*exp(V[x][y]*(0.45-1)/RTONF))/(1+0.0003*(Cai[x][y]*pow(140,3)+2.5*pow(Nai[x][y],3)));  
   
//SAP
INaCa = 1.5*(0.02*(pow(Nai[x][y],3)*7.5*exp(0.35*V[x][y]/RTONF)-pow(140,3)*Cai[x][y]*exp(V[x][y]*(0.45-1)/RTONF))/(1+0.0003*(Cai[x][y]*pow(140,3)+2.5*pow(Nai[x][y],3)))); //1.5

   if(fabs(V[x][y]+44.4) < 0.0001)
      Am=460.*12.673;
   else    
     Am=-460*(V[x][y]+44.4)/(exp(-(V[x][y]+44.4)/12.673)-1);
   Bm=18400.0*exp(-(V[x][y]+44.4)/12.673);
   m[x][y] = Am/(Am+Bm);
//   m[x][y]+=HT*(Am*(1-m[x][y])-Bm*m[x][y]);

   Ah=44.9*exp(-(V[x][y]+66.9)/5.57);
   Bh=1491.0/(1+323.3*exp(-(V[x][y]+94.6)/12.9));

   th1=0.03/(1+exp((V[x][y]+40)/6.0))+0.00015;  //0.00035 - Lindblad   
   th2=0.12/(1+exp((V[x][y]+60)/2.0))+0.00045;  //0.00295 - Lindblab    

   hm=Ah/(Ah+Bh);
   h1[x][y]+=HT*(hm-h1[x][y])/th1; 
   h2[x][y]+=HT*(hm-h2[x][y])/th2;

   if(fabs(V[x][y]) > 0.0001) 
      INa=0.75*0.0014*pow(m[x][y],3)*(0.635*h1[x][y]+0.365*h2[x][y])*140*V[x][y]*(Faraday/RTONF)*(exp((V[x][y]-ENa)/RTONF)-1)/(exp(V[x][y]/RTONF)-1); // multiply by 0.75 for instant activation
   else INa=0.75*0.0014*pow(m[x][y],3)*(0.635*h1[x][y]+0.365*h2[x][y])*140*Faraday*(exp((V[x][y]-ENa)/RTONF)-1);   

   ICap=9.509*(Cai[x][y]/(Cai[x][y]+0.0002));   
   
   Adl=-16.72*(V[x][y]+45)/(exp(-(V[x][y]+45)/2.5)-1)-50.0*(V[x][y]+10)/(exp(-(V[x][y]+10)/4.808)-1);
   if(fabs(V[x][y]+45) < 0.0001)
     Adl=16.72*2.5-50.0*(V[x][y]+10)/(exp(-(V[x][y]+10)/4.808)-1);
   if(fabs(V[x][y]+10) < 0.0001) 
     Adl=-16.72*(V[x][y]+45)/(exp(-(V[x][y]+45)/2.5)-1)+50.0*4.808;
       
   if(fabs(V[x][y]+5.) < 0.0001)
     Bdl=4.48*2.5;
   else
     Bdl=4.48*(V[x][y]+5)/(exp((V[x][y]+5)/2.5)-1.);
   tdl=1/(Adl+Bdl);

   dlm=1./(1+exp(-(V[x][y]+7.5)/5.5));
   dL[x][y]+=HT*(dlm-dL[x][y])/tdl;


  flm=1./(1+exp((V[x][y]+23.4)/5.2));
  tfl1=(0.0336/(1+exp((V[x][y]+18)/4.0)))+0.0142;
  tfl2=(3.6908/(1+exp((V[x][y]+18)/4)))+0.0534;
  fL1[x][y]+=HT*(flm-fL1[x][y])/tfl1;
  fL2[x][y]+=HT*(flm-fL2[x][y])/tfl2;

  dfca=1.0/(1+Cai[x][y]/0.00001);  //0.00025
  fca[x][y]+=HT*(dfca-fca[x][y])/0.04;

//   ICaL=18*(dL[x][y]*((0.80*fL1[x][y])+(0.2*fL2[x][y])))*fca*(V[x][y]-61.4);   
     ICaL=25*dL[x][y]*fL1[x][y]*fca[x][y]*(V[x][y]-58.3);
 
   Adt=674.173*exp((V[x][y]+23.3)/30.0);
   Bdt=674.173*exp(-(V[x][y]+23.3)/30.0); //23.3 
   tdt=1/(Adt+Bdt);
   dtm=1./(1+exp(-(V[x][y]+23.0)/5.1)); //23 6.1
   dT[x][y]+=HT*(dtm-dT[x][y])/tdt;

   Aft=9.637*exp(-(V[x][y]+75)/83.3);// 75
   Bft=9.637*exp((V[x][y]+75)/16.3);// 15.38
   tft=1.0/(Aft+Bft);
   ftm=Aft/(Aft+Bft);
   fT[x][y]+=HT*(ftm-fT[x][y])/tft;

   ICaT=1.65*8.0*dT[x][y]*fT[x][y]*(V[x][y]-38.0);

   IbNa=0.02*(V[x][y]-ENa);  // 0.02 - CT, 0.03 - PM

   IbCa=0.02*(V[x][y]-ECa);  // 0.02 - CT, 0.03 - PM


   Nai[x][y]+=HT*(-3*Ip-3*INaCa-IbNa-INa)/(Faraday*Vi);

   Iup=3*2800.0*(((Cai[x][y]/0.0003)-(0.4*0.4*Caup[x][y]/0.5))/(((Cai[x][y]+0.0003)/0.0003)+(0.4*(Caup[x][y]+0.5)/0.5)));

   Irel=4*200000.0*pow(F2[x][y]/(F2[x][y]+0.25),2)*(Carel[x][y]-Cai[x][y]);    

   Itr=(Caup[x][y]-Carel[x][y])*2.*Faraday*Vup/0.01;

   dfac=200000.0*Cai[x][y]*(1-fac[x][y])-476*fac[x][y];
   dfaTc=78400*Cai[x][y]*(1-faTc[x][y])-392*faTc[x][y];
   dfaTmgc=200000*Cai[x][y]*(1-faTmgc[x][y]-faTmgm[x][y])-6.6*faTmgc[x][y];
   dfaTmgm=2000*2.5*(1-faTmgc[x][y]-faTmgm[x][y])-666*faTmgm[x][y];
   dfaCalse=480*Carel[x][y]*(1-faCalse[x][y])-400*faCalse[x][y];       
   dfab=0.08*dfaTc+0.16*dfaTmgc+0.045*dfac;

   Caup[x][y]+=HT*(Iup-Itr)/(2.*Faraday*Vup);
   Carel[x][y]+=HT*((Itr-Irel)/(2.*Faraday*Vrel)-31.0*(480*Carel[x][y]*(1-faCalse[x][y])-400*faCalse[x][y]));   
   Cai[x][y]+=HT*((2.*INaCa-ICaL-ICaT-ICap-IbCa-Iup+Irel)/(2.*VCa*Faraday)-dfab);

   fac[x][y]+=HT*dfac;
   faTc[x][y]+=HT*dfaTc;
   faTmgc[x][y]+=HT*dfaTmgc;
   faTmgm[x][y]+=HT*dfaTmgm;
   faCalse[x][y]+=HT*dfaCalse;

   Kc[x][y]+=HT*(-2.*Ip+IKf+IKs+Ito+Ik1)/(Faraday*Vc);
   Ki[x][y]+=HT*(2.*Ip-IKf-IKs-Ito-Ik1)/(Faraday*Vi);

   ract=240.0*exp((V[x][y]-20.)/12.5)+203.8*pow(Cai[x][y]/(Cai[x][y]+0.0003),4); 
   rinact=33.96+339.6*pow(Cai[x][y]/(Cai[x][y]+0.0003),4);

   dF1[x][y]=0.815*F3[x][y]-ract*F1[x][y];
   dF2[x][y]=ract*F1[x][y]-rinact*F2[x][y];
   dF3[x][y]=rinact*F2[x][y]-0.815*F3[x][y];

   F1[x][y]+=HT*dF1[x][y];
   F2[x][y]+=HT*dF2[x][y];
   F3[x][y]+=HT*dF3[x][y];


   return ((IKf+Ik1+Ito+Ip+INaCa+INa+IbNa+ICap+ICaL+ICaT+IbCa+IKH)/Cm);

}
