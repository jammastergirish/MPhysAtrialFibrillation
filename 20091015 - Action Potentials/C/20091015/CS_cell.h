#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Atria_Variables.h"

float I_tot_CS ()
{
   float IKf,IKs,IK,Ik1,Ito,IbNa,Ip,INaCa,IbCa,INa, ICap,ICaL,ICaT, /*Iup,Itr,Irel, */ Isus;

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
   float dfca;
   float Adt, Bdt;
   float dtm, tdt;
   float Aft, Bft;
   float ftm, tft;
   float E0=V-Ek+3.6;

   IKf=3.5*Pa*Pi*(V-Ek);     

   Apa=9.0*exp(V/25.371);
   Bpa=1.3*exp(-V/13.026);
   pam=1./(1+exp(-(V+5.1)/7.4));
   tpa=1./(Apa+Bpa);
   Pa+=HT*(pam-Pa)/tpa;
      
   Api=100.*exp(-V/54.645);
   Bpi=656.*exp(V/106.157);
   pim=1./(1+exp((V+47.3921)/18.6603));
   tpi=1.0/(Api+Bpi);
   Pi+=HT*(pim-Pi)/tpi;

   IKs=2.5*n*(V-Ek);      

   An=1.66*exp(V/69.452);
   Bn=0.3*exp(-V/21.826);
   nm=1.0/(1+exp(-(V-0.9)/13.8));
   tn=1./(An+Bn)+0.06;   
   n+=HT*(nm-n)/tn;

   IK=IKf+IKs;
   
   Ik1=(5.08*pow(5.0/(5.0+0.59),3)/(1+exp(1.393*(V-Ek-12)/RTONF))+0.05)*(V-Ek-12);

   Ito=0.31*50.02*r1*(0.590*pow(s1,3)+0.410*pow(s2,3))*(0.600*pow(s3,6)+0.4)*(V-Ek); 
   Isus = 0.4*(V+85)+0.2*(V-5);
   Ito+=Isus;
   
   Ar=386.6*exp(V/12.0);
   Br=8.011*exp(-V/7.2);
   rm=1/(1+exp(-(V+15.0)/5.633));   
   tr=1./(Ar+Br)+0.0004;
   r1+=HT*(rm-r1)/tr;

   s1m=1./(1+exp((V+24.8)/5.6));
   ts1=0.0496/(1+exp((V+32.8)/0.1))+0.0158; // t1_r-t1_d / (Exp) + t1_d
   s1+=HT*(s1m-s1)/ts1;

   s2m= 1./(1+exp((V+24.8)/5.6)); 
   ts2=2.9536/(1+exp((V+32.8)/0.1))+0.1245; // t2_r+t2_d / (Exp) + t2_d
   s2+=HT*(s2m-s2)/ts2;
    
   s3m=((1./(1+exp((V+50.67)/27.38)))+0.666)/1.666;
   ts3=(7.5/(1+exp((V+23.0)/0.5)))+0.5;      
   s3+=HT*(s3m-s3)/ts3;

   Ip=64.41*Kc/(Kc+1)*(pow(Nai,1.5)/(pow(Nai,1.5)+pow(11,1.5)))*(1.6/(1.5+exp(-(V+60)/40.)));

   INaCa=0.02*(pow(Nai,3)*2.5*exp(0.450*V/RTONF)-pow(140,3)*Cai*exp(V*(0.45-1)/RTONF))/(1+0.0003*(Cai*pow(140,3)+2.5*pow(Nai,3)));  

   if(fabs(V+44.4) < 0.0001)
      Am=460.*12.673;
   else    
     Am=-460*(V+44.4)/(exp(-(V+44.4)/12.673)-1);
   Bm=18400.0*exp(-(V+44.4)/12.673);
   m = Am/(Am+Bm);
//   m+=HT*(Am*(1-m)-Bm*m);

   Ah=44.9*exp(-(V+66.9)/5.57);
   Bh=1491.0/(1+323.3*exp(-(V+94.6)/12.9));

   th1=0.03/(1+exp((V+40)/6.0))+0.00015;  //0.00035 - Lindblad   
   th2=0.12/(1+exp((V+60)/2.0))+0.00045;  //0.00295 - Lindblab    

   hm=Ah/(Ah+Bh);
   h1+=HT*(hm-h1)/th1; 
   h2+=HT*(hm-h2)/th2;

   if(fabs(V) > 0.0001) 
      INa=0.0014*pow(m,3)*(0.635*h1+0.365*h2)*140*V*(Faraday/RTONF)*(exp((V-ENa)/RTONF)-1)/(exp(V/RTONF)-1); // multiply by 0.75 for instant activation
   else INa=0.0014*pow(m,3)*(0.635*h1+0.365*h2)*140*Faraday*(exp((V-ENa)/RTONF)-1);   

   ICap=9.509*(Cai/(Cai+0.0002));   
   
   Adl=-16.72*(V+45)/(exp(-(V+45)/2.5)-1)-50.0*(V+10)/(exp(-(V+10)/4.808)-1);
   if(fabs(V+45) < 0.0001)
     Adl=16.72*2.5-50.0*(V+10)/(exp(-(V+10)/4.808)-1);
   if(fabs(V+10) < 0.0001) 
     Adl=-16.72*(V+45)/(exp(-(V+45)/2.5)-1)+50.0*4.808;
       
   if(fabs(V+5.) < 0.0001)
     Bdl=4.48*2.5;
   else
     Bdl=4.48*(V+5)/(exp((V+5)/2.5)-1.);
   tdl=1/(Adl+Bdl);

   dlm=1./(1+exp(-(V-0.5)/5.7));
   dL+=HT*(dlm-dL)/tdl;


  flm=1./(1+exp((V+17.9)/5.2));
  tfl1=(0.0335/(1+exp((V+18)/4.0)))+0.0117;
  tfl2=(2.2989/(1+exp((V+18)/4.0)))+0.0478;
  fL1+=HT*(flm-fL1)/tfl1;
  fL2+=HT*(flm-fL2)/tfl2;
/*
  dfca=1.0/(1+Cai/0.00015);
  fca+=HT*(dfca-fca)/0.02;
*/
  ICaL=0.6*29.7*(dL*((0.80*fL1)+(0.2*fL2)))*(V-58.3);   

   Adt=674.173*exp((V+23.3)/30.);
   Bdt=674.173*exp(-(V+23.3)/30.);
   tdt=1/(Adt+Bdt);
   dtm=1./(1+exp(-(V+23.0)/6.1));
   dT+=HT*(dtm-dT)/tdt;

   Aft=9.637*exp(-(V+75)/83.3);
   Bft=9.637*exp((V+75)/15.38);
   tft=1.0/(Aft+Bft);
   ftm=Aft/(Aft+Bft);
   fT+=HT*(ftm-fT)/tft;

   ICaT=6*dT*fT*(V-38.0);

   IbNa=0.02*(V-ENa);  // 0.02 - CT, 0.03 - PM

   IbCa=0.02*(V-ECa);  // 0.02 - CT, 0.03 - PM
/*
   Nai+=HT*(-3*Ip-3*INaCa-IbNa-INa)/(Faraday*Vi);

  Iup=2800.0*(((Cai/0.0003)-(0.4*0.4*Caup/0.5))/(((Cai+0.0003)/0.0003)+(0.4*(Caup+0.5)/0.5)));

  Irel=200000.0*pow(F2/(F2+0.25),2)*(Carel-Cai);    

  Itr=(Caup-Carel)*2.*Faraday*Vup/0.01;

   dfac=200000.0*Cai*(1-fac)-476*fac;
   dfaTc=78400*Cai*(1-faTc)-392*faTc;
   dfaTmgc=200000*Cai*(1-faTmgc-faTmgm)-6.6*faTmgc;
   dfaTmgm=2000*2.5*(1-faTmgc-faTmgm)-666*faTmgm;
   dfaCalse=480*Carel*(1-faCalse)-400*faCalse;       
   dfab=0.08*dfaTc+0.16*dfaTmgc+0.045*dfac;

   Caup+=HT*(Iup-Itr)/(2.*Faraday*Vup);
   Carel+=HT*((Itr-Irel)/(2.*Faraday*Vrel)-31.0*(480*Carel*(1-faCalse)-400*faCalse));   
   Cai+=HT*((2.*INaCa-ICaL-ICaT-ICap-IbCa-Iup+Irel)/(2.*VCa*Faraday)-dfab);

   fac+=HT*dfac;
   faTc+=HT*dfaTc;
   faTmgc+=HT*dfaTmgc;
   faTmgm+=HT*dfaTmgm;
   faCalse+=HT*dfaCalse;

   Kc+=HT*(-2.*Ip+IKf+IKs+Ito+Ik1)/(Faraday*Vc);
   Ki+=HT*(2.*Ip-IKf-IKs-Ito-Ik1)/(Faraday*Vi);

   ract=240.0*exp((V-20.)/12.5)+203.8*pow(Cai/(Cai+0.0003),4); 
   rinact=33.96+339.6*pow(Cai/(Cai+0.0003),4);

   dF1=0.815*F3-ract*F1;
   dF2=ract*F1-rinact*F2;
   dF3=rinact*F2-0.815*F3;

   F1+=HT*dF1;
   F2+=HT*dF2;
   F3+=HT*dF3;
*/

   return ((IKf+IKs+Ik1+Ito+Ip+INaCa+INa+IbNa+ICap+ICaL+ICaT+IbCa)/Cm);

}
