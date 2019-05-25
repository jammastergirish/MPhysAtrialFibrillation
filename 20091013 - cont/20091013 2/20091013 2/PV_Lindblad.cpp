 // Rabbit atrial crista terminalis (CT) and pectinate muscle (PM) cell models
// Created by Oleg Aslanidi, University of Manchester, April 2008

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "Atria_Variables.h"
//#include "PV_cell.h"

#define BCL 10
#define beats 10

float IKf,IKs,IK,Ik1,Ito,IbNa,Ip,INaCa,IbCa,INa, ICap,ICaL,ICaT,IbK, Itot,Iup,Itr,Irel,Isus, IfNa, IfK, If, IKH, IKH_K, IKH_Na;
float Vcell,Vi,VCa,Vc,Vup,Vrel;
float Ek,ENa,ECa,ECl,EbCl;
float ract,rinact,dF1,dF2,dF3;
float dfac,dfaTc,dfaTmgc,dfaTmgm,dfaCalse,dfab;

float HT=0.000005;    
float TMAX=2.0;

float V=-80.0;
float Pa=0.00016;
float Pi=0.76898;
float n=0.02032;
float r1=0.00006;
float s1=0.5753;
float s2=0.39871;
float s3=0.57363;
float r2=0.064;
float m=0.01309;
float h1=0.706;
float h2=0.61493;
float dL=0.00003;
float fL2=0.99981;
float fL1=0.9981;
float fL=0.9981;
float dT=0.00046;
float fT=0.30752;
float Nai=8.4;
float Caup=0.730866;
float Carel=0.726776;
float Cai=0.000071; 
float fac=0.029108;
float faTc=0.014071;
float faTmgc=0.214036;
float faTmgm=0.693565;
float faCalse=0.465921;
float Kc=5.0;
float Ki=140.0;    // 100.0 - Lindblad
float F1=0.288039;
float F2=0.002262;
float F3=0.612697;
float fca=0.7755;
float xsus=0.01;
float y2=0.0927; //0.0184

const float Temp=308;
const float Cm=0.05;
const float Faraday=96487;
const float R=8314;
float RTONF;

void comp_ikh()
{
   float Ay, By;
   float ym, ty;
   
  Ay = exp(-(2.9 + (0.04*V)));
  By  = exp((3.6 + (0.11*V)));

   ym = 1/(1 + exp((V+100)/20));
      
  ty = 0.4/(Ay+By);

  y2 += HT*(ym-y2)/ty;

  IKH_K= 0.1*ym*(V-Ek);
  IKH_Na = 0.008*ym*(V-ENa);
  IKH = (IKH_K + IKH_Na);
   
}


/*
void comp_if()
{
   float Ay, By;
   float ym, ty;
   
   Ay = exp(-(2.9 + (0.04*V)));
   By  = exp((3.6 + (0.11*V)));

   ym = 1/(1 + exp((V+100)/10));
   ty = 0.4/(Ay+By);

   y2 += HT*(ym-y2)/ty;

   IfNa = 0.00072827*y2*(V-ENa); //0.72827
   IfK = 0.00117173*y2*(V-Ek); //1.17173
   If = 4400*(IfNa + IfK); //10*
}*/

void comp_ikf ()
{
   float Apa, Bpa;
   float pam, tpa;
   float Api, Bpi;
   float pim, tpi;

   IKf=4*0.0035*Pa*Pi*(V-Ek);    

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
}

void comp_iks ()
{
   float An, Bn;
   float nm, tn;

   IKs=2.5*n*(V-Ek);      

   An=1.66*exp(V/69.452);
   Bn=0.3*exp(-V/21.826);
   nm=1.0/(1+exp(-(V-0.9)/13.8));
   tn=1./(An+Bn)+0.06;   
   n+=HT*(nm-n)/tn;
}

void comp_ik ()
{
   IK=IKf+IKs;
}

void comp_ik1 ()
{
   float E0=V-Ek+3.6;

//   Ik1=0.52*5.08*pow(Kc/(Kc+0.59),3)*(V-Ek)/(1+exp(1.393*E0/RTONF));     // CT - 2.0, PM - 2.5
     Ik1=(0.9*5.08*pow(5.0/(5.0+0.59),3)/(1+exp(1.393*(V-Ek-12)/RTONF))+0.05)*(V-Ek-12);

}

void comp_ito ()
{
   float Ar, Br;
   float rm, tr;
   float s1m, ts1;
   float s2m, ts2;
   float s3m, ts3;
   float xsusm;

   Ito=0.45*0.515*50.02*r1*(0.8*pow(s1,3)+0.2*pow(s2,3))*/*(0.600*pow(s3,6)+0.4)*/(V-Ek);     // CT - 0.2, PM - 0.35
//   Isus = 0.6*(V+70);  // 1.4, 2.4
  Isus = 0.47*(V+85)+0.0*(V-5);
   
/*   xsusm = 1.0/(1 + exp(-(V+25)/5));
   xsus+=HT*(xsusm-xsus)/0.001;
   Isus = 0.05*xsus*(V + 80);
*/   
   Ito+=Isus;
   
   Ar=386.6*exp(V/12.0);
   Br=8.011*exp(-V/7.2);
   rm=1/(1+exp(-(V+15.0)/5.633));   
   tr=1./(Ar+Br)+0.0004;
   r1+=HT*(rm-r1)/tr;

   s1m=1./(1+exp((V+29.4)/5.1));
   ts1=0.1155/(1+exp((V+32.8)/0.1))+0.0147; // t1_r-t1_d / (Exp) + t1_d
   s1+=HT*(s1m-s1)/ts1;

   s2m= 1./(1+exp((V+29.4)/5.1)); 
   ts2=3.0967/(1+exp((V+32.8)/0.1))+0.0926; // t2_r+t2_d / (Exp) + t2_d
   s2+=HT*(s2m-s2)/ts2;
    
   s3m=((1./(1+exp((V+50.67)/27.38)))+0.666)/1.666;
   ts3=(7.5/(1+exp((V+23.0)/0.5)))+0.5;      
   s3+=HT*(s3m-s3)/ts3;
}


void comp_ip ()
{
   Ip=64.41*Kc/(Kc+1)*(pow(Nai,1.5)/(pow(Nai,1.5)+pow(11,1.5)))*(1.6/(1.5+exp(-(V+60)/40.)));
}

void comp_inaca ()
{
//Regular AP
//INaCa=0.02*(pow(Nai,3)*2.5*exp(0.450*V/RTONF)-pow(140,3)*Cai*exp(V*(0.45-1)/RTONF))/(1+0.0003*(Cai*pow(140,3)+2.5*pow(Nai,3)));  
   
//SAP
INaCa = 1.5*(0.02*(pow(Nai,3.0)*7.5*exp(0.35*V/RTONF)-pow(140,3.0)*Cai*exp(V*(0.45-1)/RTONF))/(1+0.0003*(Cai*pow(140,3.0)+2.5*pow(Nai,3.0)))); //1.5

}

void comp_ina ()
{
   float Am, Bm;
   float Ah, Bh;
   float hm;
   float th1, th2;

   if(fabs(V+44.4) < 0.0001)
      Am=460.*12.673;
   else    
     Am=-460*(V+44.4)/(exp(-(V+44.4)/12.673)-1);
   Bm=18400.0*exp(-(V+44.4)/12.673);
//   m = Am/(Am+Bm);
   m+=HT*(Am*(1-m)-Bm*m);

   Ah=44.9*exp(-(V+66.9)/5.57);
   Bh=1491.0/(1+323.3*exp(-(V+94.6)/12.9));

   th1=0.03/(1+exp((V+40)/6.0))+0.00015;  //0.00035 - Lindblad   
   th2=0.12/(1+exp((V+60)/2.0))+0.00045;  //0.00295 - Lindblab    

   hm=Ah/(Ah+Bh);
   h1+=HT*(hm-h1)/th1; 
   h2+=HT*(hm-h2)/th2;

   if(fabs(V) > 0.0001) 
      INa=0.0014*pow(m,3)*(0.635*h1+0.365*h2)*140*V*(Faraday/RTONF)*(exp((V-ENa)/RTONF)-1)/(exp(V/RTONF)-1); // multiply by 0.8 for instant activation
   else INa=0.0014*pow(m,3)*(0.635*h1+0.365*h2)*140*Faraday*(exp((V-ENa)/RTONF)-1);   
}

void comp_icap ()
{
   ICap=9.509*(Cai/(Cai+0.0002));   
}

void comp_ical ()
{
   float Adl, Bdl;
   float dlm, tdl;
   float Afl, Bfl;
   float flm, tfl1, tfl2;
   float dfca;
   
   Adl=-16.72*(V+45)/(exp(-(V+45)/2.5)-1)-50.0*(V+10)/(exp(-(V+10)/4.808)-1);
   if(fabs(V+35) < 0.0001)
     Adl=16.72*2.5-50.0*(V+10)/(exp(-(V+10)/4.808)-1);
   if(fabs(V+10) < 0.0001) 
     Adl=-16.72*(V+45)/(exp(-(V+45)/2.5)-1)+50.0*4.808;
       
   if(fabs(V+5.) < 0.0001)
     Bdl=4.48*2.5;
   else
     Bdl=4.48*(V+5)/(exp((V+5)/2.5)-1.);
   tdl=1/(Adl+Bdl);

   dlm=1./(1+exp(-(V+7.5)/5.5));
   dL+=HT*(dlm-dL)/tdl;


  flm=1./(1+exp((V+23.4)/5.2));
  tfl1=(0.0336/(1+exp((V+18)/4.0)))+0.0142;
  tfl2=(3.6908/(1+exp((V+18)/4)))+0.0534;
  fL1+=HT*(flm-fL1)/tfl1;
  fL2+=HT*(flm-fL2)/tfl2;

  dfca=1.0/(1+Cai/0.0001); //was 0.0001, must be higher than 0.000071 Cai(?)
  fca+=HT*(dfca-fca)/0.02;

//   ICaL=18*(dL*((0.80*fL1)+(0.2*fL2)))*fca*(V-61.4);   
     ICaL=25*dL*fL1*fca*(V-58.3);
 
}

void comp_icat ()
{
   float Adt, Bdt;
   float dtm, tdt;
   float Aft, Bft;
   float ftm, tft;

   Adt=674.173*exp((V+23.3)/30.0);
   Bdt=674.173*exp(-(V+23.3)/30.0); //23.3 
   tdt=1/(Adt+Bdt);
   dtm=1./(1+exp(-(V+23.0)/5.1)); //23 6.1
   dT+=HT*(dtm-dT)/tdt;

   Aft=9.637*exp(-(V+75)/83.3);// 75
   Bft=9.637*exp((V+75)/16.3);// 15.38
   tft=1.0/(Aft+Bft);
   ftm=Aft/(Aft+Bft);
   fT+=HT*(ftm-fT)/tft;

   ICaT=1.65*8.0*dT*fT*(V-38.0);
}

void comp_ibna ()
{ 
  IbNa=0.02*(V-ENa);  // 0.02 - CT, 0.03 - PM
}

void comp_ibca ()
{
  IbCa=0.02*(V-ECa);  // 0.02 - CT, 0.03 - PM
}

void comp_itot ()
{
   Itot=IKf+Ik1+Ito+Ip+INaCa+INa+IbNa+ICap+ICaL+ICaT+IbCa+IKH;//If;
}

void comp_nai ()
{
   Nai+=HT*(-3*Ip-3*INaCa-IbNa-INa)/(Faraday*Vi);
}

void comp_iup ()
{
   Iup=2*2800.0*(((Cai/0.0003)-(0.4*0.4*Caup/0.5))/(((Cai+0.0003)/0.0003)+(0.4*(Caup+0.5)/0.5)));
}    
 
void comp_irel ()
{
   Irel=2*200000.0*pow(F2/(F2+0.25),2)*(Carel-Cai);    
}

void comp_itr ()
{
   Itr=(Caup-Carel)*2.*Faraday*Vup/0.01;
}

void comp_ca ()
{   
   dfac=200000.0*Cai*(1-fac)-476*fac;
   dfaTc=78400*Cai*(1-faTc)-392*faTc;
   dfaTmgc=200000*Cai*(1-faTmgc-faTmgm)-6.6*faTmgc;
   dfaTmgm=2000*2.5*(1-faTmgc-faTmgm)-666*faTmgm;
   dfaCalse=(480*Carel*(1-faCalse)-400*faCalse);       
   dfab=0.08*dfaTc+0.16*dfaTmgc+0.045*dfac;

   Caup+=HT*(Iup-Itr)/(2.*Faraday*Vup);
   Carel+=HT*((Itr-Irel)/(2.*Faraday*Vrel)-31.0*(480*Carel*(1-faCalse)-400*faCalse));   
   Cai+=HT*((2.*INaCa-ICaL-ICaT-ICap-IbCa-Iup+Irel)/(2.*VCa*Faraday)-dfab);

   fac+=HT*dfac;
   faTc+=HT*dfaTc;
   faTmgc+=HT*dfaTmgc;
   faTmgm+=HT*dfaTmgm;
   faCalse+=HT*dfaCalse;
}

void comp_ki ()
{
   Kc+=HT*(-2.*Ip+IKf+IKs+Ito+Ik1)/(Faraday*Vc);
   Ki+=HT*(2.*Ip-IKf-IKs-Ito-Ik1)/(Faraday*Vi);
}

void comp_act ()
{
   ract=240.0*exp((V-20.)/12.5)+203.8*pow(Cai/(Cai+0.0003),4); 
   rinact=33.96+339.6*pow(Cai/(Cai+0.0003),4);

   dF1=0.815*F3-ract*F1;
   dF2=ract*F1-rinact*F2;
   dF3=rinact*F2-0.815*F3;

   F1+=HT*dF1;
   F2+=HT*dF2;
   F3+=HT*dF3;
}
   

int main ()
{
   float t=0.0;
   float Istim;
   float dvdt;
   int cnt=0;
   int cnt2=0;
   int cnt3=0;
   int i=0;

   float tstim=0.01;
   float stimtime;

   float vmax[beats+1];
   float dvdtmax[beats+1];
   float apd[beats+1];
   float toneapd[beats+1];
   float ttwoapd[beats+1];
   float rmbp[beats+1];

   FILE *in;
/*
   Vcell=0.003497;  
   Vi=0.465*Vcell;
   VCa=Vi;
   Vc=0.1*Vcell;
   Vup=0.0283*Vcell;
   Vrel=0.00315*Vcell;
*/
   Vi=0.0126;
   VCa=0.005884;
   Vc=0.0025;
   Vup=0.0003969;
   Vrel=0.000044;
//
   RTONF=Temp*0.08554;

   in = fopen ("PV.txt", "w");

   TMAX=BCL*beats+tstim;
   for (i=beats; i>0; i--) {
     vmax[i] = -100;
	 dvdtmax[i] = 0;
   }

   while (t < TMAX) {
	 t+=HT;
     
	 Ek=RTONF*log(Kc/Ki);
     ENa=RTONF*log(140./Nai);     
     ECa=(RTONF/2.)*log(2.5/Cai);  
     ECl=RTONF*exp(30./132.); 
     EbCl=ECl-0.49*(ECl+30.59);

     comp_ikf ();
//	 comp_iks ();
	 comp_ik ();
	 comp_ik1 ();
	 comp_ito ();
	 comp_ina ();
	 comp_ibna ();
	 comp_ip ();
	 comp_inaca ();
	 comp_ical ();
	 comp_icat ();
	 comp_icap ();
	 comp_ibca ();
	 comp_itot ();
     comp_ikh();
//comment this out to switch off Ca2+ handling
     comp_nai ();
     comp_ki ();
     comp_iup ();
     comp_irel ();
     comp_itr ();
     comp_act ();
     comp_ca ();


     if (t>=tstim && t<(tstim+HT)) {
       stimtime = 0;
       i = i+1;

       tstim = tstim+BCL;
       printf ("Stimulus %d applied, ", i);
       printf ("Time = %f\n", t);
       rmbp[i]=V;
     }

     if(stimtime>=0 && stimtime<0.0002)    
       Istim=-15000.0; // 15000
     else Istim=0.0;
      
     stimtime = stimtime+HT;
     
	 dvdt=-(Itot+Istim)/Cm;

	 V+=HT*dvdt;

	 if (V>vmax[i])
       vmax[i] = V;
     if (dvdt>dvdtmax[i] && t > (tstim-BCL+0.0003)) {
       dvdtmax[i] = dvdt;
       toneapd[i] = t;
	 }
     if (V>=(vmax[i]-0.9*(vmax[i]-rmbp[i])))
       ttwoapd[i] = t;

	 cnt++;
	 cnt2++;
	 cnt3++;
     if (cnt == 200) { 
       fprintf (in, "%lf\t%lf\t%lf\n", t, V, 1000*Cai);
	   cnt = 0; 
        if(cnt2==200) {printf("."); cnt2=0;}
        if(cnt3 == 300000) {      printf (" Working on Stimulus %d, ", i);
       printf ("Time = %f ", t); cnt3 = 0;}
       
	 }
   } 
   fclose (in);

   apd[beats]=ttwoapd[beats]-toneapd[beats];
   printf ("APD = %8.2lf, dVdtmax = %8.2lf\n", 1000*apd[beats], 0.001*dvdtmax[beats]);
   
   system("pause");

   return 0;
}


