// Rabbit atrial crista terminalis (CT) and pectinate muscle (PM) cell models
// Created by Oleg Aslanidi, University of Manchester, April 2008

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctime>

#define BCL 10
#define beats 1

double IKH, IKH_K, IKH_Na,IKf,IKs,IK,Ik1,Ito,IbNa,Ip,INaCa,IbCa,INa, ICap,ICaL,ICaT,IbK, Itot,Iup,Itr,Irel,Isus, IfNa, IfK, If, IClCa_ss, IClCa, IbNa_b, IbK_b, IbCl_b, Isus_b, IClCa_jss;
double Vcell,Vi,VCa,Vc,Vup,Vrel;
double Ek,ENa,ECa,ECl,EbCl;
double ract,rinact,dF1,dF2,dF3;
double dfac,dfaTc,dfaTmgc,dfaTmgm,dfaCalse,dfab;

double HT=0.000005;    
double TMAX=2.0;

double V=-80.0;
double Pa=0.00016;
double Pi=0.76898;
double n=0.02032;
double r1=0.00006;
double s1=0.5753;
double s2=0.39871;
double s3=0.57363;
double r2=0.064;
double m=0.01309;
double h1=0.706;
double h2=0.61493;
double dL=0.00003;
double fL3=0.999981;
double fL2=0.99981;
double fL1=0.9981;
double fL=0.9981;
double dT=0.00046;
double fT=0.30752;
double Nai=8.4;
double fac=0.029108;
double Kc=5.0;
double Ki=140.0;    // 100.0 - Lindblad
double F1=0.288039;
double F2=0.002262;
double F3=0.612697;
double fca=0.7755;
double xsus=0.01;
double y2=0.0927; //0.0184
double Kf=13315.6/20;		
double Kb=2000.0/20;
double ss1=0.021468; //0.01309
double ss2=0.032789;


const double Temp=308;
const double Cm=0.05;
const double Faraday=96487;
const double R=8314;
double RTONF;

//------------------------------------------------------------------------------
// Constants
//------------------------------------------------------------------------------

float Ca_buffer___Bmax_Calsequestrin;   // millimolar
float Ca_buffer___Bmax_SLB_SL;   // millimolar
float Ca_buffer___Bmax_SLB_jct;   // millimolar
float Ca_buffer___Bmax_SLHigh_SL;   // millimolar
float Ca_buffer___Bmax_SLHigh_jct;   // millimolar
float Ca_buffer___koff_Calsequestrin;   // per_millisecond
float Ca_buffer___koff_SLB;   // per_millisecond
float Ca_buffer___koff_SLHigh;   // per_millisecond
float Ca_buffer___kon_Calsequestrin;   // per_millimolar_per_millisecond
float Ca_buffer___kon_SL;   // per_millimolar_per_millisecond
float Jleak_SR___KSRleak;   // per_millisecond
float Jpump_SR___H;   // dimensionless
float Jpump_SR___Kmf;   // millimolar
float Jpump_SR___Kmr;   // millimolar
float Jpump_SR___Q10_SRCaP;   // dimensionless
float Jpump_SR___V_max;   // millimolar_per_millisecond
float Jrel_SR___EC50_SR;   // millimolar
float Jrel_SR___HSR;   // dimensionless
float Jrel_SR___Max_SR;   // dimensionless
float Jrel_SR___Min_SR;   // dimensionless
float Jrel_SR___kiCa;   // per_millimolar_per_millisecond
float Jrel_SR___kim;   // per_millisecond
float Jrel_SR___koCa;   // per_millimolar2_per_millisecond
float Jrel_SR___kom;   // per_millisecond
float Jrel_SR___ks;   // per_millisecond
float cytosolic_Ca_buffer___Bmax_Calmodulin;   // millimolar
float cytosolic_Ca_buffer___Bmax_Myosin_Ca;   // millimolar
float cytosolic_Ca_buffer___Bmax_Myosin_Mg;   // millimolar
float cytosolic_Ca_buffer___Bmax_SRB;   // millimolar
float cytosolic_Ca_buffer___Bmax_TroponinC;   // millimolar
float cytosolic_Ca_buffer___Bmax_TroponinC_Ca_Mg_Ca;   // millimolar
float cytosolic_Ca_buffer___Bmax_TroponinC_Ca_Mg_Mg;   // millimolar
float cytosolic_Ca_buffer___koff_Calmodulin;   // per_millisecond
float cytosolic_Ca_buffer___koff_Myosin_Ca;   // per_millisecond
float cytosolic_Ca_buffer___koff_Myosin_Mg;   // per_millisecond
float cytosolic_Ca_buffer___koff_SRB;   // per_millisecond
float cytosolic_Ca_buffer___koff_TroponinC;   // per_millisecond
float cytosolic_Ca_buffer___koff_TroponinC_Ca_Mg_Ca;   // per_millisecond
float cytosolic_Ca_buffer___koff_TroponinC_Ca_Mg_Mg;   // per_millisecond
float cytosolic_Ca_buffer___kon_Calmodulin;   // per_millimolar_per_millisecond
float cytosolic_Ca_buffer___kon_Myosin_Ca;   // per_millimolar_per_millisecond
float cytosolic_Ca_buffer___kon_Myosin_Mg;   // per_millimolar_per_millisecond
float cytosolic_Ca_buffer___kon_SRB;   // per_millimolar_per_millisecond
float cytosolic_Ca_buffer___kon_TroponinC;   // per_millimolar_per_millisecond
float cytosolic_Ca_buffer___kon_TroponinC_Ca_Mg_Ca;   // per_millimolar_per_millisecond
float cytosolic_Ca_buffer___kon_TroponinC_Ca_Mg_Mg;   // per_millimolar_per_millisecond
float model_parameters___F;   // coulomb_per_mole
float model_parameters___R;   // joule_per_kilomole_kelvin
float model_parameters___T;   // kelvin
float model_parameters___cell_length;   // micrometre
float model_parameters___cell_radius;   // micrometre

//------------------------------------------------------------------------------
// Computed variables
//------------------------------------------------------------------------------

float Ca_buffer___dCa_SLB_SL;   // millimolar_per_millisecond
float Ca_buffer___dCa_SLB_jct;   // millimolar_per_millisecond
float Ca_buffer___dCa_SLHigh_SL;   // millimolar_per_millisecond
float Ca_buffer___dCa_SLHigh_jct;   // millimolar_per_millisecond
float Ca_buffer___dCa_SL_tot_bound;   // millimolar_per_millisecond
float Ca_buffer___dCa_jct_tot_bound;   // millimolar_per_millisecond
float Ca_buffer___dCalsequestrin;   // millimolar_per_millisecond
float Ca_buffer___i_Ca_SL_tot;   // microA_per_microF
float Ca_buffer___i_Ca_jct_tot;   // microA_per_microF
float Jleak_SR___j_leak_SR;   // millimolar_per_millisecond
float Jpump_SR___Q_SRCaP;   // dimensionless
float Jpump_SR___j_pump_SR;   // millimolar_per_millisecond
float Jrel_SR___RI;   // dimensionless
float Jrel_SR___j_rel_SR;   // millimolar_per_millisecond
float Jrel_SR___kCaSR;   // dimensionless
float Jrel_SR___kiSRCa;   // per_millimolar_per_millisecond
float Jrel_SR___koSRCa;   // per_millimolar2_per_millisecond
float cytosolic_Ca_buffer___dCa_Calmodulin;   // millimolar_per_millisecond
float cytosolic_Ca_buffer___dCa_Myosin;   // millimolar_per_millisecond
float cytosolic_Ca_buffer___dCa_SRB;   // millimolar_per_millisecond
float cytosolic_Ca_buffer___dCa_TroponinC;   // millimolar_per_millisecond
float cytosolic_Ca_buffer___dCa_TroponinC_Ca_Mg;   // millimolar_per_millisecond
float cytosolic_Ca_buffer___dCa_cytosol_tot_bound;   // millimolar_per_millisecond
float cytosolic_Ca_buffer___dMg_Myosin;   // millimolar_per_millisecond
float cytosolic_Ca_buffer___dMg_TroponinC_Ca_Mg;   // millimolar_per_millisecond
float ion_diffusion___J_Ca_SL_cytosol;   // millimole_per_millisecond
float ion_diffusion___J_Ca_SL_cytosol_old;   // millimole_per_millisecond
float ion_diffusion___J_Ca_jct_SL;   // millimole_per_millisecond
float ion_diffusion___J_Ca_jct_SL_old;   // millimole_per_millisecond
float model_parameters___Vol_Cell;   // litre
float model_parameters___Vol_SL;   // litre
float model_parameters___Vol_SR;   // litre
float model_parameters___Vol_cytosol;   // litre
float model_parameters___Vol_jct;   // litre

//------------------------------------------------------------------------------
// State variables
//------------------------------------------------------------------------------

float Y[19];
float dY[19];

//------------------------------------------------------------------------------
// Initialisation
//------------------------------------------------------------------------------

void init()
{
	   //---------------------------------------------------------------------------
   // State variables
   //---------------------------------------------------------------------------

   Y[0] = 1.121211963;   // Ca_buffer___Ca_Calsequestrin (millimolar)
   Y[1] = 0.000071563;   // Ca_buffer___Ca_SL (millimolar)
   Y[2] = 0.006662134;   // Ca_buffer___Ca_SLB_SL (millimolar)
   Y[3] = 0.005341615;   // Ca_buffer___Ca_SLB_jct (millimolar)
   Y[4] = 0.087317089;   // Ca_buffer___Ca_SLHigh_SL (millimolar)
   Y[5] = 0.061253947;   // Ca_buffer___Ca_SLHigh_jct (millimolar)
   Y[6] = 0.492827646;   // Ca_buffer___Ca_SR (millimolar)
   Y[7] = 0.000119433;   // Ca_buffer___Ca_jct (millimolar)
   Y[8] = 0.000071642;   // Ca_buffer___Cai (millimolar)
   Y[9] = 0.000000077;   // Jrel_SR___I (dimensionless)
   Y[10] = 0.000000267;   // Jrel_SR___O (dimensionless)
   Y[11] = 0.775155733;   // Jrel_SR___R (dimensionless)
   Y[12] = 0.000243254;   // cytosolic_Ca_buffer___Ca_Calmodulin (millimolar)
   Y[13] = 0.001615156;   // cytosolic_Ca_buffer___Ca_Myosin (millimolar)
   Y[14] = 0.001826852;   // cytosolic_Ca_buffer___Ca_SRB (millimolar)
   Y[15] = 0.007520170;   // cytosolic_Ca_buffer___Ca_TroponinC (millimolar)
   Y[16] = 0.115353534;   // cytosolic_Ca_buffer___Ca_TroponinC_Ca_Mg (millimolar)
   Y[17] = 0.137869474;   // cytosolic_Ca_buffer___Mg_Myosin (millimolar)
   Y[18] = 0.011586448;   // cytosolic_Ca_buffer___Mg_TroponinC_Ca_Mg (millimolar)

   //---------------------------------------------------------------------------
   // Constants
   //---------------------------------------------------------------------------

   Ca_buffer___Bmax_Calsequestrin = 0.14;   // millimolar
   Ca_buffer___Bmax_SLB_SL = 0.0374;   // millimolar
   Ca_buffer___Bmax_SLB_jct = 0.0046;   // millimolar
   Ca_buffer___Bmax_SLHigh_SL = 0.0134;   // millimolar
   Ca_buffer___Bmax_SLHigh_jct = 0.00165;   // millimolar
   Ca_buffer___koff_Calsequestrin = 65.0;   // per_millisecond
   Ca_buffer___koff_SLB = 1.3;   // per_millisecond
   Ca_buffer___koff_SLHigh = 30.0e-3;   // per_millisecond
   Ca_buffer___kon_Calsequestrin = 100.0;   // per_millimolar_per_millisecond
   Ca_buffer___kon_SL = 100.0;   // per_millimolar_per_millisecond
   Jleak_SR___KSRleak = 5.348e-6;   // per_millisecond
   Jpump_SR___H = 1.787;   // dimensionless
   Jpump_SR___Kmf = 0.000246;   // millimolar
   Jpump_SR___Kmr = 1.7;   // millimolar
   Jpump_SR___Q10_SRCaP = 2.6;   // dimensionless
   Jpump_SR___V_max = 286.0e-6;   // millimolar_per_millisecond
   Jrel_SR___EC50_SR = 0.45;   // millimolar
   Jrel_SR___HSR = 2.5;   // dimensionless
   Jrel_SR___Max_SR = 15.0;   // dimensionless
   Jrel_SR___Min_SR = 1.0;   // dimensionless
   Jrel_SR___kiCa = 0.5;   // per_millimolar_per_millisecond
   Jrel_SR___kim = 0.005;   // per_millisecond
   Jrel_SR___koCa = 10.0;   // per_millimolar2_per_millisecond
   Jrel_SR___kom = 0.06;   // per_millisecond
   Jrel_SR___ks = 25.0;   // per_millisecond
   cytosolic_Ca_buffer___Bmax_Calmodulin = 0.024;   // millimolar
   cytosolic_Ca_buffer___Bmax_Myosin_Ca = 0.14;   // millimolar
   cytosolic_Ca_buffer___Bmax_Myosin_Mg = 0.14;   // millimolar
   cytosolic_Ca_buffer___Bmax_SRB = 0.0171;   // millimolar
   cytosolic_Ca_buffer___Bmax_TroponinC = 0.07;   // millimolar
   cytosolic_Ca_buffer___Bmax_TroponinC_Ca_Mg_Ca = 0.14;   // millimolar
   cytosolic_Ca_buffer___Bmax_TroponinC_Ca_Mg_Mg = 0.14;   // millimolar
   cytosolic_Ca_buffer___koff_Calmodulin = 238.0e-3;   // per_millisecond
   cytosolic_Ca_buffer___koff_Myosin_Ca = 0.46e-3;   // per_millisecond
   cytosolic_Ca_buffer___koff_Myosin_Mg = 0.057e-3;   // per_millisecond
   cytosolic_Ca_buffer___koff_SRB = 60.0e-3;   // per_millisecond
   cytosolic_Ca_buffer___koff_TroponinC = 19.6e-3;   // per_millisecond
   cytosolic_Ca_buffer___koff_TroponinC_Ca_Mg_Ca = 0.032e-3;   // per_millisecond
   cytosolic_Ca_buffer___koff_TroponinC_Ca_Mg_Mg = 3.33e-3;   // per_millisecond
   cytosolic_Ca_buffer___kon_Calmodulin = 34.0;   // per_millimolar_per_millisecond
   cytosolic_Ca_buffer___kon_Myosin_Ca = 13.8;   // per_millimolar_per_millisecond
   cytosolic_Ca_buffer___kon_Myosin_Mg = 15.7e-3;   // per_millimolar_per_millisecond
   cytosolic_Ca_buffer___kon_SRB = 100.0;   // per_millimolar_per_millisecond
   cytosolic_Ca_buffer___kon_TroponinC = 32.7;   // per_millimolar_per_millisecond
   cytosolic_Ca_buffer___kon_TroponinC_Ca_Mg_Ca = 2.37;   // per_millimolar_per_millisecond
   cytosolic_Ca_buffer___kon_TroponinC_Ca_Mg_Mg = 3.0e-3;   // per_millimolar_per_millisecond
   model_parameters___F = 96486.7;   // coulomb_per_mole
   model_parameters___R = 8314.3;   // joule_per_kilomole_kelvin
   model_parameters___T = 308.0;   // kelvin
   model_parameters___cell_length = 100.0;   // micrometre
   model_parameters___cell_radius = 10.25;   // micrometre

   //---------------------------------------------------------------------------
   // Computed variables
   //---------------------------------------------------------------------------

   model_parameters___Vol_Cell = 0.5*3.141592654*pow(model_parameters___cell_radius/1000.0, 2.0)*model_parameters___cell_length/pow(1000.0, 3.0);  // 0.4
   model_parameters___Vol_cytosol = 0.65*model_parameters___Vol_Cell;
   model_parameters___Vol_SR = 0.035*model_parameters___Vol_Cell;
   model_parameters___Vol_SL = 0.02*model_parameters___Vol_Cell;
   model_parameters___Vol_jct = 0.00051*model_parameters___Vol_Cell;
/*
   printf ("%f\n", 1000000000*model_parameters___Vol_Cell);
   printf ("%f\n", 1000000000*model_parameters___Vol_cytosol);
   printf ("%f\n", 1000000000*model_parameters___Vol_SR);
   printf ("%f\n", 1000000000*model_parameters___Vol_SL);
   printf ("%f\n", 1000000000*model_parameters___Vol_jct);
*/
}

void comp_cai ()
{
   int i;

   Ca_buffer___dCalsequestrin = 1000*(Ca_buffer___kon_Calsequestrin*Y[6]*(Ca_buffer___Bmax_Calsequestrin*model_parameters___Vol_cytosol/model_parameters___Vol_SR-Y[0])-Ca_buffer___koff_Calsequestrin*Y[0]);
   dY[0] = Ca_buffer___dCalsequestrin;
   Ca_buffer___dCa_SLB_SL = 1000*(Ca_buffer___kon_SL*Y[1]*(Ca_buffer___Bmax_SLB_SL*model_parameters___Vol_cytosol/model_parameters___Vol_SL-Y[2])-Ca_buffer___koff_SLB*Y[2]);
   Ca_buffer___dCa_SLB_jct = 1000*(Ca_buffer___kon_SL*Y[7]*(Ca_buffer___Bmax_SLB_jct*0.1*model_parameters___Vol_cytosol/model_parameters___Vol_jct-Y[3])-Ca_buffer___koff_SLB*Y[3]);
   Ca_buffer___dCa_SLHigh_SL = 1000*(Ca_buffer___kon_SL*Y[1]*(Ca_buffer___Bmax_SLHigh_SL*model_parameters___Vol_cytosol/model_parameters___Vol_SL-Y[4])-Ca_buffer___koff_SLHigh*Y[4]);
   Ca_buffer___dCa_SLHigh_jct = 1000*(Ca_buffer___kon_SL*Y[7]*(Ca_buffer___Bmax_SLHigh_jct*0.1*model_parameters___Vol_cytosol/model_parameters___Vol_jct-Y[5])-Ca_buffer___koff_SLHigh*Y[5]);
   dY[2] = Ca_buffer___dCa_SLB_SL;
   dY[3] = Ca_buffer___dCa_SLB_jct;
   dY[4] = Ca_buffer___dCa_SLHigh_SL;
   dY[5] = Ca_buffer___dCa_SLHigh_jct;
   Ca_buffer___dCa_jct_tot_bound = Ca_buffer___dCa_SLB_jct+Ca_buffer___dCa_SLHigh_jct;
   Ca_buffer___dCa_SL_tot_bound = Ca_buffer___dCa_SLB_SL+Ca_buffer___dCa_SLHigh_SL;

   Ca_buffer___i_Ca_jct_tot = ICaL+ICaT;
   Ca_buffer___i_Ca_SL_tot = -2.0*INaCa+IbCa+ICap;
   
   Jpump_SR___j_pump_SR = 5*(1000*Jpump_SR___V_max*model_parameters___Vol_cytosol/model_parameters___Vol_SR*(pow(Y[8]/Jpump_SR___Kmf, Jpump_SR___H)-pow(Y[6]/Jpump_SR___Kmr, Jpump_SR___H))/(1.0+pow(Y[8]/Jpump_SR___Kmf, Jpump_SR___H)+pow(Y[6]/Jpump_SR___Kmr, Jpump_SR___H)));
   Jleak_SR___j_leak_SR = 5*1000*Jleak_SR___KSRleak*(Y[6]-Y[7]);
   Jrel_SR___j_rel_SR = 5*1000*Jrel_SR___ks*Y[10]*(Y[6]-Y[7]); 
   
   dY[6] = Jpump_SR___j_pump_SR-(Jleak_SR___j_leak_SR*model_parameters___Vol_cytosol/model_parameters___Vol_SR+Jrel_SR___j_rel_SR)-Ca_buffer___dCalsequestrin;
   ion_diffusion___J_Ca_jct_SL = 1000*(Y[7]-Y[1])*8.2413e-13;
   dY[7] = -0.000000001*Ca_buffer___i_Ca_jct_tot/(model_parameters___Vol_jct*2.0*model_parameters___F)-ion_diffusion___J_Ca_jct_SL/model_parameters___Vol_jct+Jrel_SR___j_rel_SR*model_parameters___Vol_SR/model_parameters___Vol_jct+Jleak_SR___j_leak_SR*model_parameters___Vol_cytosol/model_parameters___Vol_jct-1.0*Ca_buffer___dCa_jct_tot_bound;
   ion_diffusion___J_Ca_SL_cytosol = 1000*(Y[1]-Y[8])*3.7243e-12;
   dY[1] = -0.000000001*Ca_buffer___i_Ca_SL_tot/(model_parameters___Vol_SL*2.0*model_parameters___F)+(ion_diffusion___J_Ca_jct_SL-ion_diffusion___J_Ca_SL_cytosol)/model_parameters___Vol_SL-1.0*Ca_buffer___dCa_SL_tot_bound;
   cytosolic_Ca_buffer___dCa_TroponinC = 1000*(cytosolic_Ca_buffer___kon_TroponinC*Y[8]*(cytosolic_Ca_buffer___Bmax_TroponinC-Y[15])-cytosolic_Ca_buffer___koff_TroponinC*Y[15]);
   cytosolic_Ca_buffer___dCa_TroponinC_Ca_Mg = 1000*(cytosolic_Ca_buffer___kon_TroponinC_Ca_Mg_Ca*Y[8]*(cytosolic_Ca_buffer___Bmax_TroponinC_Ca_Mg_Ca-(Y[16]+Y[18]))-cytosolic_Ca_buffer___koff_TroponinC_Ca_Mg_Ca*Y[16]);
   cytosolic_Ca_buffer___dMg_TroponinC_Ca_Mg = 1000*(cytosolic_Ca_buffer___kon_TroponinC_Ca_Mg_Mg*(cytosolic_Ca_buffer___Bmax_TroponinC_Ca_Mg_Mg-(Y[16]+Y[18]))-cytosolic_Ca_buffer___koff_TroponinC_Ca_Mg_Mg*Y[18]);
   cytosolic_Ca_buffer___dCa_Calmodulin = 1000*(cytosolic_Ca_buffer___kon_Calmodulin*Y[8]*(cytosolic_Ca_buffer___Bmax_Calmodulin-Y[12])-cytosolic_Ca_buffer___koff_Calmodulin*Y[12]);
   cytosolic_Ca_buffer___dCa_Myosin = 1000*(cytosolic_Ca_buffer___kon_Myosin_Ca*Y[8]*(cytosolic_Ca_buffer___Bmax_Myosin_Ca-(Y[13]+Y[17]))-cytosolic_Ca_buffer___koff_Myosin_Ca*Y[13]);
//   cytosolic_Ca_buffer___dMg_Myosin = 1000*(cytosolic_Ca_buffer___kon_Myosin_Mg*(cytosolic_Ca_buffer___Bmax_Myosin_Mg-(Y[13]+Y[17]))-cytosolic_Ca_buffer___koff_Myosin_Mg*Y[17]);
   cytosolic_Ca_buffer___dCa_SRB = 1000*(cytosolic_Ca_buffer___kon_SRB*Y[8]*(cytosolic_Ca_buffer___Bmax_SRB-Y[14])-cytosolic_Ca_buffer___koff_SRB*Y[14]);
   cytosolic_Ca_buffer___dCa_cytosol_tot_bound = cytosolic_Ca_buffer___dCa_TroponinC+cytosolic_Ca_buffer___dCa_TroponinC_Ca_Mg+cytosolic_Ca_buffer___dMg_TroponinC_Ca_Mg+cytosolic_Ca_buffer___dCa_Calmodulin+cytosolic_Ca_buffer___dCa_Myosin+cytosolic_Ca_buffer___dMg_Myosin+cytosolic_Ca_buffer___dCa_SRB;
   dY[8] = -Jpump_SR___j_pump_SR*model_parameters___Vol_SR/model_parameters___Vol_cytosol+ion_diffusion___J_Ca_SL_cytosol/model_parameters___Vol_cytosol-1.0*cytosolic_Ca_buffer___dCa_cytosol_tot_bound;

   Jrel_SR___kCaSR = Jrel_SR___Max_SR-(Jrel_SR___Max_SR-Jrel_SR___Min_SR)/(1.0+pow(Jrel_SR___EC50_SR/Y[6], Jrel_SR___HSR));
   Jrel_SR___koSRCa = Jrel_SR___koCa/Jrel_SR___kCaSR;
   Jrel_SR___kiSRCa = Jrel_SR___kiCa*Jrel_SR___kCaSR;
   Jrel_SR___RI = 1.0-Y[11]-Y[10]-Y[9];
   dY[11] = 275*(Jrel_SR___kim*Jrel_SR___RI-Jrel_SR___kiSRCa*Y[7]*Y[11]-(Jrel_SR___koSRCa*pow(Y[7], 2.0)*Y[11]-Jrel_SR___kom*Y[10])); //1000 1
   
   //dY[11]= 275*Y[7]*(1-Y[11])-2.9*Y[11]; //Rakan's
   dY[10] = 1000*(Jrel_SR___koSRCa*pow(Y[7], 2.0)*Y[11]-Jrel_SR___kom*Y[10]-(Jrel_SR___kiSRCa*Y[7]*Y[10]-Jrel_SR___kim*Y[9]));
   dY[9] = 1000*(Jrel_SR___kiSRCa*Y[7]*Y[10]-Jrel_SR___kim*Y[9]-(Jrel_SR___kom*Y[9]-Jrel_SR___koSRCa*pow(Y[7], 2.0)*Jrel_SR___RI));

   dY[15] = cytosolic_Ca_buffer___dCa_TroponinC;
   dY[16] = cytosolic_Ca_buffer___dCa_TroponinC_Ca_Mg;
   dY[18] = cytosolic_Ca_buffer___dMg_TroponinC_Ca_Mg;
   dY[12] = cytosolic_Ca_buffer___dCa_Calmodulin;
   dY[13] = cytosolic_Ca_buffer___dCa_Myosin;
//   dY[17] = cytosolic_Ca_buffer___dMg_Myosin;
   dY[14] = cytosolic_Ca_buffer___dCa_SRB;

   for (i = 0; i <= 18; i++)
	 Y[i] += HT*dY[i];	  
}

void comp_iclca()
{
     
double jss_dist;
double Kd;
double Ass1, Bss1;
double Ass2, Bss2;
double CICF, CICF6;

jss_dist = 0.015; //distribution in junctional subsarcommel space
  

Kd=Kb/Kf;

Ass1 = Kf*pow(Y[1], 1.7);
Bss1 = Kb*pow(Kd, 1.7);

ss1 = Ass1/(Ass1+Bss1);
ss1 += HT*(Ass1*(1-ss1)-Bss1*ss1);

CICF = (V/RTONF)*((Y[1]- 2.0)*exp(V/RTONF))/(1-exp(V/RTONF));

IClCa_ss = CICF * pow(ss1, 1.7);

Ass2 = Kf*pow(Y[7], 1.7);
Bss2 = Kb*pow(Kd, 1.7);
ss2 = Ass2/(Ass2+Bss2);
ss2 += HT*(Ass2*(1-ss2)-Bss2*ss2);

CICF6 = (V/RTONF)*((Y[7]- 2.0)*exp(V/RTONF))/(1-exp(V/RTONF));

IClCa_jss = CICF6 * ss2;


IClCa = 435.93 * (IClCa_ss*(1-jss_dist) + IClCa_jss*(jss_dist));

}

void comp_ikh()
{
   double Ay, By;
   double ym, ty;
   
  Ay = exp(-(2.9 + (0.04*V)));
  By  = exp((3.6 + (0.11*V)));

   ym = 1/(1 + exp((V+100)/20));
      
  ty = 0.4/(Ay+By);

  y2 += HT*(ym-y2)/ty;

  IKH_K= 0.1*ym*(V-Ek);
  IKH_Na = 0.008*ym*(V-ENa);
  IKH = (IKH_K + IKH_Na);
   
}

/*void comp_if()
{
   double Ay, By;
   double ym, ty;
   
   Ay = exp(-(2.9 + (0.04*V)));
   By  = exp((3.6 + (0.11*V)));

   ym = 1/(1 + exp((V+100)/10));
   ty = 0.4/(Ay+By);
   //ty=0.4;
 
   y2 += HT*(ym-y2)/ty;

   IfNa = 0.00072827*ym*(V-ENa); //0.72827
   IfK = 0.00117173*ym*(V-Ek); //1.17173
   If = 4400*(IfNa + IfK); //10*
}*/

void comp_ikf ()
{
   double Apa, Bpa;
   double pam, tpa;
   double Api, Bpi;
   double pim, tpi;

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
   double An, Bn;
   double nm, tn;

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
   double E0=V-Ek+3.6;

//   Ik1=0.52*5.08*pow(Kc/(Kc+0.59),3)*(V-Ek)/(1+exp(1.393*E0/RTONF));     // CT - 2.0, PM - 2.5
     Ik1=(0.9*5.08*pow(5.0/(5.0+0.59),3)/(1+exp(1.393*(V-Ek-12)/RTONF))+0.05)*(V-Ek-12);

}

void comp_ito ()
{
   double Ar, Br;
   double rm, tr;
   double s1m, ts1;
   double s2m, ts2;
   double s3m, ts3;
   double xsusm;

   Ito=0.45*0.515*50.02*r1*(0.8*pow(s1,3)+0.2*pow(s2,3))*/*(0.600*pow(s3,6)+0.4)*/(V-Ek);     // CT - 0.2, PM - 0.35
//   Isus = 0.6*(V+70);  // 1.4, 2.4
  Isus = 0.46*(V+85)+0.0*(V-5);
   
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
//INaCa=0.02*(pow(Nai,3)*2.5*exp(0.450*V/RTONF)-pow(140,3)*Y[8]*exp(V*(0.45-1)/RTONF))/(1+0.0003*(Y[8]*pow(140,3)+2.5*pow(Nai,3)));  
   
//SAP
INaCa = 1.5*(0.02*(pow(Nai,3)*7.5*exp(0.35*V/RTONF)-pow(140,3)*Y[8]*exp(V*(0.45-1)/RTONF))/(1+0.0003*(Y[8]*pow(140,3)+2.5*pow(Nai,3)))); //1.5

}

void comp_ina ()
{
   double Am, Bm;
   double Ah, Bh;
   double hm;
   double th1, th2;

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
   ICap=9.509*(Y[8]/(Y[8]+0.0002));   
}

void comp_ical ()
{
   double Adl, Bdl;
   double dlm, tdl;
   double Afl, Bfl;
   double flm, tfl1, tfl2;
   double dfca;
   
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

//ultraslow
  /*flm=1./(1+exp((V+30.153)/5.122));
  tfl1=(0.0336/(1+exp((V+18)/4.0)))+0.0142;
  tfl2=(3.6908/(1+exp((V+18)/4)))+0.0534;
  fL1+=HT*(flm-fL1)/tfl1;
  fL2+=HT*(flm-fL2)/tfl2;*/
  
  //Boyden et al
  flm=1./(1+exp((V+23.4)/5.2));
  tfl1=(0.0336/(1+exp((V+18)/4.0)))+0.0142;
  tfl2=(3.6908/(1+exp((V+18)/4)))+0.0534;
  fL1+=HT*(flm-fL1)/tfl1;
  fL2+=HT*(flm-fL2)/tfl2;
  
  
  dfca=1.0/(1+Y[7]/0.00025);  //0.000119
  fca+=HT*(dfca-fca)/0.02;

//ultraslow
//   ICaL=30*(dL*((0.8*fL1)+(0.2*fL2)))*fca*(V-58.3);   
     ICaL=25*dL*fL1*fca*(V-58.3);
 
}

void comp_icat ()
{
   double Adt, Bdt;
   double dtm, tdt;
   double Aft, Bft;
   double ftm, tft;

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
  IbNa=0.02*(V-ENa);  // 0.02
}

void comp_ibca ()
{
  IbCa=0.02*(V-ECa);  // 0.02
}

void comp_itot ()
{
   Itot=IKf+Ik1+Ito+Ip+INaCa+INa+IbNa+ICap+ICaL+ICaT+IbCa+If+IKH+IClCa;
}

/*void comp_nai ()
{
   Nai+=HT*(-3*Ip-3*INaCa-IbNa-INa)/(Faraday*Vi);
}*/

main ()
{
      time_t start, end;
     time(&start);
     printf("%s", ctime ( &start ));
       
   double t=0.0;
   double Istim;
   double dvdt;
   int cnt=0;
   int cnt2=0;
   int cnt3=0;
   int i=0;

   double tstim=0.01;
   double stimtime;
   char str[80];
   char str2[10];

   double vmax[beats+1];
   double dvdtmax[beats+1];
   double apd[beats+1];
   double toneapd[beats+1];
   double ttwoapd[beats+1];
   double rmbp[beats+1];

   FILE *in;
   
   Vcell=0.003497;   
   Vi=0.465*Vcell;
   Vc=0.1*Vcell;
   
   RTONF=Temp*0.08554;
   
   init();
   
   in = fopen ("Shannon_PV.txt","w"); 
             
    TMAX=BCL*beats+tstim;
   for (i=beats; i>0; i--) {
     vmax[i] = -100;
	 dvdtmax[i] = 0;
   }

   while (t < TMAX) {
	 t+=HT;
     
	 Ek=RTONF*log(Kc/Ki);
     ENa=RTONF*log(140./Nai);     // 140.
     ECa=(RTONF/2.)*log(2.5/Y[8]);  
     ECl=RTONF*exp(2/Y[2]); 
     EbCl=ECl-0.49*(ECl+30.59);

     comp_ikh () ;
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
	 comp_iclca ();
	 comp_cai ();
   //  comp_nai ();
   


     if (t>=tstim && t<(tstim+HT)) {
       stimtime = 0;
       i = i+1;

       tstim = tstim+BCL;
       printf ("\nStimulus %d applied, ", i);
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
     if (cnt == 100) { 
       //fprintf (in, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", t, V, 1000*Y[8], ICaL/Cm, ICaT/Cm, Ik1/Cm, IKf/Cm, Ito/Cm, IClCa/Cm, INaCa/Cm, Ip/Cm, IbNa/Cm, IbCa/Cm);
	   fprintf (in, "%lf\t%lf\t%lf\n", t, V, Y[8]*1000);
       cnt = 0; 
       
       if(cnt2==200) {printf("."); cnt2=0;}
        if(cnt3 == 150000) {      printf (" CALCULATING ");
       printf ("t = %f ", t); cnt3 = 0;}
	 }
   } 

   fclose (in);

   time(&end);
   double diff;
   diff = difftime(end, start);
   apd[beats]=ttwoapd[beats]-toneapd[beats];
   printf ("APD = %8.2lf, dVdtmax = %8.2lf, Runtime = %0.2lf\n", 1000*apd[beats], 0.001*dvdtmax[beats], diff);
  
   system("pause");
}


