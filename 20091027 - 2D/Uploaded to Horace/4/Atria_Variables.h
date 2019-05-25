#include <stdlib.h>

#define X   300
#define Y   600
#define D   100


float HT, HX;

float Ek,ENa,ECa,EbCl,ECl;
float Vcell,Vi,VCa,Vc,Vup,Vrel;
float rinact;
float dfac,dfaTc,dfaTmgc,dfaTmgm,dfaCalse,dfab;

//const float pv_radius = 50.0; //moved up just in case

//const float 314 = 2*3*pv_radius; //put in to define the extra 314ength on x axis for the 'new' sheet, a continuation of the o314d
//const int l = 314;
float V[X+314+1][Y+1];
float Pa[X+314+1][Y+1];
float Pi[X+314+1][Y+1];
float n[X+314+1][Y+1];
float r1[X+314+1][Y+1];
float s1[X+314+1][Y+1];
float s2[X+314+1][Y+1];
float s3[X+314+1][Y+1];
float m[X+314+1][Y+1];
float h1[X+314+1][Y+1];
float h2[X+314+1][Y+1];
float dL[X+314+1][Y+1];
float fL1[X+314+1][Y+1];
float fL2[X+314+1][Y+1];
float y2[X+314+1][Y+1];
float fca[X+314+1][Y+1];
float dT[X+314+1][Y+1];
float fT[X+314+1][Y+1];
float F1[X+314+1][Y+1];
float F2[X+314+1][Y+1];
float F3[X+314+1][Y+1];
float dF1[X+314+1][Y+1];
float dF2[X+314+1][Y+1];
float dF3[X+314+1][Y+1];
float Kc[X+314+1][Y+1];
float Ki[X+314+1][Y+1];
float Nai[X+314+1][Y+1];

float Cai[X+314+1][Y+1];
float Caup[X+314+1][Y+1];
float Carel[X+314+1][Y+1];

float ract;
float fac;
float faTc;
float faTmgc;
float faTmgm;
float faCalse;

const float Temp;
const float Cm;
const float Faraday;
const float R;
float RTONF;
   
