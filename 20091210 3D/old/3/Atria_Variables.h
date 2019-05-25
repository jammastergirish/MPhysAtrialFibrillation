#include <stdlib.h>

#define X   300
#define Y   600
#define Z   50
#define D   100


float HT, HX;

float Ek,ENa,ECa,EbCl,ECl;
float Vcell,Vi,VCa,Vc,Vup,Vrel;
float rinact;
float dfac,dfaTc,dfaTmgc,dfaTmgm,dfaCalse,dfab;

float V[X+1][Y+1][Z+1];
float Pa[X+1][Y+1][Z+1];
float Pi[X+1][Y+1][Z+1];
float n[X+1][Y+1][Z+1];
float r1[X+1][Y+1][Z+1];
float s1[X+1][Y+1][Z+1];
float s2[X+1][Y+1][Z+1];
float s3[X+1][Y+1][Z+1];
float m[X+1][Y+1][Z+1];
float h1[X+1][Y+1][Z+1];
float h2[X+1][Y+1][Z+1];
float dL[X+1][Y+1][Z+1];
float fL1[X+1][Y+1][Z+1];
float fL2[X+1][Y+1][Z+1];
float y2[X+1][Y+1][Z+1];
float fca[X+1][Y+1][Z+1];
float dT[X+1][Y+1][Z+1];
float fT[X+1][Y+1][Z+1];
float F1[X+1][Y+1][Z+1];
float F2[X+1][Y+1][Z+1];
float F3[X+1][Y+1][Z+1];
float dF1[X+1][Y+1][Z+1];
float dF2[X+1][Y+1][Z+1];
float dF3[X+1][Y+1][Z+1];
float Kc[X+1][Y+1][Z+1];
float Ki[X+1][Y+1][Z+1];
float Nai[X+1][Y+1][Z+1];

float Cai[X+1][Y+1][Z+1];
float Caup[X+1][Y+1][Z+1];
float Carel[X+1][Y+1][Z+1];

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
   
