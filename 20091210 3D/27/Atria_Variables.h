#include <stdlib.h>

#define X   300
#define Y   600
#define Z   50
#define D   100

#define HX  0.1
#define HT  0.00001

#define Temp    308
#define Cm      0.05
#define Faraday 96487
#define R       8314

float Ek,ENa,ECa,EbCl,ECl;
float Vcell,Vi,VCa,Vc,Vup,Vrel;

int g[X+1][Y+1][Z+1];
int h[X+1][Y+1][Z+1];

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
float fac[X+1][Y+1][Z+1];
float faTc[X+1][Y+1][Z+1];
float faTmgc[X+1][Y+1][Z+1];
float faTmgm[X+1][Y+1][Z+1];
float faCalse[X+1][Y+1][Z+1];

float RTONF;
   
