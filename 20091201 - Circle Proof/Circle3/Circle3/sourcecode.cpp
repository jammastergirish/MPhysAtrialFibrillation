#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define X   300
#define Y   600
#define D   100
#define lucy 180 // 2*pi()*radius of PV (r=50)

int main()
{

double y, x, xmap, ymap, pv_radius;

pv_radius = 50;
FILE *in;

   in = fopen ("AP.txt", "w");


for (y=Y+1; y<=Y+lucy; y++)
{    //for the bottom of the sheet to match the pulmonary vein area of the main sheet
 // new_V[X][y]=new_V[X-1][y];

  ymap = 100+pv_radius*cos(2*3.1415926*(y-Y)/lucy);
  xmap = 100+pv_radius*sin(2*3.1415926*(y-Y)/lucy);

  //cout << y << " " << xmap << " " << ymap << "\n";
  fprintf (in, "%lf %lf\n", xmap, ymap);

 // ymap1 = 100+(pv_radius+1)*cos(2*3.1415926*(y-Y)/lucy);
 // xmap1 = 100+(pv_radius+1)*sin(2*3.1415926*(y-Y)/lucy);

//  new_V[0][y] = new_V[xmap][ymap]; //changed from V
  //new_V[xmap1][ymap1] = new_V[1][y];
}
//cin >> y;
fclose(in);
return 0;
}
