#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "LA_cell.h"
#include "RA_cell.h"
#include "PM_cell.h"
#include "BB_cell.h"
#include "CS_cell.h"
#include "PV_cell.h"
#include "CT_cell.h"
#include "Atria_Variables.h"

float HX=0.1;  
float HT=0.00001;
 
float fac=0.029108;
float faTc=0.014071;
float faTmgc=0.214036;
float faTmgm=0.693565;
float faCalse=0.465921;


const float Temp=308;
const float Cm=0.05;
const float Faraday=96487;
const float R=8314;
float RTONF;


#define BLOCK_LOW(id,p,n) ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id+1),p,n)-1)
#define BLOCK_SIZE(id,p,n) (BLOCK_LOW((id+1),p,n)-BLOCK_LOW(id,p,n))

int main(int argc, char **argv)
{

  int oneangularstep = 360/lucy; //we need this to be an integer for the steps, however, it clearly isn't 360/314

  int pv_radius = 50; 

  int xmap, ymap;
  int xmap1, ymap1;


   float global_V[X+1][Y+lucy+1]; // added +l to x component to extend sheet as above, didn't need to for y axis as no extension in this direction
   float new_V[X+1][Y+lucy+1];
   float dvdt[X+1][Y+lucy+1];
   int geometry[X+1][Y+lucy+1];

   float dvdx2, dvdy2;

   float t=0.0;

   int x, y;
   int cnt=0;
   int num=0;
   
   int cnt2 = 0;

   FILE *out;
   char *str;
   
   RTONF=Temp*0.08554;

   ECl=RTONF*log(30./132.); 
   EbCl=ECl-0.49*(ECl+30.59);
   
   Vi=0.0126;
   VCa=0.005884;
   Vc=0.0025;
   Vup=0.0003969;
   Vrel=0.000044;
   

   int rank = 0;
   int size = 1;

   int *recv_cnts, *recv_disp ;
   MPI_Status status ;

   MPI_Init(&argc, &argv) ;
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   recv_cnts = calloc (size, sizeof(int)) ;
   recv_disp = calloc (size, sizeof(int)) ;

   for (int i = 0 ; i < size ; i++) {
     recv_cnts[i] = BLOCK_SIZE(i, size, X+1)*(Y+lucy+1);
     recv_disp[i] = BLOCK_LOW(i, size, X+1)*(Y+lucy+1);
   } 


// ** Initialise all variables at all co-ordinated ** //
   for (x=0; x<=X; x++) {
       
     for (y=0; y<=Y+lucy; y++) {
      
       V[x][y]=-80.0; 
	   new_V[x][y]=-73.8; //?
       Pa[x][y]=0.00016;
       Pi[x][y]=0.76898;
       n[x][y]=0.02032;
       r1[x][y]=0.00006;
       s1[x][y]=0.5753;
       s2[x][y]=0.39871;
       s3[x][y]=0.57363;
       m[x][y]=0.01309;
       h1[x][y]=0.706;
       h2[x][y]=0.61493;
       dL[x][y]=0.00003;
       fL1[x][y]=0.99981;
       fL2[x][y]=0.99981;
       y2[x][y]=0.0927;
       fca[x][y]=0.7755; //couldn't find in 1
       dT[x][y]=0.00046;
       fT[x][y]=0.30752;
       F1[x][y]=0.288039; // F1,2,3 not in 1
       F2[x][y]=0.002262;
       F3[x][y]=0.612697;
       Cai[x][y]=0.000071;
       Caup[x][y]=0.730866;
       Carel[x][y]=0.726776;
       Kc[x][y]=5.0;
       Ki[x][y]=140.0;
       Nai[x][y]=8.4; //not in 1
      }
      }
            
// ** End Initialisation of variables ** //

// ** Set the stimulus sites ** //

   for (x=0; x<=5; x++)
	 for (y=0; y<=Y; y++)
       V[x][y] = 20.0;

// ** Begin voltage calculation ** //

  for (x=0;x<=X;x++)
  {
       for (y=0; y<=Y+lucy; y++) {
                
//LA Background 
      if (y<300)
         geometry[x][y]=6;
//RA Background      
      else
         geometry[x][y]=7;
               
//PM Coordinates
		 if ((y>350 && y<= 450) && ((x >= 40 && x < 60) || (x >= 80 && x <100) || (x >= 120 && x < 140) ||
                                    (x >= 160 && x < 180) || (x >= 200 && x <220 ) || (x >= 240 && x < 260))){
                                    geometry[x][y]=1; 
                                    }
//CT Coordinates
         if((y>=325 && y<= 350) && (x>=25 && x<=275)){
                    geometry[x][y]=2;
                    }
//BB Coordinates
         if((y>=250 && y<= 350) && (x>=0 && x<25)){
                    geometry[x][y]=3;}
//CS Coordinates
         if((y>=250 && y<= 350) && (x>=275 && x<X)){
                    geometry[x][y]=4;}

//PV Coordinates
       if(y>=Y) 
         geometry[x][y]=5;
	  }
  }
	 
   while (t < 100.0) { 
              
	 t+=HT;
     cnt++;
     
     for (x = BLOCK_LOW(rank, size, X+1); x <= BLOCK_HIGH(rank, size, X+1); x++) {
       if (x == 0 || x == X) continue;

       for (y=1; y<Y+lucy; y++) 
	   
	   if (((x-100)*(x-100)+(y-100)*(y-100))>pv_radius*pv_radius) // only do this if you're outside the hole where PV connects

	   {              
   
		   
	     Ek=RTONF*log(Kc[x][y]/Ki[x][y]);
         ENa=RTONF*log(140./Nai[x][y]);     
         ECa=(RTONF/2.)*log(2.5/Cai[x][y]);

// ** Define the differential approximation ** //
	     dvdx2 = (V[x+1][y]+V[x-1][y]-2*V[x][y])/(HX*HX);
	     dvdy2 = (V[x][y+1]+V[x][y-1]-2*V[x][y])/(HX*HX);
	    
	 switch (geometry[x][y]) {
       case 1: dvdt[x][y]=D*(dvdx2+dvdy2)-I_tot_PM(x, y); break;
       case 2: dvdt[x][y]=D*(dvdx2+dvdy2)-I_tot_CT(x, y); break;
       case 3: dvdt[x][y]=D*(dvdx2+dvdy2)-I_tot_BB(x, y); break; 
       case 4: dvdt[x][y]=D*(dvdx2+dvdy2)-I_tot_CS(x, y); break;
       case 5: dvdt[x][y]=0.5*D*(dvdx2+dvdy2)-I_tot_PV(x, y); break; //D is the diffusion coefficient, notice smaller for PV.  now halved to 0.25
       case 6: dvdt[x][y]=D*(dvdx2+dvdy2)-I_tot_LA(x, y); break;
       case 7: dvdt[x][y]=D*(dvdx2+dvdy2)-I_tot_RA(x, y); break;
//     default: dvdt[x][y]=D*(dvdx2+dvdy2)-I_tot_RA(x, y); 
	}

       new_V[x][y]=V[x][y]+HT*dvdt[x][y];
	   }
   }

// ** Set the boundary conditions on the tissue ** //
	   
for (x = BLOCK_LOW(rank, size, X+1); x <= BLOCK_HIGH(rank, size, X+1); x++) {

// ** Set values at edge cells to the value of the cells at (edge -1) ** //

	   new_V[x][0]=new_V[x][1];
       new_V[x][Y]=new_V[x][Y-1];

// ** Zero-flux conditions at the atrial septum ** //
	   if (x >= 25 && x <= 275) {
	     new_V[x][300]=new_V[x][299];
	     new_V[x][301]=new_V[x][302];
	   }
}

// ** Set values at edge cells to the value of the cells at (edge -1) ** //

for (y=0; y<=Y; y++)
{   
	   new_V[0][y]=new_V[1][y];
       new_V[X][y]=new_V[X-1][y];
}


for (y=Y+1; y<=Y+lucy; y++)
{    //for the bottom of the sheet to match the pulmonary vein area of the main sheet
  new_V[X][y]=new_V[X-1][y];

  ymap = 100+pv_radius*cos(2*3.1415926*(y-Y)/lucy);
  xmap = 100+pv_radius*sin(2*3.1415926*(y-Y)/lucy);

  ymap1 = 100+(pv_radius+1)*cos(2*3.1415926*(y-Y)/lucy);
  xmap1 = 100+(pv_radius+1)*sin(2*3.1415926*(y-Y)/lucy);

  V[0][y] = V[xmap][ymap];
  V[xmap1][ymap1] = V[1][y];
}

for (x=0; x<=X; x++) 
// boundary condition on the x axis forcing a cylinder shape where both ends are of equal V
{
  new_V[x][Y+1] = new_V[x][Y+lucy-1];
  new_V[x][Y+lucy] = new_V[x][Y+2];
}


// ** Replace the old Voltage values with the new ones ** //
  for (x = BLOCK_LOW(rank, size, X+1); x <= BLOCK_HIGH(rank, size, X+1); x++) 
    for (y=0; y<=Y; y++)
	if (((x-100)*(x-100)+(y-100)*(y-100))>pv_radius*pv_radius)
    { 
		V[x][y]=new_V[x][y];
    }
	 
   if (t>=0.04 && t<(0.04+HT))
	 for (x = BLOCK_LOW(rank, size, X+1); x <= BLOCK_HIGH(rank, size, X+1); x++)
	 {
         for (y=0; y<=200; y++)
         {
		   V[x][y] = -73.8;
		   new_V[x][y] = -73.8;
		   
		 }
		 for (y=400; y<=Y; y++)
         {
		   V[x][y] = -73.8;
		   new_V[x][y] = -73.8;
		 }
     }
     
	 
	 if (rank < (size-1))
       MPI_Send(V[BLOCK_HIGH(rank, size, X+1)], Y+1, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD);
     if (rank > 0)
       MPI_Recv(V[BLOCK_HIGH((rank-1), size, X+1)], Y+1, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, &status);
     if (rank > 0)
       MPI_Send(V[BLOCK_LOW(rank, size, X+1)], Y+1, MPI_FLOAT, rank-1, 1, MPI_COMM_WORLD);
     if (rank < (size-1))
       MPI_Recv(V[BLOCK_LOW((rank+1), size, X+1)], Y+1, MPI_FLOAT, rank+1, 1, MPI_COMM_WORLD, &status);

     if (cnt % 1000 == 0) {     //took off a zero
       MPI_Gatherv (&(((float *)V)[recv_disp[rank]]), recv_cnts[rank], MPI_FLOAT, &global_V[0][0], recv_cnts, recv_disp, MPI_FLOAT, 0, MPI_COMM_WORLD);
     }

     if (!rank)
 
		 if (cnt == 1000) //ditto 
         { 
//         cnt2++;
//         printf("Step %4d\n", cnt2);
	     cnt=0;         
         
	     str = malloc (8*sizeof(char));
	     sprintf (str, "new%d.vtk", num++);
         out = fopen (str, "wt");
         
         fprintf (out, "# vtk DataFile Version 3.0\n");   
         fprintf (out, "vtk output\n");
         fprintf (out, "ASCII\n");
         fprintf (out, "DATASET STRUCTURED_POINTS\n");
         fprintf (out, "DIMENSIONS %d %d %d\n", X+1, Y+lucy+1, 1);
         fprintf (out, "SPACING 1 1 1\n");
         fprintf (out, "ORIGIN 0 0 0\n");
         fprintf (out, "POINT_DATA %d\n", (X+1)*(Y+lucy+1));
         fprintf (out, "SCALARS ImageFile float 1\n");
         fprintf (out, "LOOKUP_TABLE default\n"); 

	     for (y=0; y<=Y+lucy; y++)
         {
           for (x=0; x<=X; x++)
              fprintf (out, "%2.2lf ", global_V[x][y]);
              
	       fprintf (out, "\n");
	       
         }  
	     fclose (out);
	     free (str);
         } 
}
return 0;

}