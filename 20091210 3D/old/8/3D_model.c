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

#define BLOCK_LOW(id,p,n) ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id+1),p,n)-1)
#define BLOCK_SIZE(id,p,n) (BLOCK_LOW((id+1),p,n)-BLOCK_LOW(id,p,n))


int main(int argc, char **argv)
{   
   float new_V[X+1][Y+1][Z+1];
   float global_V[X+1][Y+1][Z+1]; // added this line 17th dec as variable wasn't defined
   float dvdt[X+1][Y+1][Z+1];

   float dvdx2, dvdy2, dvdz2;
   float t=0.0;

   int x, y, z;
   int cnt = 0;
   int cnt2 = 0;
   int dd;
  
   int isbound;
   int i, j, k;
   int num = 0;
   int gg[27];

   float root2 = sqrt(2.0);
   float root3 = sqrt(3.0);
   float ic, ir, il, imax;
   float tflt;

   FILE *in, *out;
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
     recv_cnts[i] = BLOCK_SIZE(i, size, X+1)*(Y+1)*(Z+1);
     recv_disp[i] = BLOCK_LOW(i, size, X+1)*(Y+1)*(Z+1);
   } 

  if (!rank) {

    in = fopen ("geometry.c", "r");	
	for (z = 0; z <= Z; z++) { 
      for (y = 0; y <= Y; y++) { 
		for (x = 0; x <= X; x++) {
		  fscanf(in, "%d ", &dd);
          h[x][y][z] = dd;     
          if (dd > 0)
			g[x][y][z] = 1;
		  else
            g[x][y][z] = 0;
		}
        fscanf(in, "\n"); 
      } 
      fscanf(in, "\n"); 
    }
    fclose (in);

    num = 1;      
    for (x = 1; x < X; x++) 
      for (y = 1; y < Y; y++) 
        for (z = 1; z < Z; z++)
          if (g[x][y][z] > 0) {
            g[x][y][z] = num;
            num++;
		  }   

	for (x = 1; x < X; x++) 
     for (y = 1; y < Y; y++) 
      for (z = 1; z < Z; z++) {
        gg[1] = g[x - 1][y - 1][z - 1];
        gg[2] = g[x - 1][y - 1][z];
        gg[3] = g[x - 1][y - 1][z+1];
        gg[4] = g[x - 1][y][z - 1];
        gg[5] = g[x - 1][y][z];
        gg[6] = g[x - 1][y][z + 1];
        gg[7] = g[x - 1][y + 1][z - 1];
        gg[8] = g[x - 1][y + 1][z];
        gg[9] = g[x - 1][y + 1][z + 1];

        gg[10] = g[x][y - 1][z - 1];
        gg[11] = g[x][y - 1][z];
        gg[12] = g[x][y - 1][z + 1];
        gg[13] = g[x][y][z - 1];
        gg[14] = g[x][y][z + 1];
        gg[15] = g[x][y + 1][z - 1];
        gg[16] = g[x][y + 1][z];
        gg[17] = g[x][y + 1][z + 1];

        gg[18] = g[x + 1][y - 1][z - 1];
        gg[19] = g[x + 1][y - 1][z];
        gg[20] = g[x + 1][y - 1][z + 1];
        gg[21] = g[x + 1][y][z - 1];
        gg[22] = g[x + 1][y][z];
        gg[23] = g[x + 1][y][z + 1];
        gg[24] = g[x + 1][y + 1][z - 1];
        gg[25] = g[x + 1][y + 1][z];
        gg[26] = g[x + 1][y + 1][z + 1];

        isbound = 0;
        for(i = 1; i <= 26; i++) { 
          if (gg[i] > 0) {gg[i] = 1; isbound++;} 
          else gg[i] = 0;
        }

        if (g[x][y][z] == 0 && isbound > 0) {
          ic = (gg[3]/root3) - (gg[1]/root3) + (gg[6]/root2) +
               (gg[9]/root3) - (gg[7]/root3) - (gg[4]/root2) +
               (gg[12]/root2) - (gg[10]/root2) + gg[14] +
               (gg[17]/root2) - (gg[15]/root2) - gg[13] +
               (gg[20]/root3) - (gg[18]/root3) + (gg[23]/root2) +
               (gg[26]/root3) - (gg[24]/root3) - (gg[21]/root2);

          ir = (gg[9]/root3) - (gg[2]/root2) - (gg[3]/root3) - 
               (gg[1]/root3) + (gg[8]/root2) + (gg[7]/root3) +
               (gg[17]/root2) - gg[11] - (gg[12]/root2) -
               (gg[10]/root2) + gg[16] + (gg[15]/root2) +
               (gg[26]/root3) - (gg[19]/root2) - (gg[20]/root3) -
               (gg[18]/root3) + (gg[25]/root2) + (gg[24]/root3);

          il = (gg[18]/root3) + (gg[19]/root2) + (gg[20]/root3) +
               (gg[21]/root2) + gg[22] + (gg[23]/root2) +
               (gg[24]/root3) + (gg[25]/root2) + (gg[26]/root3) -
               (gg[1]/root3) - (gg[2]/root2) - (gg[3]/root3) -
               (gg[4]/root2) - gg[5] - (gg[6]/root2) - 
               (gg[7]/root3) - (gg[8]/root2) - (gg[9]/root3);

          imax = fabs(ic);
          if (fabs(ir) > imax) imax = fabs(ir);
          if (fabs(il) > imax) imax = fabs(il);

          i = 0; j = 0; k = 0;

          tflt = ir / fabs(imax);
          if (tflt <= 0.5 && tflt >= -0.5) i = 0;
          else if (tflt > 0.5) i = 1;
          else if (tflt < -0.5) i = -1;

          tflt = ic / fabs(imax);
          if (tflt <= 0.5 && tflt >= -0.5) j = 0;
          else if (tflt > 0.5) j = 1;
          else if (tflt < -0.5) j = -1;

          tflt = il / fabs(imax);
          if (tflt <= 0.5 && tflt >= -0.5) k = 0;
          else if (tflt > 0.5) k = 1;
          else if (tflt < -0.5) k = -1;

          if (imax == 0) { i = 0; j = 0; k = 0; }

          if (g[x + k][y + i][z + j] > 0)    
            g[x][y][z] = -1 * g[x + k][y + i][z + j];   
          else
            g[x][y][z] = g[x + k][y + i][z + j];  
         }
      }
  }

  MPI_Bcast (g, (X+1)*(Y+1)*(Z+1), MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (h, (X+1)*(Y+1)*(Z+1), MPI_INT, 0, MPI_COMM_WORLD);

  // ** Initialise all variables at all co-ordinates ** //
   for (x=0; x<=X; x++)       
     for (y=0; y<=Y; y++) 
	   for (z=0; z<=Z; z++) {    
       V[x][y][z]=-100; // changed to -100 so that we can isolate it in paraview
	   new_V[x][y][z]=-73.8;
       Pa[x][y][z]=0.000028;
       Pi[x][y][z]=0.933953;
       n[x][y][z]=0.00925;
       r1[x][y][z]=0.0016;
       s1[x][y][z]=0.5753;
       s2[x][y][z]=0.5300;
       s3[x][y][z]=0.5844;
       m[x][y][z]=0.007788;
       h1[x][y][z]=0.87277;
       h2[x][y][z]=0.84549;
       dL[x][y][z]=0.00016;
       fL1[x][y][z]=0.9981;
       fL2[x][y][z]=0.9981;
       y2[x][y][z]=0.0927;
       fca[x][y][z]=0.7755;
       dT[x][y][z]=0.00046;
       fT[x][y][z]=0.30752;
       F1[x][y][z]=0.288039;
       F2[x][y][z]=0.002262;
       F3[x][y][z]=0.612697;
       Cai[x][y][z]=0.000073;
       Caup[x][y][z]=0.730866;
       Carel[x][y][z]=0.726776;
       Kc[x][y][z]=5.0;
       Ki[x][y][z]=140.0;
       Nai[x][y][z]=8.4;
	   fac[x][y][z]=0.029108;
       faTc[x][y][z]=0.014071;
       faTmgc[x][y][z]=0.214036;
       faTmgm[x][y][z]=0.693565;
       faCalse[x][y][z]=0.465921;
      }
            
// ** End Initialisation of variables ** //

// ** Set the stimulus sites ** //

   for (x=0; x<=10; x++)
	 for (y=0; y<=Y; y++)
	   for (z=0; z<=10; z++)
		 if (h[x][y][z] > 0) {              
           V[x][y][z] = 20.0; 
		   new_V[x][y][z] = 20.0;
		 }

		 
   for (x=0; x<=X; x++)
	 for (y=0; y<=Y; y++)
	   for (z=40; z<=Z; z++)
		 if (h[x][y][z] > 0) {              
           V[x][y][z] = 20.0; 
		   new_V[x][y][z] = 20.0;
		 }

// ** Begin voltage calculation ** //
	 
   while (t < 50.0) { 
              
	 t+=HT;
     cnt++;

	 for (x = BLOCK_LOW(rank, size, X+1); x <= BLOCK_HIGH(rank, size, X+1); x++) {
      if (x == 0 || x == X) continue;

      for (y = 1; y < Y; y++) 
		for (z = 1; z < Z; z++) 
		  if (g[x][y][z] < 0)
            for (i = -1; i <= 1; i++)
              for (j = -1; j <= 1; j++)
                for (k = -1; k <= 1; k++)
				  if (g[x][y][z] == -g[x + i][y + j][z + k]) 
                    V[x][y][z] = V[x + i][y + j][z + k];
	 }
     
     for (x = BLOCK_LOW(rank, size, X+1); x <= BLOCK_HIGH(rank, size, X+1); x++) {
      if (x == 0 || x == X) continue;

      for (y=1; y<Y; y++) 
       for (z=1; z<Z; z++) 
	    if (h[x][y][z] > 0) {              
		   
	     Ek=RTONF*log(Kc[x][y][z]/Ki[x][y][z]);
         ENa=RTONF*log(140./Nai[x][y][z]);     
         ECa=(RTONF/2.)*log(2.5/Cai[x][y][z]);

// ** Define the differential approximation ** //
	     dvdx2 = (V[x+1][y][z]+V[x-1][y][z]-2*V[x][y][z])/(HX*HX);
	     dvdy2 = (V[x][y+1][z]+V[x][y-1][z]-2*V[x][y][z])/(HX*HX);
	     dvdz2 = (V[x][y][z-1]+V[x][y][z+1]-2*V[x][y][z])/(HX*HX);

	 switch (h[x][y][z]) {
       case 1: dvdt[x][y][z]=D*(dvdx2+dvdy2+dvdz2)-I_tot_PM(x, y, z); break;
       case 2: dvdt[x][y][z]=D*(dvdx2+dvdy2+dvdz2)-I_tot_CT(x, y, z); break;
       case 3: dvdt[x][y][z]=D*(dvdx2+dvdy2+dvdz2)-I_tot_BB(x, y, z); break; 
       case 4: dvdt[x][y][z]=D*(dvdx2+dvdy2+dvdz2)-I_tot_CS(x, y, z); break;
       case 5: dvdt[x][y][z]=0.5*D*(dvdx2+dvdy2+dvdz2)-I_tot_PV(x, y, z); break; //D is the diffusion coefficient, notice smaller for PV.  now halved to 0.25
       case 6: dvdt[x][y][z]=D*(dvdx2+dvdy2+dvdz2)-I_tot_LA(x, y, z); break;
       case 7: dvdt[x][y][z]=D*(dvdx2+dvdy2+dvdz2)-I_tot_RA(x, y, z); break;
//     default: dvdt[x][y][z]=D*(dvdx2+dvdy2)-I_tot_RA(x, y, z); 
	}

       new_V[x][y][z]=V[x][y][z]+HT*dvdt[x][y][z];
	   }
   }


// ** Replace the old Voltage values with the new ones ** //
  for (x = BLOCK_LOW(rank, size, X+1); x <= BLOCK_HIGH(rank, size, X+1); x++) 
    for (y=0; y<=Y; y++)
	  for (z=0; z<=Z; z++)
	    if (h[x][y][z] > 0) { 
		  V[x][y][z]=new_V[x][y][z];
		}

   if (t>=0.04 && t<(0.04+HT))
	 for (x = BLOCK_LOW(rank, size, X+1); x <= BLOCK_HIGH(rank, size, X+1); x++)
	 {
       for (y=0; y<=200; y++)
		 for (z=0; z<=10; z++)
         {
		   V[x][y][z] = -73.8;
		   new_V[x][y][z] = -73.8;
		   
		 }
	   for (y=400; y<=Y; y++)
	 	 for (z=0; z<=10; z++)
         {
		   V[x][y][z] = -73.8;
		   new_V[x][y][z] = -73.8;
		 }
     }     
	 
	 if (rank < (size-1))
       MPI_Send(V[BLOCK_HIGH(rank, size, X+1)], (Y+1)*(Z+1), MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD);
     if (rank > 0)
       MPI_Recv(V[BLOCK_HIGH((rank-1), size, X+1)], (Y+1)*(Z+1), MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, &status);
     if (rank > 0)
       MPI_Send(V[BLOCK_LOW(rank, size, X+1)], (Y+1)*(Z+1), MPI_FLOAT, rank-1, 1, MPI_COMM_WORLD);
     if (rank < (size-1))
       MPI_Recv(V[BLOCK_LOW((rank+1), size, X+1)], (Y+1)*(Z+1), MPI_FLOAT, rank+1, 1, MPI_COMM_WORLD, &status);

     if (cnt % 1000 == 0) {     //took off a zero
       MPI_Gatherv (&(((float *)V)[recv_disp[rank]]), recv_cnts[rank], MPI_FLOAT, &global_V[0][0][0], recv_cnts, recv_disp, MPI_FLOAT, 0, MPI_COMM_WORLD);
     }

     if (!rank)
 
		 if (cnt == 1000) //ditto 
         { 
	     cnt=0;         
         
	     str = malloc (8*sizeof(char));
	     sprintf (str, "a%d.vtk", cnt2++);
         out = fopen (str, "wt");
         
         fprintf (out, "# vtk DataFile Version 3.0\n");   
         fprintf (out, "vtk output\n");
         fprintf (out, "ASCII\n");
         fprintf (out, "DATASET STRUCTURED_POINTS\n");
         fprintf (out, "DIMENSIONS %d %d %d\n", X+1, Y+1, Z+1);
         fprintf (out, "SPACING 1 1 1\n");
         fprintf (out, "ORIGIN 0 0 0\n");
         fprintf (out, "POINT_DATA %d\n", (X+1)*(Y+1)*(Z+1));
         fprintf (out, "SCALARS ImageFile float 1\n");
         fprintf (out, "LOOKUP_TABLE default\n"); 

		 for (z=0; z<=Z; z++)
		 {
	       for (y=0; y<=Y; y++)
           {
             for (x=0; x<=X; x++)
			   fprintf (out, "%2.2f ", global_V[x][y][z]);
              
	         fprintf (out, "\n");	       
           }
           fprintf (out, "\n");	
		 }
	     fclose (out);
	     free (str);
         } 
   }
   return 0;
}
