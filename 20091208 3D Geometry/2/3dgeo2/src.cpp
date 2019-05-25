#include <stdio.h>
#include <stdlib.h>
#include <math.h>

	const int X = 300; // one side of either left or right atrium
	const int Y = 600; // left and right atrium length of one side
	const int Z = 50; // height of the whole lot
	const int pv_radius = 30; // radius of pulmonary vein
	const int pv_thickness = 5;
	const int offset_x = 250;
	const int offset_y = 200;

	int garray[X+1][Y+1][Z+1]; // the array that will later be outputted to a file

int main()
{
	int x, y, z;

	FILE *in;
	in = fopen("geometry.c", "w");

	for (z=0; z<=10; z++) // layer number one-10, mostly 1s but also with a hole
	{
		for (y=0; y<=Y; y++)
		{
			for (x=0; x<=X; x++)
			{
				if ((((x-offset_x)*(x-offset_x)+(y-offset_y)*(y-offset_y))>=pv_radius*pv_radius))
				{

					//LA/RA
					if (y<300)
					{
						garray[x][y][z]=6;
					}
					else
					{
						garray[x][y][z]=7;
					}

					//PM
					if ((y>350 && y<= 450) && ((x >= 40 && x < 60) || (x >= 80 && x <100) || (x >= 120 && x < 140) || (x >= 160 && x < 180) || (x >= 200 && x <220 ) || (x >= 240 && x < 260)))
					{
						garray[x][y][z] = 1;
                    }

					//CT
					if((y>=325 && y<= 350) && (x>=25 && x<=275))
					{
						garray[x][y][z] = 2;
                    }

					//BB
					if((y>=250 && y<= 350) && (x>=0 && x<25))
					{
						garray[x][y][z] = 3;
					}

					//CS
					if((y>=250 && y<= 350) && (x>=275 && x<X))
					{
						garray[x][y][z] = 4;
					}
				
				}
				else
				{
					garray[x][y][z] = 0;
				}
			}
		}
	}

	for (z=11; z<=Z; z++) // layer number 11-50, mostly 0s but also with a pul vein
	{
		for (y=0; y<=Y; y++)
		{
			for (x=0; x<=X; x++)
			{
				if ((((x-offset_x)*(x-offset_x)+(y-offset_y)*(y-offset_y))>=pv_radius*pv_radius)&&(((x-offset_x)*(x-offset_x)+(y-offset_y)*(y-offset_y))<=(pv_radius+pv_thickness)*(pv_radius+pv_thickness)))
				{
					garray[x][y][z] = 5; // 5 being pul vein
				}
				else
				{
					garray[x][y][z] = 0;
				}
			}
		}
	}

	for (z=0; z<=Z; z++) // all edges are zero
	{
		for (y=0; y<=Y; y++)
		{
			for (x=0; x<=X; x++)
			{
				if (x==0 || x==X || y==0 || y==Y || z==0 || z==Z || x==1 || x==X-1 || y==1 || y==Y-1 || z==1 || z==Z-1)
				{
					garray[x][y][z] = 0;
				}
			}
		}
	}

	for (z=0; z<=Z; z++) // septum is zero
	{
		for (y=0; y<=Y; y++)
		{
			for (x=0; x<=X; x++)
			{
				if (y==300)
				{
					if (garray[x][y][z]!=4 && garray[x][y][z]!=5)
					{
						garray[x][y][z] = 0;
					}
				}
			}
		}
	}


    for (z=0; z<=Z; z++)
	{
		for (y=0; y<=Y; y++) // printing it all
		{
			for (x=0; x<=X; x++)
			{
				fprintf(in, "%d ", garray[x][y][z]);
			}
			fprintf(in, "\n");
		}
		fprintf(in, "\n");
	}


 fclose(in);
 return 0;
}
