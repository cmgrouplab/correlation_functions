//***************************************************************
//Sampling Fss using lattice-point method from Large Finite Image
//***************************************************************

//updated Nov, 30, 2011....
//we still use periodic boundary conditions, otherwise the lattice-point method doesn't work...

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


int MAXX; //this is the linear size of the image
int Nt; //the sample length, since no periodic condition is used, Nt is not necessarily 0.5*MAXX
int NP; //the total number of black pixels

int N_sur = 0; //number of surface pixels
int N_void;

double f1; // (double)NP/(double)(MAXX*MAXX);
double f2; // 1-f1;

int* x; //coordiantes of the black pixels, with size NP
int* y;

int* x_v; //coordinates of the void pixels 
int* y_v;

int* x_s; //coordinates of the surface pixels
int* y_s;

double* mindist; //size N_void; the value of min dist from void pixel to surface pixel

int* PFN; //# void pxiel to the closest interface black pixel, size Nt
double* P; //the pore-size distribution function

int** config; //the digitized configuration


//*******************************************

void read_config()
{ 
 
  FILE* fp;
  if((fp = fopen("bconfig.txt","r"))==NULL)
    {
      printf("Can not open Tconfig.txt file! Abort!\n");
      exit(1);
    }

  fscanf(fp, "%d", &MAXX); //read in image size
  fscanf(fp, "%d", &NP); //read in # of black pixels

  MAXX = MAXX+1; //the actual square image is 1 pixel large than MAXX...
  Nt = (int)floor(MAXX/2); //sample length

  N_void = MAXX*MAXX-NP;

  config = new int*[MAXX];
  mindist = new double[N_void];
 

  for(int i=0; i<MAXX; i++)
    {
      config[i] = new int [MAXX];

      for(int j=0; j<MAXX; j++)
	config[i][j] = 0;
    }
  
  x = new int [NP];
  y = new int [NP];

  x_s = new int [NP];
  y_s = new int [NP];

  x_v = new int [N_void];
  y_v = new int [N_void];

  PFN = new int [Nt];

  P = new double [Nt];
  
  int indx;
  int indy;

  //printf("Here!\n");

  for(int i=0; i<NP; i++)
    {

      fscanf(fp, "%d", &indx);
      fscanf(fp, "%d", &indy);

      x[i] = indx;
      y[i] = indy;

      config[indx][indy] = 1;
    }

  
  int temp_void_ct = 0;

  for(int i=0; i<MAXX; i++)
    for(int j=0; j<MAXX; j++)
      {
	if(config[i][j] == 0)
	  {
	    x_v[temp_void_ct] = i;
	    y_v[temp_void_ct] = j;

	    temp_void_ct++;
	  }
      }

  if(temp_void_ct != N_void)
    {
      printf("temp_void_ct = %d; N_void = %d\n", temp_void_ct, N_void);
      exit(1);
    }

  fclose(fp);
}


void Clear_All()
{
  for(int i=0; i<MAXX; i++)
    delete [] config[i];
  delete [] config;

  delete [] PFN;
  delete [] mindist;

  delete [] P;
      
  delete [] x;
  delete [] y;
  delete [] x_v;
  delete [] y_v;
  delete [] x_s;
  delete [] y_s;
  
}


//simply the real distance, no periodic boundary condition is used
double pixeldistance(int indi, int indj)
{
  int dx = abs(x[indi] - x[indj]);
  if(dx > MAXX/2) dx = MAXX - dx;
 
  int dy = abs(y[indi] - y[indj]);
  if(dy > MAXX/2) dy = MAXX - dy;

  double d = sqrt(dx*dx + dy*dy);

  return d;
}


//get the distance between a void pixel and a surface pixel
double pixeldistance_v2s(int indi, int indj)
{
  int dx = abs(x_v[indi] - x_s[indj]);
  if(dx > MAXX/2) dx = MAXX - dx;
 
  int dy = abs(y_v[indi] - y_s[indj]);
  if(dy > MAXX/2) dy = MAXX - dy;

  double d = sqrt(dx*dx + dy*dy);

  return d;
}


void check_interface()
{
  for(int i=0; i<NP; i++)
    {
      int temp_xind;
      int temp_yind;

      int SE = 0;

      for(int m=-1; m<=1; m++)
	for(int n=-1; n<=1; n++)
	  {
	    temp_xind = x[i] + m;
	    if(temp_xind>=MAXX) temp_xind = temp_xind-MAXX;
	    else if(temp_xind<0) temp_xind = temp_xind+MAXX;
	    temp_yind = y[i] + n;
	    if(temp_yind>=MAXX) temp_yind = temp_yind-MAXX;
	    else if(temp_yind<0) temp_yind = temp_yind+MAXX;

	    if(config[temp_xind][temp_yind]>0) SE++;

	  }

      //this means the pixel is on the interface ...
      if(SE<9) 
	{
	  config[x[i]][y[i]] = 2;

	  x_s[N_sur] = x[i];
	  y_s[N_sur] = y[i];

	  N_sur++; //to compute volume fraction of surface pixels for normalization
	}

    }
    
}


double Getmind()
{
  for(int i=0; i<N_void; i++)
    {
      double min_d = 2.0*MAXX;
      double temp_d;

      for(int j=0; j<N_sur; j++)
	{
	  temp_d = pixeldistance_v2s(i, j);

	  if(temp_d <min_d)
	    min_d = temp_d;
	}

      mindist[i] = min_d;
    }
}

void init_data()
{ 
  check_interface();

  printf("Computing the minimum void-to-surface distances...\n");

  Getmind();

  //obtain the bin values...
  for(int i=0; i<Nt; i++)
    {
      PFN[i] = 0;
    }

  for(int i=0; i<N_void; i++)
    {
      if(mindist[i]>=0) //only bin the meaningful distances...
	{
  	  int nd = (int)floor(mindist[i]);
          double dd = mindist[i] - nd;
  	  if(dd<0.5)
	    {
	      if(nd<Nt) 
		{    
                  PFN[nd] = PFN[nd] + 2;
		}
	      
	    }
	  else
	    {
	      if(nd<(Nt-1))
		{ 
                  PFN[nd+1] = PFN[nd+1] + 2;
		}
	    }
	}
    }
 
  //compare and bin the initial configuration of black pixels
}

void sample_P(double* P)
{
  for(int i=0; i<Nt; i++)
    {
      P[i] = (double)PFN[i]/(double)N_void;

      //printf("%d  %d \t %d\n", i, BN[i], SN[i]/(MAXX*MAXX));
    }
}



void print_P()
{
  FILE* fp = fopen("Sampled_P.txt","w");
  for(int r=0; r<Nt; r++)
   {
     //printf("%d \t %f \n", r, S2[r]);
     fprintf(fp, "%d \t %f \n", r, P[r]);
   }
  fclose(fp);
}


main()
{
  //First, need to do a lot of initialization....
  //***********************************************
  //init_config();

  read_config();

  //printf("Here!\n");

  init_data();

  printf("PFN[0] = %d\t PFN[1] = %d \n", PFN[0], PFN[1]);

  printf("*************************************************\n");
 
  sample_P(P);
   
  print_P();

  Clear_All();

}
