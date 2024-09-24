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

double f1; // (double)NP/(double)(MAXX*MAXX);
double f2; // 1-f1;

int* x; //coordiantes of the pixels, with size NP
int* y;

int* SN;//# of site pairs in each bin, with size Nt
int* BN; //for S2, # of black pixel pairs separted at specific dist, with size Nt
int* FN; //# of interface black pixel pairs in each bin


double* S2; //the sampled S2, with size Nt
double* Fss; //the sampled Fss with size Nt

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

  config = new int*[MAXX];
  
 


  for(int i=0; i<MAXX; i++)
    {
      config[i] = new int [MAXX];

      for(int j=0; j<MAXX; j++)
	config[i][j] = 0;
    }


  
  x = new int [NP];
  y = new int [NP];

  SN = new int [Nt];
  BN = new int [Nt];
  FN = new int [Nt];

  S2 = new double [Nt];
  Fss = new double [Nt];


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

  fclose(fp);
}


void Clear_All()
{
  for(int i=0; i<MAXX; i++)
    delete [] config[i];
  delete [] config;

 
  delete [] SN;
  delete [] BN;
  delete [] FN;

  delete [] S2;
  delete [] Fss;
      

 
  delete [] x;
  delete [] y;
  
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

	  N_sur++; //to compute volume fraction of surface pixels for normalization
	}

    }
    
}

void init_data()
{ 
  check_interface();

  //obtain the bin values...
  for(int i=0; i<Nt; i++)
    {
      SN[i] = 0;
      BN[i] = 0;
      FN[i] = 0;
    }

 
  printf("Initializing the site seperation...\n");

  //printf("Here!\n");
  double temp_dist;
  int nd, tempI, tempJ;
  double dd;

  for(int i=0; i<MAXX; i++)//every lattice is equivalent
    for(int j=0; j<MAXX; j++)
      {
	//printf("Here!\n");
	tempI = i;
	tempJ = j;
	
	if(tempI>MAXX/2) tempI = MAXX - tempI;
	if(tempJ>MAXX/2) tempJ = MAXX - tempJ;
	

	temp_dist = sqrt(tempI*tempI+tempJ*tempJ);
	nd = (int)floor(temp_dist);
	dd = temp_dist - nd;
	//printf("nd = %d\n", nd);
	if(dd<0.5)
	  {
	    if(nd<Nt) SN[nd] = SN[nd]+1;
	  }
	else
	  {
	    if(nd<(Nt-1)) SN[nd+1] = SN[nd+1]+1;
	  }
      }
  
  //printf("SN[62] = %d\n", SN[62]);


  //multiplying MAXX*MAXX would be too large for int type, and WILL CAUSE OVERFLOW!!!!!!!
  //must use long int OR take care of this using double type manipulation, as we do here
  /*
  for(int i=0; i<Nt; i++)
    {
      //SN[i] = SN[i]*MAXX*MAXX;

      printf("%d  %d \t %d\n", i, BN[i], SN[i]);
    }
  */
  SN[0] = SN[0]*2;
  //compute and bin all the site pair seperations

  printf("Initializing the pixel seperation...\n");

  for(int i=0; i<NP; i++)
    for(int j=0; j<i+1; j++)
      {
	temp_dist = pixeldistance(i,j);

	nd = (int)floor(temp_dist);
	dd = temp_dist - nd;
	if(dd<0.5)
	  {
	    if(nd<Nt) 
	      {    
		BN[nd] = BN[nd] + 2;

		//only store the pixel separation on the interface...
		if(config[x[i]][y[i]]==2 && config[x[j]][y[j]]==2)
		  FN[nd] = FN[nd] + 2;
	      }
	  }
	else
	  {
	    if(nd<(Nt-1))
	      { 
		BN[nd+1] = BN[nd+1] + 2;
		
		//only store the pixel separation on the interface...
		if(config[x[i]][y[i]]==2 && config[x[j]][y[j]]==2)
		  FN[nd+1] = FN[nd+1] + 2;
	      }
	  }
      }
  //compare and bin the initial configuration of black pixels
}

void sample_S2(double* S)
{
  for(int i=0; i<Nt; i++)
    {
      S[i] = (double)BN[i]/(double)SN[i];

      //printf("%d  %d \t %d\n", i, BN[i], SN[i]/(MAXX*MAXX));

      S[i] = S[i]/(MAXX*MAXX);
    }
}


void sample_Fss(double* F)
{
 for(int i=0; i<Nt; i++)
    {
      F[i] = (double)FN[i]/(double)SN[i];

      F[i] = F[i]/(MAXX*MAXX);
    }
}


void get_autocovar()
{
  f1 = (double)NP/(double)(MAXX*MAXX);
  f2 = 1.0 - f1;

  FILE* fp = fopen("Sampled_AutoCovar.txt","w");
  for(int r=0; r<Nt; r++)
    {
      //printf("%d \t %f \n", r, (S2[r]-f1*f1)/(f1*f2));
      fprintf(fp, "%d \t %f \n", r, (S2[r]-f1*f1)/(f1*f2));
    }
  fclose(fp);
}


void get_scaled_Fss()
{
  f1 = (double)N_sur/(double)(MAXX*MAXX);
  f2 = 1.0 - f1;

  FILE* fp = fopen("Sampled_Scaled_Fss.txt","w");
  for(int r=0; r<Nt; r++)
    {
      //printf("%d \t %f \n", r, (S2[r]-f1*f1)/(f1*f2));
      fprintf(fp, "%d \t %f \n", r, (Fss[r]-f1*f1)/(f1*f2));
    }
  fclose(fp);
}

void print_S2()
{
  FILE* fp = fopen("Sampled_S2.txt","w");
  for(int r=0; r<Nt; r++)
   {
     //printf("%d \t %f \n", r, S2[r]);
     fprintf(fp, "%d \t %f \n", r, S2[r]);
   }
  fclose(fp);
}

void print_Fss()
{
  FILE* fp = fopen("Sampled_Fss.txt","w");
  for(int r=0; r<Nt; r++)
   {
     //printf("%d \t %f \n", r, S2[r]);
     fprintf(fp, "%d \t %f \n", r, Fss[r]);
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

  printf("SN[0] = %d\t BN[0] = %d \n", SN[0], BN[0]);

  printf("*************************************************\n");

  //get_BN();

  //sample_S2(S2);

  //print_S2();

  //get_autocovar();

  sample_Fss(Fss);
   
  print_Fss();

  get_scaled_Fss();

  Clear_All();

}
