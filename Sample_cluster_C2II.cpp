//****************************************************
//****** Compute Two-Point Cluster Function C2 ******
//****************************************************

//The algoritm tells the pixels belongs to which cluster
//The complexity is O(3^d*N^2), d- dimension, N- # of pixels
//The algorithm can be easily generalized to point process

//New function added:
//compute the number size of each cluster
//compute the largest linear size of the clusters

//Jan. 10, 2008
//Modification is made to make sure the "majority rule" works when binning the distances;
//Also, the ordered cluster list is important!!!

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
//#include "const.h"


#define MAXX 180
#define NP 11427 //# of black pixels
#define Nt 80
#define N MAXX*MAXX



int config[MAXX][MAXX];
int x[NP];
int y[NP];
int sitex[N];
int sitey[N];

int cluster[NP];
int local_cluster[NP][9];//store the local cluster situation; the maximum # in 2d is 9, 3^d is the number for dimension d...
int Nc; //total number of clusters


double distance[NP][NP];
double separation[N];//separation of pixels

int SN[Nt];//# of site pairs in each bin
int CN[Nt];//# of black pixel cluster pairs in each bin
int BN[Nt]; //for S2...

double C2[Nt];  
double S2[Nt];

struct node
{
  int index;
  struct node* next;
};

struct node* newlist[NP];//ordered final cluster list

int Lsize[NP];//store the linear size of each cluster
//int Nsize[NP];//store the number size of each cluster

//int LCsite[NP][2];//store the points that have the largest distance

int NclusterX[MAXX];//number of pixels in each slice with the same X for a particular cluster
int NclusterY[MAXX];//number of pixels in each slice with the same Y for a particluar cluster

//***************************************************



void read_config()
{ 
  for(int i=0; i<MAXX; i++)
    for(int j=0; j<MAXX; j++)
      config[i][j] = 0;

  FILE* fp;
  if((fp = fopen("Fconfig.txt","r"))==NULL)
    {
      printf("Can not open Fconfig.txt file! Abort!\n");
      exit(1);
    }

  int indx;
  int indy;

  for(int i=0; i<NP; i++)
    {
      fscanf(fp, "%d", &indx);
      fscanf(fp, "%d", &indy);

      config[indx][indy] = 1;
      x[i] = indx;
      y[i] = indy;
    }

  fclose(fp);
}



double pdistance(int indi, int indj)
{
  double dx = abs(x[indi]-x[indj]);
  if(dx>=MAXX/2) dx = MAXX - dx;

  double dy = abs(y[indi]-y[indj]);
  if(dy>=MAXX/2) dy = MAXX - dy;

  double dd = sqrt(dx*dx+dy*dy);

  return dd;
}

double sdistance(int indi, int indj)
{
  int dxs = abs(sitex[indi] - sitex[indj]);
  if(dxs >= MAXX/2) dxs = MAXX - dxs;
  int dys = abs(sitey[indi] - sitey[indj]);
  if(dys >= MAXX/2) dys = MAXX - dys;

  double ds = sqrt(dxs*dxs + dys*dys);

  return ds;

}



void init_config()
{
  for(int i=0; i<MAXX; i=i+2)
    for(int j=0; j<MAXX; j++)
    {
      x[i*MAXX/2+j] = i;
      y[i*MAXX/2+j] = j;
    }

}


void identi_cluster()
{
  //initialize the local cluster list
  printf("Initializing the local cluster list....\n");

  for(int i=0; i<NP; i++)
    {
      int pt = 0;

      for(int j=0; j<NP; j++)
        {
            double dt = distance[i][j];
            if(dt<2)
	      {
                  local_cluster[i][pt] = j;
        	  pt++;
	      }
	}
      }
  //****************************************************
    
  //Now we search the local cluster list to identify global clusters

  int addpt = 1;//store the # of global clusters

  for(int i=0; i<NP; i++)
    {

      //printf("Calculate the %d the Cluster... \n", i+1);

      int indicator = 0;
      int cluster_ind = 100;

      //*************************************************************
      //first we investigate the list
      for(int j=0; j<9; j++)
	{
	  if(local_cluster[i][j]!= -1)
	    {
	      int k = local_cluster[i][j];

	      //identify a new cluster
	      if(cluster[i] == 0 && cluster[k] == 0)
		{
		  indicator = 1;
		}
	      //for the following two cases, we need to combine small clusters to form bigger ones....
	      else if(cluster[i] == 0 && cluster[k] !=0)
		{
		  indicator = 2;
		  cluster_ind = cluster[k];
		  break;
		}
	      else if(cluster[i]!=0)
		{
		  indicator = 3;
		  cluster_ind = cluster[i];
		  break;
		}
	     
	    }

	  else if(local_cluster[i][j]==-1) break;
	}
     

      //***************************************************************
      if(indicator == 1)
	{
	  for(int j=0; j<9; j++)
	    {
	      if(local_cluster[i][j]!=-1)
		{
		  int k = local_cluster[i][j];
		  cluster[k] = addpt;
		}
	      else if(local_cluster[i][j]==-1) break;
	    }
	  addpt++;

	  //printf("%d \n", addpt);
	}

      else if(indicator == 2)
	{
	  for(int j=0; j<9; j++)
	    {
	      if(local_cluster[i][j]!=-1)
		{
		  int k = local_cluster[i][j];

		  int ci = cluster[k];
                 
		  if(ci!= 0)
		    {
		       for(int m=0; m<NP; m++)
		        {
		          if(cluster[m] == ci) cluster[m] = cluster_ind;
		        }
		    }
		  else
		    {
		      cluster[k] = cluster_ind;
		    }
		  
		}
	      else if(local_cluster[i][j] == -1) break;
	    }
	}
      else if(indicator == 3)
	{
	  for(int j=0; j<9; j++)
	    {
	      if(local_cluster[i][j]!=-1)
		{
		  int k = local_cluster[i][j];
		  int ci = cluster[k];

		  if(ci!=0)
		    {
	       	       for(int m=0; m<NP; m++)
	         	 {
		           if(cluster[m] == ci) cluster[m] = cluster_ind;
		         }
		    }
		  else
		    {
		      cluster[k] = cluster_ind;
		    }
		}
	      else if(local_cluster[i][j] == -1) break;
	    }
	}
    }
  //Now we put them in a better list
  int temp_list[NP];
  for(int i=0; i<NP; i++)
    {
      temp_list[i] = -1;
    }

  //Count how many clusters 
  int pt = 0;

  for(int i=0; i<NP; i++)
    {
      int ck = cluster[i];//indicating belonging to which cluster
      int cj = 0;

      for(int j=0; j<=pt; j++)
	{
	  if(ck == temp_list[j]) break;
	  else cj++;
	}

      if(cj == pt+1)
	{
	  temp_list[pt] = ck;
	  pt++;
	}
    }
  Nc = pt; //get the number of clusters

  //Now we generate the newlist...
  for(int i=0; i<NP; i++)
    {
      newlist[i] = NULL;
    }
  for(int i=0; i<Nc; i++)
    {
      int ck = temp_list[i];
      struct node* pointer = NULL;

      for(int j=0; j<NP; j++)
	{
	  if(cluster[j] == ck)
	    {
	      pointer = (struct node*)malloc(sizeof(struct node));
	      pointer->index = j;
	      pointer->next = newlist[i];
	      newlist[i] = pointer;
	      cluster[j] = i;
	    }
	}
    }
}

/*
void free_list()
{
  for(int i=0; i<Nc; i++)
    {
      struct node* pointer = newlist[i];

      while(pointer!=NULL)
	{
	  pointer = newlist[i];
          newlist[i] = newlist[i]->next;
	  free(pointer);
	}
    }
}
*/

void init_data()
{
 for(int i=0; i<NP; i++)
    {
      cluster[i] = 0;
    }

  for(int i=0; i<NP; i++)
    for(int j=0; j<9; j++)
      {
	local_cluster[i][j] = -1;
      }

 for(int i=0; i<NP; i++)
    for(int j=0; j<=i; j++)
      {
	double d = pdistance(i,j);
	distance[i][j] = distance[j][i] = d;
      }

 for(int i=0; i<NP; i++)
   {
     newlist[i] = NULL;
   }

 //************************************************
//compute and bin all the site pair seperations
 for(int i=0; i<Nt; i++)
    {
      SN[i] = 0;
      CN[i] = 0;
      BN[i] = 0;
    }

  for(int i=0; i<MAXX; i++)
    for(int j=0; j<MAXX; j++)
      {
	sitex[MAXX*i+j] = j;
	sitey[MAXX*i+j] = i;
      }

  printf("Initializing the site seperation...\n");

  for(int i=0; i<N; i++)//every lattice is equivalent
      {
	separation[i] = sdistance(i,0);

	int nd = floor(separation[i]);
	double dd = separation[i] - nd;
	if(dd<0.5)
	  {
	     if(nd<Nt) SN[nd] = SN[nd] + 1;
	  }
	else
	  {
	    if(nd<Nt-1) SN[nd+1] = SN[nd+1] + 1;
	  }
      }

  for(int i=0; i<Nt; i++)
    {
      SN[i] = SN[i]*N;
    }
  SN[0] = SN[0]*2;

 //*************************************************
}

void get_CN()
{
  int Plist[NP];
  for(int i=0; i<NP; i++)
    {
      Plist[i] = -1;
    }

  for(int i=0; i<Nc; i++)
    {
      struct node* pointer = newlist[i];
      int counter = 0;

      while(pointer!=NULL)
	{
	  Plist[counter] = pointer->index;
	  pointer = pointer->next;
	  counter++;
	}
      for(int j=0; j<counter; j++)
	for(int k=0; k<j+1; k++)
	  {

	    double dpp = distance[Plist[j]][Plist[k]];
	    int ndpp = floor(dpp);
	    double ddpp = dpp - ndpp;
	    if(ddpp<0.5)
	      {
		if(ndpp<Nt) CN[ndpp] = CN[ndpp]+2;
	      }
	    else
	      {
		if(ndpp<Nt-1) CN[ndpp+1] = CN[ndpp+1]+2; 
	      }
	  }

      for(int j=0; j<counter+2; j++)
	{
	  Plist[j] = -1;//resume the list for further use
	}

    }

}




void sample_C2()
{
  for(int r=0; r<Nt; r++)
    {
      C2[r] = (double)CN[r]/(double)SN[r];
    }
}


void get_BN()
{
 printf("Initializing the pixel seperation...\n");

  for(int i=0; i<NP; i++)
    for(int j=0; j<i+1; j++)
      {
	//distance[i][j] = distance[j][i] = pixeldistance(i,j);

	int nd = floor(distance[i][j]);
	double dd = distance[i][j] - nd;
	if(dd<0.5)
	  {
	     if(nd<Nt) BN[nd] = BN[nd] + 2;
	  }
	else
	  {
	    if(nd<Nt-1) BN[nd+1] = BN[nd+1] + 2;
	  }
      }
}


void sample_S2()
{
  for(int r=0; r<Nt; r++)
    {
      S2[r] = (double)BN[r]/(double)SN[r];
    }
}


int max(int x, int y)
{
  if(x>y) return x;
  else return y;
}

void get_Lsize()
{
  for(int i=0; i<NP; i++)
    Lsize[i] = 0;


  for(int i=0; i<Nc; i++) //loop over all clusters
    {
      double ls = 0;

      for(int j=0; j<MAXX; j++)
	{
	  NclusterX[j] = 0;
	  NclusterY[j] = 0;
	}

      struct node* pt1 = newlist[i];

      while(pt1 != NULL)
	{
	  int indm = pt1->index;
	  pt1 = pt1->next;

	  NclusterX[x[indm]]++;
	  NclusterY[y[indm]]++;
	 
	}//loop over pt1

      int lx = 0;
      int ly = 0;

      for(int j=0; j<MAXX; j++)
	{
	  if(NclusterX[j]>0) lx++;
	  if(NclusterY[j]>0) ly++;
	}

      Lsize[i] = max(lx, ly);
    }
}


double cluster_size()
{
  double lcs = 0;

  for(int i=0; i<Nc; i++)
    {
      if(Lsize[i]>lcs) lcs = Lsize[i];
    }

  return lcs;
}


void print_cluster(int indm)
{
  FILE* fp = fopen("cluster.txt","w");

  struct node* pt = newlist[indm];

  while(pt!=NULL)
    {
      int indn = pt->index;
      pt = pt->next;

      fprintf(fp, "%d \t %d\n", x[indn], y[indn]);
    }

  fclose(fp);
}



//need to further check the criterion
int percolate()
{
  double lcs = 0;
  int index;
  int a = 0;
  int b = 0;

  for(int i=0; i<Nc; i++)
    {
      if(Lsize[i]>lcs)
	{
	  lcs = Lsize[i];
	  index = i;
	}
    }

  print_cluster(index);

  if(lcs>(MAXX-1)) b++;

  if(b==1) return 1;
  else return 0;
		     
}

void print_C2()
{
  FILE* fp = fopen("Sampled_C2.txt","w");
  for(int r=0; r<Nt; r++)
   {
     printf("%d \t %f \n", r, C2[r]);
     fprintf(fp, "%d \t %f \n", r, C2[r]);
   }
  fclose(fp);
}


void print_S2()
{
  FILE* fp = fopen("Sampled_S2.txt","w");
  for(int r=0; r<Nt; r++)
   {
     printf("%d \t %f \n", r, S2[r]);
     fprintf(fp, "%d \t %f \n", r, S2[r]);
   }
  fclose(fp);
}


void print_IndCluster()
{
  FILE *fp = fopen("identi_cluster.txt","w");
  fprintf(fp, "Cluster-Point Identification List: \n");

  for(int i=0; i<NP; i++)
   {
     printf("Pixel: %d  =>  Cluster: %d \n", i, cluster[i]);

     fprintf(fp, "Pixel: %d  =>  Cluster: %d \n", i, cluster[i]);

   }

  fclose(fp);
}

void print_LocalCluster()
{
  FILE* fp = fopen("Local_Cluster.txt","w");
  fprintf(fp, "Local Cluster List: \n");
 
  for(int i=0; i<NP; i++)
   {
       printf("Local Cluster: %d \t ", i);
       fprintf(fp, "Local Cluster: %d \t ", i);

       for(int j=0; j<9; j++)
        {
	  printf("-> %d ", local_cluster[i][j]);
          fprintf(fp, "-> %d ", local_cluster[i][j]);
        }

       printf("\n");
       fprintf(fp, "\n");
   }

 /*
 for(int i=0; i<NP; i++)
   {
     printf("distance[1][%d] = %f \n", i, distance[1][i]);
   }
 */

  fclose(fp);
}


void print_config()
{
  FILE* fp = fopen("config.txt","w");
  for(int i=0; i<NP; i++)
   {
     fprintf(fp, "%d \t %d \n", x[i], y[i]);
   }
  fclose(fp);
}

void print_newlist()
{
  printf("***************************************************\n");

  printf("The number of clusters is:%d \n", Nc);

  FILE* fp = fopen("new_list.txt","w");
  fprintf(fp, "The number of clusters is:%d \n", Nc);

  for(int i=0; i<Nc; i++)
   {

     printf("Cluster %d := { ", i);
     fprintf(fp,"Cluster %d := { ", i);

     struct node* pointer = newlist[i];
     while(pointer!=NULL)
       {
	 int ind = pointer->index;
	 printf(" %d ", ind);
	 fprintf(fp, " %d ", ind);

	 pointer = pointer->next;
       }

     printf(" } \n");
     fprintf(fp, " } \n");
   }

  fclose(fp);
}

main()
{
  //First, need to do a lot of initialization....
  //***********************************************
  //init_config();

  read_config();

  init_data();

  //**********************************************

  identi_cluster();

  get_CN();

  sample_C2();

  print_C2();

  get_BN();

  sample_S2();

  print_S2();

  //get_Lsize();

 //Print out the results...
 //************************************************

  //double CS = cluster_size();
 
  /*
  //print the linear size of the cluster....
  FILE* fp = fopen("size.txt","a");
  fprintf(fp, "%f \t %d \n", CS, percolate());
  fclose(fp);
  */

  //printf("The linear size of the largest cluster is: %f; Percolate = %d\n", CS, percolate());

  //print_newlist();

  print_config();

 //*************************************************

}
