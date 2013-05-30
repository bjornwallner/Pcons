#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "nets.h"

 
int read_net(char* filename,network* net)
{
  char	buff[512];	/* Input string */
  char	word[100];	/* PDB file line mode */
  char  tmp[100];
  FILE  *fp;
  int temp;
  long int pos;
  int i=0;
  int j=0;
  int get_w1=0;
  int get_b1=0;
  int get_w2=0;
  int get_b2=0;
  char a;

  //net=(network*)malloc(sizeof(net));


  //net->w1=&temp2;
  //printf("%lf\n",*(net->w1));
  // return net;
  //net->nin=12;
  //printf("%d\n",net->nin);
  
  fp=fopen(filename,"r");	/* Does file exist? */
  if (fp!=NULL)	/* If yes, read in coordinates */
    {
      while(fscanf(fp,"%s",buff)!=EOF)
      {
	
	if(strcmp("nin",buff)==0)
	  {
	    fscanf(fp,"%d",&net->nin);
	    
	    //	    printf("%d",net->nin);
	  }
	  else if(strcmp("nhidden",buff)==0)
	    {
	      if(fscanf(fp,"%d",&net->nhidden)!=1) {
		fprintf(stderr,"error reading nhidden\n");
	      }
	      
	    }
	  else if(strcmp("nout",buff)==0)
	    {
	      fscanf(fp,"%d",&net->nout);
	    }
	  else if(strcmp("w1",buff)==0)
	    {
	      //Store values in nhidden x nin matrix
	      for(j=0;j<net->nin;j++)
		{
		  for(i=0;i<net->nhidden;i++)
		    {
		      fscanf(fp,"%lf",&net->w1[i][j]); 
		    }
		}
	    }
	  else if(strcmp("b1",buff)==0)
	    {
	      for(i=0;i<net->nhidden;i++)
		{
		  fscanf(fp,"%lf",&net->b1[i]);
		  //printf("%lf ",net->b1[i]);
		  //		}
		}
	    }
	  else if(strcmp("w2",buff)==0)
	    {
	      for(i=0;i<net->nhidden;i++)
		{
		  fscanf(fp,"%lf",&net->w2[i]);
		}
	    }
	  else if(strcmp("b2",buff)==0)
	    {
	      fscanf(fp,"%lf",&net->b2);
	    }
//	  //zprintf("%s\n",buff);
      }
      fclose(fp);
//	  printf("nin %d\n",net->nin);
//	  printf("nhidden %d\n",net->nhidden);
//	  printf("nout %d\n", net->nout);
//	  printf("w1\n");
//	  for(j=0;j<net->nin;j++)
//	    {
//	      for(i=0;i<net->nhidden;i++)
//		{
//		  printf("%lf ", net->w1[i][j]); 
//		}
//	      printf("\n");
//	    }
//	  printf("b1\n");
//	  for(i=0;i<net->nhidden;i++)
//	    {
//	      printf("%lf ", net->b1[i]); 
//	    }
//	  printf("\nw2\n");
//	  for(i=0;i<net->nhidden;i++)
//	    {
//	      printf("%lf ", net->w2[i]); 
//	    }
//	  printf("\nb2\n");
//	  printf("%lf\n", net->b2); 
    }
  else
    {
      printf("Couldn't open file %s\n",filename);
      exit(1);
    }
//  for(i=0;i<net->nhidden;i++)
//    {
//	for(j=0;j<net->nin;j++)
//	  {
//	    printf("%6.5lf ",net->w1[i][j]);
//	  }
//	printf("\n");
//    }
  
  return 0;
}

//void free_net(network *net)
//{
//  //  free(net->w1);
//  free(net->b1);
//  free(net->w2);
//  free(net->b2);
//}


double netfwd(double* values,network* net)
{
  int i,j;
  double nodsum=0;
  double output=0;

  for(i=0;i<net->nhidden;i++)
    {
      nodsum=0;
      for(j=0;j<net->nin;j++)
	{
	  nodsum+=net->w1[i][j]*values[j];
	}
      nodsum+=net->b1[i];
      nodsum=tanh(nodsum);
      output+=nodsum*net->w2[i];
    }
  output+=net->b2;
  //  printf("%lf\n",output);
  return output;
}

