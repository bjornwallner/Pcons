#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <sys/types.h>
#include <dirent.h>
#include <sys/times.h>
#include <time.h>
#include <limits.h>
#include <math.h>
#ifdef _OPENMP		      
#include <omp.h>
#endif
//#include "/work/bjornw/source/c/pdb/molecule.h"
//#include "/work/bjornw/source/c/pdb/quality_measure/lgscore.h"
//#include "/afs/pdc.kth.se/home/b/bjornw/source/c/pdb/molecule.h"
//#include "/afs/pdc.kth.se/home/b/bjornw/source/c/pdb/quality_measure/lgscore.h"
#include "molecule.h"
#include "lgscore.h"

//#include "/afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/quality_molecule.h"
#define MAXFILE 10000
#define STRING_BUFFER 1000
void usage();
int SameMethod(char *file1, char *file2);
int Rank(char *file);
int ispdb(char *filename);
char *basename(char *fname);

int main(int argc,char *argv[])             /* Main routine */
{
  DIR            *dip;
  FILE           *fp;
  struct dirent  *dit;
  char           dir[PATH_MAX]="";
  char           file2[PATH_MAX];
  char           tempfilename[PATH_MAX]="";
  char           infilename[PATH_MAX]="";
  char           matrix_outfile[PATH_MAX]="";
  int            rank[MAXFILE];
  // char           filenames[MAXFILE][100];
  char           *filenames_with_path[MAXFILE]={NULL}; //[PATH_MAX];
  char           *filenames[MAXFILE]={NULL}; //[MAXFILE][100];
   // char           **filenames; //[MAXFILE]; //[MAXFILE][100];
  //char           **filenames_with_path='\0';//[MAXFILE][STRING_BUFFER];
  int            i,j,k,l;
  int            files=0;
  int            L=4;
  int            step=2;
  double         w1=0; //relative weight for first ranked 
  double         minsim=121.0;
  double         d0=2.23606798;  //sqrt(5);
  double         factor=0.5;
  double         rank_weight1=1.0;
  double         rank_weight2=1.0;
  char           temp[PATH_MAX];
  //  lgscore        LG[1];
  //  double         LG_average[MAXFILE]={0};
  //double         S_average[MAXFILE]={0};
  //double         Sstr[MAXFILE][MAXRES]={{0}};
  //  int            number_of_comparisons[MAXFILE]={0};

	  
  double         *LG_average,*S_average,*LG_average1,*S_average1,*pcons_google,*pcons_google_previous,*google_weight;
  double          **Sstr,**sim_mat;
  // double          **sim_mat;
  //double *Sstr[MAXFILE]={NULL};
  double            *number_of_comparisons,*number_of_comparisons1;
  int            total_number_of_comparisons=0;
  int nofree=0;
  int            maxlen=0;
  int            normL=1;
  int            userdef_len=-1;
  int            use_rank_weight=0; 
  int            rank1=0; 
  double         rmsd=0;
  double        sqrt5=2.23606798;  //sqrt(5);
  char          target[STRING_BUFFER]="T0XXX";
  int           int_buffer;
  int           max_int;
  int           verbose=0;
  int           fastmode=0;
  int           memorymode=1;
  int           bytes_freed=0;
  int           target_defined_by_user=0;
  int           caspoutput=0;
  int           lgscoreoutput=0;
  int           superimpose_all=0;
  int           force_compare_all=0;
  int           only_pairwise=0;
  int           google=0;
  int           google_iter=0;
  double        google_weight_cut=-1;
  int           output_similarity_matrix=0;
  int           maxfile=MAXFILE;
  int           temp_counter=0;
  //molecule      m[MAXFILE];
  dyn_molecule  *dm;
  int           ignore_res[2000]={0};
  char          ignore_res_file[PATH_MAX]="undef";
  //ignore_res=(int*)malloc(sizeof(int)*maxlen+1);
  //for(i=0;i<=maxlen;i++)
  //  ignore_res[i]=0;

  static clock_t st_time;
  static clock_t en_time;
  static struct tms st_cpu;
  static struct tms en_cpu;
  /* Parse command line for PDB files and options */
  //i=1;
  
  
  

  //printf("%d %d %d %d\n",sizeof(dm),sizeof(dyn_molecule),sizeof(atm),sizeof(double));
  //exit(0);
  
  i=1;
  if(argc > 2)
    {

      while (i<argc)
	{
	  if(strcmp(argv[i],"-d")==0)
	    {
	      i++;
	      strcpy(dir,argv[i]);
	      if(dir[strlen(dir)-1] != '/')
		{
		  strcat(dir,"/");
		}
	      if ((dip = opendir(dir)) == NULL)
		{
		  perror("opendir");
		  exit(0);
		}
	      //printf("Opened \"%s\" directory stream\n",dir);
	      
	      /*  struct dirent *readdir(DIR *dir);
	       *
	       * Read in the files from argv[1] and print */
	      files=0;
	      while ((dit = readdir(dip)) != NULL)
		{
		  //  printf("%s\n", dit->d_name);
		  //	  strcpy(filenames[files],dir);
		  if(dit->d_name[0] != '.')
		    {
		      
		      strcpy(tempfilename,dir);
		      strcat(tempfilename,dit->d_name);
		      
		      if(ispdb(tempfilename))
			{
			  //printf("Allocating memory for %d files\n",files);
			  //filenames=(char**)realloc((char**)filenames,sizeof(char*)*(files+1)); //(files+1));
			  //printf("done\n");
			  //			  filenames[files]=(char*)realloc(filenames[files],sizeof(char)*(strlen(dit->d_name)+1));
			  //printf("Allocating memory for %ld chars\n",sizeof(char)*100);


			  filenames[files]=(char*)realloc(filenames[files],sizeof(char)*PATH_MAX); //(strlen(dit->d_name)+1));
			  filenames_with_path[files]=(char*)realloc(filenames_with_path[files],sizeof(char)*PATH_MAX); //(strlen(dit->d_name)+1));
			  //printf("done\n");
			  //filenames[files]=(char*)realloc(filenames[files],sizeof(char)*20000);
			  //printf("%s\n","first");
			  //
			  //filenames_with_path=(char**)realloc(filenames_with_path,sizeof(char)*20); //(files+1));
			  //// filenames_with_path[files]=(char*)realloc(filenames_with_path[files],sizeof(char)*(strlen(tempfilename)+1));
			  //filenames_with_path[files]=(char*)realloc(filenames_with_path[files],sizeof(char)*2000);
			  //
			  
			  
			  // filenames_with_path[files]=malloc(sizeof(char)*PATH_MAX);
			  strcpy(filenames[files],dit->d_name);
			  //strncpy_NULL(filenames[files],dit->d_name,(strlen(dit->d_name)));

//void strncpy_NULL(char *dest, char *src, size_t n)
//{
//  strncpy(dest, src, n);
//  dest[n]='\0';
//}
			  
			  rank[files]=Rank(filenames[files]);
			  //strncpy_NULL(filenames_with_path[files],tempfilename,(strlen(tempfilename)));
			  strcpy(filenames_with_path[files],tempfilename);
			  //			   printf("%s %d\n",filenames[files],rank[files]);
			  //strcat(filenames_with_path[files],dir);
			  //strcat(filenames_with_path[files],filenames[files]);
			  if(filenames[files][0] == 'T' &&
			     filenames[files][1] == '0' &&
			     !target_defined_by_user)
			    {
			      strncpy(target,filenames[files],5);
			      target[5]='\0';
			    }
		  
			  files++;
			  if(files>MAXFILE)
			    {
			      fprintf(stderr,"Maximum number of files %d exceded. Increase MAXFILE.\n",MAXFILE);
			      exit(1);
			    }

			}
		    }
		} 
	      //exit(1);
	   

	      if (closedir(dip) == -1)
		{
		  perror("closedir");
		  return 0;
		}
	    }
	  if(strcmp(argv[i],"-i")==0)
	    {
	      i++;
	      
	      strcpy(infilename,argv[i]);
	      // printf("%s\n",infilename);
	      fp=fopen(infilename,"r");
	      if(fp!=NULL)
		{
		  files=0;
		  //while(fscanf(fp,"%s",tempfilename)!=EOF)

		  while(fscanf(fp,"%s",tempfilename)!=EOF)
		    {
		      //printf("%s %d\n",tempfilename,strlen(tempfilename));
		      if(ispdb(tempfilename))
			{

			  filenames[files]=(char*)realloc(filenames[files],sizeof(char)*PATH_MAX); //(strlen(dit->d_name)+1));
			  filenames_with_path[files]=(char*)realloc(filenames_with_path[files],sizeof(char)*PATH_MAX);
			  
			  strcpy(filenames_with_path[files],tempfilename);
			  // printf("%s %d\n",filenames_with_path[files],strlen(filenames_with_path[files]));
			  ///			  printf("%s %d\n",filenames_with_path[files],strlen(filenames_with_path[files]));
			  
			  strcpy(filenames[files],basename(tempfilename));
			  rank[files]=Rank(filenames[files]);
			  //			  printf("%s %d\n",filenames[files],rank[files]);
			  // printf("%s\n", basename(tempfilename));
			  //			  filenames[files]=tempfilename[basename(tempfilename)];

			    // strcpy(filenames[files],tempfilename);
			//		strcpy(filenames[i],buff);
		      //tmp=ftell(fp);

			  //filenames[files];
			  //  printf("%s %s %d %s\n",filenames_with_path[files],filenames[files],files,basename(tempfilename));
		      

			  files++;
			}
		    }
		  //fclose(fp);
		}
	      else
		{
		  printf("Cannot open %s",tempfilename);
		  exit(0);
		}
	      fclose(fp);
	    }
	  if(strcmp(argv[i],"-m")==0)
	    {
	      i++;
	      strcpy(matrix_outfile,argv[i]);
	      output_similarity_matrix=1;
	      force_compare_all=1;
	    }
	  if(strcmp(argv[i],"-L")==0)
	    {
	      i++;
	      userdef_len=atoi(argv[i]);
	    }
	    
	  if(strcmp(argv[i],"-t")==0)
	    {
	      i++;
	      target_defined_by_user=1;
	      strcpy(target,argv[i]);
	    }
	   if(strcmp(argv[i],"-f")==0)
	    {
	        fastmode=1;
	        L=4;
		minsim=289.0;
		factor=2;
	    }
	   if(strcmp(argv[i],"-m")==0)
	    {
	      i++;
	      
	      memorymode=atoi(argv[i]);
	      if(memorymode!=0 || memorymode!=1)
		{
		  fprintf(stderr,"Illegal -m switch only -m 0 or -m 1 is allowed!\n\n");
		  usage();
		}

	      
	    }
	   if(strcmp(argv[i],"-A")==0)
	     {
	       force_compare_all=1;
	     }
	   if(strcmp(argv[i],"-only_pairwise")==0)
	     {
	       only_pairwise=1;
	     }
	   if(strcmp(argv[i],"-google")==0)
	     {
	       google=1;
	       google_iter=1;
	     }
	   if(strcmp(argv[i],"-google_iter")==0)
	    {
	      i++;
	      google_iter=atoi(argv[i]);
	    }
	   if(strcmp(argv[i],"-google_weight_cut")==0)
	    {
	      i++;
	      google_weight_cut=atof(argv[i]);
	    }
	   if(strcmp(argv[i],"-casp")==0)
	    {
	      caspoutput=1;
	    }
	   if(strcmp(argv[i],"-lgscore")==0)
	    {
	      lgscoreoutput=1;
	    }
	   if(strcmp(argv[i],"-superimpose_all")==0)
	    {
	      superimpose_all=1;
	    }
	    if(strcmp(argv[i],"-use_rank_weight")==0)
	    {
	      use_rank_weight=1;
	    }
	     if(strcmp(argv[i],"-w1")==0)
	    {
	      i++;
	      w1=atof(argv[i]);
	      if(w1>1){
		fprintf(stderr,"-w1 should be between 0 and 1!\n");
		exit(0);
	      }
		
	    }
	    if(strcmp(argv[i],"-rank1")==0)
	    {
	      rank1=1;
	      //	      w1=1;
	    }
	    if(strcmp(argv[i],"-nofree")==0)
	    {
	      nofree=1;
	      //	      w1=1;
	    }
	   if(strcmp(argv[i],"-normL")==0)
	    {
	       normL=1;
	    }
	   if(strcmp(argv[i],"-v")==0)
	    {
	      printf("Verbose mode on start talking...\n");
	      verbose=1;
	    }
	   if(strcmp(argv[i],"-L0")==0)
	    {
	      i++;
	      L=atoi(argv[i]);
	    }
	   if(strcmp(argv[i],"-step")==0)
	    {
	      i++;
	      step=atoi(argv[i]);
	    }
	    if(strcmp(argv[i],"-factor")==0)
	    {
	      i++;
	      factor=atof(argv[i]);
	    }
	    if(strcmp(argv[i],"-minsim")==0)
	    {
	      i++;
	      minsim=atof(argv[i]);
	    }
	    if(strcmp(argv[i],"-ignore_res")==0) 
	      {

		i++;
		strcpy(ignore_res_file,argv[i]);
		//printf("Will use %s to ignore res\n",ignore_res_file);
		fp=fopen(ignore_res_file,"r");
		if(fp!=NULL)
		  {
		    max_int=0;
		    //while(fscanf(fp,"%s",tempfilename)!=EOF)
		    while(fscanf(fp,"%d\n",&int_buffer)!=EOF)
		      {
			//printf("res : %d\n",int_buffer);
			//if(int_buffer>max_int) {
			//	ignore_res=(int*)realloc(ignore_res,sizeof(int)*int_buffer+1); 
			//	max_int=int_buffer;
			//}
			
			ignore_res[int_buffer]=1;
		      }
		  }
		else
		  {
		    printf("Cannot open %s",ignore_res_file);
		    exit(0);
		  }
		fclose(fp);
	      }
	    i++;
	}
      
      // Done readin input arguments:

      // printf("%s\n",matrix_outfile);

    
      st_time=times(&st_cpu);
      
      // Read all CA coordinates to memory
      if(memorymode)
	{
	  //	  printf("%d",files);
	  dm=malloc(sizeof(dyn_molecule)*files);
	  //for(i=0;i<totfiles;i++)
//	{
//	  dm[i].sequence=malloc(sizeof(char)*residues+1);
//	  dm[i].ss=malloc(sizeof(char)*residues+1);
//	  dm[i].atm=malloc(sizeof(atm)*residues);
//	}
	  
	  for(i=0;i<files;i++)
	    {
	      // strcpy(m[i].filename,filenames_with_path[i]);
	      strcpy(dm[i].filename,filenames_with_path[i]);
	      ///printf("%s\n",dm[i].filename);
	      if(read_molecules_dynamic(&dm[i],'c',ignore_res)==0)
	      	{
		  //printf("%s %d\n",dm[i].filename,dm[i].residues);
		  if(dm[i].residues>0) {
		    if(dm[i].atm[dm[i].residues-1].resnum>maxlen)
		      {
			maxlen=dm[i].atm[dm[i].residues-1].resnum;
		      }
		  }
		  //if(dm[i].method!=NULL)
		  //  {
		  //    printf("%s\n",dm[i].method);
		  //  }
		}
	      else
		{
		  fprintf(stderr,"error reading %s\n",dm[i].filename);
		  exit(0);
		}
	      //else
	      //	{
	      //	  fprintf(stderr,"%s\n",m[i].filename);
	      //}
	      
	    }
	}
      
      //      exit(0);


//int residues=100;
//int totfiles=;
//
//dm=malloc(sizeof(dyn_molecule)*totfiles);
//for(i=0;i<totfiles;i++)
//	{
//	  dm[i].sequence=malloc(sizeof(char)*residues+1);
//	  dm[i].ss=malloc(sizeof(char)*residues+1);
//	  dm[i].atm=malloc(sizeof(atm)*residues);
//	}


      //printf("At most %d comparisons needs to be done...\n",(files-1)*files/2);

      if(verbose)
	printf("Allocating memory for %d files...\n", files);

      //ignore_res

    
      
      
      //  for(i=0;i<=maxlen;i++)
      //	printf("%d %d\n",i,ignore_res[i]);
      LG_average=malloc(sizeof(double)*files);
      S_average=malloc(sizeof(double)*files);
      LG_average1=malloc(sizeof(double)*files);
      S_average1=malloc(sizeof(double)*files);
      number_of_comparisons=malloc(sizeof(double)*files);
      number_of_comparisons1=malloc(sizeof(double)*files);
      Sstr=(double**)malloc(sizeof(double*)*files);
      sim_mat=(double**)malloc(sizeof(double*)*files);


      // printf("%d %d\n",sizeof(double*),sizeof(double));
      //exit(1);
      
      //if(userdef_len>maxlen)
	//{
	  //maxlen=userdef_len;
	  //	}
      //  printf("%d\n",maxlen);
      //exit(1);
      for(i=0;i<files;i++)
	{
	  LG_average[i]=0;
	  S_average[i]=0;
	  LG_average1[i]=0;
	  S_average1[i]=0;
	  number_of_comparisons[i]=0;
	  number_of_comparisons1[i]=0;
	  Sstr[i]=(double*)malloc(sizeof(double)*maxlen);
	  sim_mat[i]=(double*)malloc(sizeof(double)*files);
	  //printf("%d\n",i);
	  for(j=0;j<maxlen;j++)
	    {
	      //printf("%d %d\n",i,j);
	      Sstr[i][j]=0;
	    }
	  for(j=0;j<files;j++)
	    {
	      //printf("%d %d\n",i,j);
	      sim_mat[i][j]=-1;
	    }

	}
      
      //      printf("#Maximum number of comparisons: %d (probably lower if there are models from same method)\n",(int)files*(files-1)/2);
      if(output_similarity_matrix && files>1)
	{
	  fp=fopen(matrix_outfile,"w");
	  if(fp==NULL)
	    {
	      fprintf(stderr,"Couldn't open matrix otufile for writing (%s)\n",matrix_outfile);
	      exit(0);
	    }
	}
#ifdef _OPENMP
      double stime=omp_get_wtime();
#else 
      double stime=(double)clock() / (double)CLOCKS_PER_SEC;
#endif

  
#pragma omp parallel for default(shared) firstprivate(rank_weight1,rank_weight2) shared(Sstr,sim_mat,LG_average,LG_average1,S_average,S_average1,number_of_comparisons,number_of_comparisons1,rank,dm) schedule(dynamic)// private(LG)
       for(int i=0;i<files;i++)
	{
//#ifdef _OPENMP
//int threads = omp_get_num_threads();
//#else
//int threads = 1
//#endif
//	  printf("Running on %d threads...\n",threads);
	  lgscore LG[1];
	  //	  lgscore *LG;
	  //LG=malloc(sizeof(lgscore));
	  if(verbose)
	    printf("%d / %d\n",i+1,files+1);
	  if(output_similarity_matrix)
	    fprintf(fp,"%-40s ",filenames[i]);
	  
	  for(int j=i+1;j<files;j++)
	    {
	      // printf("%s %s %s %s %d %d same:%d\n",dm[i].method,dm[j].method,filenames[i],filenames[j],i,j,!SameMethod(filenames[i],filenames[j]));
	      
	      //	      if((dm[i].method==NULL || dm[j].method==NULL) && (force_compare_all || !SameMethod(filenames[i],filenames[j])) ||
	      //	 (dm[i].method!=NULL && dm[j].method!=NULL) && (force_compare_all || strcmp(dm[i].method,dm[j].method)!=0))

	      if(dm[i].residues>0 && dm[j].residues>0) {
		if(force_compare_all || !SameMethod(filenames[i],filenames[j]))
		  {
		    //printf("IN!\n");
		    //printf("%s %s %d %d\n",dm[i].method,dm[j].method,strcmp(dm[i].method,dm[j].method),SameMethod(filenames[i],filenames[j]));
		    //printf("%d Comparing %s %s %s %s %d %d\n",total_number_of_comparisons,filenames[i],filenames[j],dm[i].filename,dm[j].filename,dm[i].residues,dm[j].residues);
		    rank_weight1=1;
		    rank_weight2=1;
		    if(use_rank_weight)
		      {
			//The weight is reversed... I think...
			rank_weight1=1/(double)(rank[j]);
			rank_weight2=1/(double)(rank[i]);
			//  printf("%f %d %f %d\n",rank_weight1,rank[j],rank_weight2,rank[i]);
		      }
		    if(rank1) {
		      
		      rank_weight1=0;
		      rank_weight2=0;
		      if(rank[j]==1) {
			rank_weight1=1;
		      }
		      if(rank[i]==1) {
			rank_weight2=1;
		      }
		    }
		    
		    if(rank_weight1 > 0 || rank_weight2 > 0) {
		      if(memorymode)
			{
			  LGscore_res_pt(&dm[i],&dm[j],LG,d0,minsim,L,factor,step); //,ignore_res);
			}
		      else
			{
			  LGscore_res(filenames_with_path[i],filenames_with_path[j],LG,d0,minsim,L,factor,step);
			}
		    
		    //		  printf("%lf\n",LG[0].LGscore);
		    //if(verbose)
		    //  printf("%d Comparing %s %s %s %s %lf %lf\n",total_number_of_comparisons,filenames[i],filenames[j],dm[i].filename,dm[j].filename,LG[0].LGscore,LG[0].Ssum);
		      if(output_similarity_matrix)
			fprintf(fp,"%6.3lf ",LG[0].LGscore);
		    
		    
		      sim_mat[i][j]=LG[0].LGscore;
		      sim_mat[j][i]=LG[0].LGscore;
#pragma omp atomic
		      LG_average[i]+=rank_weight1*LG[0].LGscore;
#pragma omp atomic

		      LG_average[j]+=rank_weight2*LG[0].LGscore;
#pragma omp atomic
		      S_average[i]+=rank_weight1*LG[0].Ssum;
#pragma omp atomic
		      S_average[j]+=rank_weight2*LG[0].Ssum;
#pragma omp atomic
		      number_of_comparisons[i]+=rank_weight1;
#pragma omp atomic

		      number_of_comparisons[j]+=rank_weight2;
		      
		      if(rank[i]==1)
			{
			//  printf("%d %s %s %8.3f %e\n",i,filenames[i],filenames[j],LG[0].LGscore,pow(10,-LG[0].LGscore));
#pragma omp atomic
			  LG_average1[j]+=rank_weight2*LG[0].LGscore;
#pragma omp atomic
			  S_average1[j]+=rank_weight2*LG[0].Ssum;
#pragma omp atomic
			  number_of_comparisons1[j]+=rank_weight2;
			}
		      
		      if(rank[j]==1)
			{
			  //  printf("%d %s %s\n",j,filenames[j],filenames[i]);
			  //printf("%d %s %s %8.3f %e \n",j,filenames[j],filenames[i],LG[0].LGscore,pow(10,-LG[0].LGscore));
#pragma omp atomic
			  LG_average1[i]+=rank_weight1*LG[0].LGscore;
#pragma omp atomic
			  S_average1[i]+=rank_weight1*LG[0].Ssum;
#pragma omp atomic
			  number_of_comparisons1[i]+=rank_weight1;
			}
		      
#pragma omp atomic
		      total_number_of_comparisons++;
		      if(verbose)
			printf("%d Comparing2 %s %s %s %s %lf %lf %lf %lf %lf %lf\n",total_number_of_comparisons,filenames[i],filenames[j],dm[i].filename,dm[j].filename,LG[0].LGscore,LG[0].Ssum,S_average1[i],S_average1[j],rank_weight1,rank_weight2);
		    
		    //if(strcmp(filenames[i],"
		    
		    //if(LG[0].residues>maxlen)
		    //  {
		    //    maxlen=LG[0].residues;
		    //  }
		    // In the case of CASP I can assume that the resnum is the correct index
		      //#pragma omp parallel for firstprivate(rank_weight1, rank_weight2) schedule(auto)
		      for(int k=0;k<LG[0].residues;k++)
			{
			  //   printf("%d %d\n",k+1,LG[0].resnum[k]);
			  //  if(i==97 || j==97 ) { 
			  //	if(LG[0].resnum[k]==145 && LG[0].S[k]>0) 
			  //	  {
			  //	    printf("%d i %d: %s %d %lf\n",k,i, filenames[i],LG[0].resnum[k],LG[0].S[k]);
			  //	    printf("%d j %d: %s %d %lf\n",k,j, filenames[j],LG[0].resnum[k],LG[0].S[k]);
			  //	    for(l=0;l<LG[0].residues;l++)
			  //	      {
			  //		printf("%d %d %lf\n",l+1,LG[0].resnum[l],LG[0].S[l]);
			  //	      }
			  //	    exit(1);
			  //	  }
			  //}
			  
			  //	printf("%d %s %s %lf\n",k,filenames[i],filenames[j],rank_weight1);
#pragma omp atomic
			  Sstr[i][LG[0].resnum[k]-1]+=rank_weight1*LG[0].S[k];
#pragma omp atomic
			  Sstr[j][LG[0].resnum[k]-1]+=rank_weight2*LG[0].S[k];
			}
		    }
		    else
		      {
			if(verbose)
			  {
			    printf("Does not compare %s %s SameMethod: %d\n",filenames[i],filenames[j],SameMethod(filenames[i],filenames[j]));
			  }
		      }

		  }
	      }
	    }
	  if(output_similarity_matrix)
	    fprintf(fp,"\n");
	  //	  free(LG);
	}
#pragma omp barrier

#ifdef _OPENMP
       double ftime=omp_get_wtime()-stime;
#else 
       double ftime=(double)clock() / (double)CLOCKS_PER_SEC-stime;
#endif


       //}
      if(output_similarity_matrix)
	{
	  fprintf(fp,"END\n");
	  printf("Done writing pairwise similarity...\n");
	  if(only_pairwise)
	    exit(1);
	}
      //sqrt(1/${$average_similarity[$i]}[$j]-1)*$sqrt5;
      
      en_time=times(&en_cpu);

     
      //Calculate initial Pcons score
      // for(i=0;i<maxlen;i++) {
      //	printf("SAM: %s %d %lf\n",filenames[6],i+1,Sstr[6][i]);
      //  }
      for(i=0;i<files;i++)
	{
	  if(number_of_comparisons[i]!=0) {
	    LG_average[i]/=number_of_comparisons[i];
	    S_average[i]/=number_of_comparisons[i];
	  }  
	  //else {
	  // LG_average[i]=0;
	  //  S_average[i]=0;
	  //} 
	  if(number_of_comparisons1[i]!=0) {
	    LG_average1[i]/=number_of_comparisons1[i];
	    S_average1[i]/=number_of_comparisons1[i];
	  } 
	  //else {
	  // LG_average1[i]=0;
	  //  S_average1[i]=0;
	  // }
	  if(normL)
	    {
	      if(userdef_len>0) {
		S_average[i]/=userdef_len;
		S_average1[i]/=userdef_len;
	      } else {
		S_average[i]/=maxlen;
		S_average1[i]/=maxlen;
	      }
	    }
	  //	  printf("%lf %lf %d\n",S_average1[i],number_of_comparisons1[i],maxlen);
	}


      //Calculate reweighted pcons scores...
      if(google){
	  lgscoreoutput=0;
	  google_weight=malloc(sizeof(double)*files);
	  pcons_google=malloc(sizeof(double)*files);
	  pcons_google_previous=malloc(sizeof(double)*files);
	  for(i=0;i<files;i++) {
	      pcons_google_previous[i]=LG_average[i];
	    }

	  for(k=0;k<google_iter;k++){
	    for(i=0;i<files;i++){
	      pcons_google[i]=0;
	      google_weight[i]=pcons_google_previous[i];
	    }
	      //	   temp_counter=0;
	      for(i=0;i<files;i++){
		  temp_counter=0;
		  for(j=0;j<files;j++){
		    if(sim_mat[i][j]!=-1){
		      if(google_weight[j] > google_weight_cut)
			{
			  pcons_google[i]+=sim_mat[i][j]*google_weight[j];
			  temp_counter++;
			}
		    }
		  }
		  //	printf("%d %d %6.3f %6.3f %6.3f\n",temp_counter,number_of_comparisons[i],pcons_google[i],pcons_google[i]/number_of_comparisons[i],LG_average[i]);
		  if(temp_counter>0)
		    pcons_google[i]/=temp_counter;
		  pcons_google_previous[i]=pcons_google[i];
	      }
	  }
	  // exit(1);
	  free(google_weight);
	  free(pcons_google);
	  free(pcons_google_previous);
	}


      if(files>1)
	{
	  if(userdef_len >= 0 && userdef_len<maxlen)
	    {
	      fprintf(stderr,"WARNING the specified target length is shorter than the longest model (%d vs %d)\n",userdef_len,maxlen);
	      fprintf(stderr,"By removing the -L switch WARNING the specified target length is shorter than the longest model (%d vs %d)\n",userdef_len,maxlen);
	    }


	  printf("PFRMAT QA\n");
	  printf("TARGET %s\n",target);
	  printf("AUTHOR 0566-5283-0576\n"); //1645-7339-3688\n");
	  if(strcmp(ignore_res_file,"undef")!=0) {

	    printf("REMARK Residues in the file: \'%s\' were ignored in the calculation\n",ignore_res_file);
	  }
	  if(caspoutput)
	    {
	      printf("REMARK Error estimate is CA-CA rms distance in Angstroms\n");
	    }
	  else
	    {
	      printf("REMARK Error estimate is a quality measure between 0 and 1 local S-score 1/(1+rms^2/5)\n");
	    }
	  if(lgscoreoutput)
	    {
	      printf("REMARK Global quality estimate is LGscore\n");
	    }
	  else{
	    if(google)
	      {
		printf("REMARK Global quality estimate is reweighted LGscore (google-style, iter=%d, cut=%5.3lf)\n",google_iter,google_weight_cut);
	      }
	    else
	      {
		if(userdef_len > 0) {
		  printf("REMARK Global quality estimate is average S-score normalised by len=%d\n",userdef_len);
		} else {
		  printf("REMARK Global quality estimate is average S-score normalised by len=%d\n",maxlen);
		}
	      }
	  }
	  if(strlen(dir)>0)
	    {
	      printf("REMARK Input directory: %s\n",dir);
	    }
	  if(strlen(infilename)>0)
	    {
	      printf("REMARK Input file: %s\n",infilename);
	    }
	  printf("REMARK Total number of pdbs: %d\n",files);
	  printf("REMARK Total number of comparisons: %d\n",total_number_of_comparisons);
	  printf("REMARK Total CPU time used: %8.2lf seconds\n",(double)(((double)en_cpu.tms_utime - (double)st_cpu.tms_utime)/100));
	  printf("REMARK Total time used: %lf seconds\n",ftime);
	  printf("REMARK LGscore parameters superimpose_all=%d L=%d minsim=%lf factor=%lf step=%d\n",superimpose_all,L,minsim,factor,step);
	  if(userdef_len>0)
	    {
	      printf("REMARK Target sequence length (user defined): %d\n",userdef_len);
	      if(userdef_len<maxlen)
		{
		  printf("REMARK WARNING the specified target length is shorter than the longest model (%d vs %d)\n",userdef_len,maxlen);
		}
	      
	    }
	  else
	    {
	      printf("REMARK Assumed target sequence length (max resnum present in any of the pdbs): %d\n",maxlen);
	    }
	  
	  //if(fastmode)
	  //  {
	  //    printf("REMARK Method run in fast mode\n");
	  //  }
	  //else
	  //  {
	  //    printf("REMARK Method run in slow mode\n");
	  //  }
	  //if(memorymode)
	  //  {
	  //    printf("REMARK Method run in memory mode\n");
	  //  }
 
	  printf("METHOD ********************************************************************************\n");
	  printf("METHOD *                               PCONS                                          *\n");
	  printf("METHOD * Calculates the model quality based on structural consensus                   *\n");
	  printf("METHOD *                                                                              *\n");
	  printf("METHOD * Reference: Bjorn Wallner and Arne Elofsson, Protein Sci., 2006 15(4):900-913 *\n");
	  printf("METHOD * For comments, please email to: bjorn@sbc.su.se                               *\n");
	  printf("METHOD *                                                                              *\n");
	  printf("METHOD *                                                                              *\n");
	  if(lgscoreoutput)
	    {
	      printf("METHOD * Global: Average LGscore                                                      *\n");
	    }
	  else
	    {
	      if(normL)
		{
		  printf("METHOD * Global: Average pairwise S-score divided by target length                    *\n");
		}
	      else
		{
		  printf("METHOD * Global: Average pairwise S-score sum                                         *\n");
		}
	    }
	  if(rank1) 
	    {
	      
	      printf("METHOD * Global: Only comparison to rank1 models                                      *\n");
	      printf("METHOD * Local:  Only comparison to rank1 models                                      *\n");
	    }  
	  else if(use_rank_weight)
	    {
	      printf("METHOD * Global: Average weighted by 1/rank                                           *\n");
	    }
	  
	  if(w1>0)
	    {
	      printf("METHOD * Global: Comparison to first ranked models included with %5.2lf:%5.2lf          *\n",w1,(1-w1));       
	    }
	   if(caspoutput)
	    {
	      printf("REMARK Error estimate is CA-CA distance in Angstroms                                  \n");
	    }
	  else
	    {
	      printf("METHOD * Local:  (Pcons-local) The average pairwise S_i=1/(1+(rms_i/5)^2) all against *\n");
	    }
	   printf("METHOD ********************************************************************************\n");
	  printf("MODEL 1\n");
	  printf("QMODE 2\n");
	  

	  //if(userdef_len>0)
	    //{
	      //maxlen=userdef_len;
	      //  }
	  for(i=0;i<files;i++)
	    {
	      //LG_average[i]/=number_of_comparisons[i];
	      //S_average[i]/=number_of_comparisons[i];
	      //S_average[i]/=maxlen;

	      if(lgscoreoutput)
		{
		  //printf("%s %5.3lf %5.3lf %5.3lf",filenames[i],(1-w1)*LG_average[i]+w1*LG_average1[i],LG_average[i],LG_average1[i]);
		  printf("%s %5.3lf ",filenames[i],(1-w1)*LG_average[i]+w1*LG_average1[i]);
		}
	      else{
		if(google)
		  {
		    printf("%s %5.3lf ",filenames[i],pcons_google[i]);
		  }
		else
		  {
		    // printf("%s %5.3lf %5.3lf %5.3lf\n",filenames[i],(1-w1)*S_average[i]+w1*S_average1[i],S_average[i],S_average1[i]);
		    printf("%s %5.3lf ",filenames[i],(1-w1)*S_average[i]+w1*S_average1[i]);
		    
		  }
	      }
	      //   printf("%s %5.3lf %5.3lf ",filenames[i],LG_average[i],S_average[i]);
	      for(j=0;j<maxlen;j++)
		{
		  Sstr[i][j]/=number_of_comparisons[i];
		  if(Sstr[i][j]>0.0001)
		    {
		      //Sstr[i][j]/=number_of_comparisons[i];
		      // Sstr[i][j]/=number_of_comparisons[i];
		      if(caspoutput==1)
			{
			  rmsd=d0*sqrt(1/Sstr[i][j]-1);
			  if(rmsd>15)
			    {
			      rmsd=15;
			    }
		      
			  printf(" %4.3lf",rmsd);
			}
		      else
			{
			  printf(" %4.3lf",Sstr[i][j]);
			}
		    }
		  else
		    {
		      printf(" X");
		    }
		}
	      for(;j<userdef_len;j++) {
		printf(" X");
	      }
	      printf("\n");
	      
	    }
	  printf("END\n");
	}
      else
	{
	  printf("Too few pdbfiles (%d) found\n",files); 
	}

      //Freeing up some allocated memory
      if(nofree) {
	exit(1);
      }
      bytes_freed=0;
      for(i=0;i<files;i++)
	 {
	   // printf("%d\n",i); 
	   free_dyn_molecule(&dm[i]);
	  //printf("%lf \n",Sstr[i][0]);
	   free(Sstr[i]);
	   free(sim_mat[i]);
	   free(filenames[i]);
	   free(filenames_with_path[i]);
	 }
// free(ignore_res);
      free(Sstr);
      free(sim_mat);
      free(LG_average);
      free(S_average);
      free(LG_average1);
      free(S_average1);
      free(number_of_comparisons);
      free(number_of_comparisons1);
      free(dm);

    
    }
  else
    {
      usage();
    
    }
}


void usage()
{
  
  fprintf(stderr,"********************************************************************************\n");
  fprintf(stderr,"*                               PCONS                                          *\n");
  fprintf(stderr,"* Calculates structural consensus for all models in a specified directory      *\n");
  fprintf(stderr,"*                                                                              *\n");
  fprintf(stderr,"* Reference: Bjorn Wallner and Arne Elofsson, Protein Sci., 2006 15(4):900-913 *\n");
  fprintf(stderr,"* For comments, please email to: bjorn@sbc.su.se                               *\n");
  fprintf(stderr,"*                                                                              *\n");
  fprintf(stderr,"* If the program Seg. faults make sure the stacksize is unlimited:             *\n");
  fprintf(stderr,"* bash: ulimit -s unlimited                                                    *\n");
  fprintf(stderr,"* tcsh: limit stacksize unlimited                                              *\n");
  fprintf(stderr,"*                                                                              *\n");
  fprintf(stderr,"********************************************************************************\n");
  fprintf(stderr,"\nUsage: pcons\t-d <directory w/ models>\n");
  fprintf(stderr,"\t\t-i <inputfile containing full path to models (one one each line) OBS the filenames w/o path must by unique>\n");
  fprintf(stderr,"\t\t-casp will output local rmsd as local quality (default is average local S-score)\n");
  fprintf(stderr,"\t\t-lgscore will output average LGscore as global quality measure (default is average S-score)\n");
  fprintf(stderr,"\t\t-L <target sequence length. default: longest sequence in the set>\n");
  fprintf(stderr,"\t\t-t <target id, default: will look for T0XXX in the beginning of filenames >\n");
  fprintf(stderr,"\t\t-m <output pairwise similarity matrix to this file, this will force compare all targets against each other>\n");
  fprintf(stderr,"\t\t-A <force compare all targets, by default targets from same method are not compared>\n");
  fprintf(stderr,"\t\t-only_pairwise <will skip the pcons evaluation>\n");
  fprintf(stderr,"\t\t-google  <do a google-like weighting>\n");
  fprintf(stderr,"\t\t-w1 <the final score is (1-w1)*sim_all+w1*sim_first>\n");
  fprintf(stderr,"\t\t-use_rank_weight <weight the average sim with 1/rank>\n");
  fprintf(stderr,"\t\t-rank1 <only compare to first ranked, this will put w1=1 as well and override any setting of w1>\n");
  fprintf(stderr,"\t\t-normL <normalise S-score by length>\n");
  fprintf(stderr,"\t\t-ignore_res <file with residues to ignore (one on each line)>\n");
  fprintf(stderr,"\t\t-superimpose_all <use only one superposition>\n");
 //fprintf(stderr,"\t\t-f <fast mode, more sloppy heuristics for finding optimal sub alignments>\n");
  //fprintf(stderr,"\t\t-m <memory options\n");
  //fprintf(stderr,"\t\t    0 = re-read from disk\n");
  //fprintf(stderr,"\t\t    1 = read all coordinates in to memory (default)>\n");

  exit(0);

}




int SameMethod(char *file1, char *file2)
{
  int ident=0;
  int diffpos=0;
  int i=0,j=0;
  int len1,len2;
  int dot_counter=0;
  char method1[100];
  char method2[100];

  len1=strlen(file1);
  len2=strlen(file2);




  if(file1[len1-3]=='T' && file1[len1-2]=='S' && file2[len2-3]=='T' && file2[len2-2]=='S' ||
     file1[len1-7]=='A' && file1[len1-6]=='L' & file1[len1-3]=='p' && file1[len1-2]=='d' && file1[len1-1]=='b' && file2[len2-7]=='A' && file2[len2-6]=='L' & file2[len2-3]=='p' && file2[len2-2]=='d' && file2[len2-1]=='b' )
    {
      //CASP naming convension ending with TS{rank}
      //OR AL{rank}.pdb 
      if(len1 != len2)
	{
	  return 0;
	}
     
      
      //All files here will have the same length
      //Check for identical matches
      
      for(i=0;i<len1;i++)
	{
	  if(file1[i] == file2[i])
	    {
	      ident++;
	    }
	  else
	    {
	      diffpos=i+1;
	    }
	}
      
      // printf("id: %d len %d pos %d ext %s\n",ident,len1,diffpos,&file1[len1-3]);
      
      if(len1-ident == 1 && (diffpos==len1 || diffpos==len1-4 && strcmp(&file1[len1-2],"pdb")==0))
	{
	  //      printf("%s same as %s\n",file1,file2);
	  
	  return 1;
	  
	}
      else
	{
	  //printf("%s not same as %s\n",file1,file2);
	  return 0;
	}
    }
   if(file1[len1-3]=='T' && file1[len1-2]=='S' && file2[len2-3]=='p' && file2[len2-2]=='d' && file2[len2-1]=='b'||
      file2[len2-3]=='T' && file2[len2-2]=='S' && file1[len1-3]=='p' && file1[len1-2]=='d' && file1[len1-1]=='b')
     {
       return 0;
     }
   else  //assume the following format target.template.method.rank.pdb 
     {
       dot_counter=0;
       j=0;
       for(i=0;i<len1;i++){
	 if(file1[i]=='.'){
	   dot_counter++;
	 }else{
	   if(dot_counter==2){
	     method1[j]=file1[i];
	     j++;
	   }
	 }
       }
       method1[j]='\0';
       j=0;
       dot_counter=0;
       for(i=0;i<len2;i++){
	 if(file2[i]=='.'){
	   dot_counter++;
	 }else{
	   if(dot_counter==2){
	     method2[j]=file2[i];
	     j++;
	   }
	 }
       }
       method2[j]='\0';
       
       
       //rintf("IN SAME: %s %s %d\n",method1,method2,strcmp(method1,method2));
       if(strcmp(method1,method2)==0)
	 {
	   return 1;
	 }
       else
	 {
	   return 0;
	 }
     }
   
}
int Rank(char *file)
{
  int i=0,j=0;
  int len1,len2;
  int dot_counter=0;
  int rank=0;
  char temp[2];
  len1=strlen(file);
 

  //  printf("R: %s\n",file);
  if(file[len1-3]=='T' && file[len1-2]=='S')
    {
      //  printf("TS ");
      temp[0]=file[len1-1];
      temp[1]='\0';
      rank=atoi(temp);
    }
  else 
    {
      if(file[len1-7]=='A' && file[len1-6]=='L' && file[len1-3]=='p' && file[len1-2]=='d' && file[len1-1]=='b')
	{
	  //  printf("AL ");
	  temp[0]=file[len1-5];
	  temp[1]='\0';
	  rank=atoi(temp);
	}
      else  //assume the following format target.template.method.rank.pdb 
	{
	  dot_counter=0;
	  j=0;
       for(i=0;i<len1;i++){
	 if(file[i]=='.'){
	   dot_counter++;
	 }else{
	   if(dot_counter==3){
	     temp[0]=file[i];
	     temp[1]='\0';
	     rank=atoi(temp);
	     j++;
	   }
	 }
       }
       if(j==2)
	 rank=10;
	}
    }
  return rank;
}

int ispdb(char *filename)
{
  molecule m[1];
  m[0].residues=0;
  strcpy(m[0].filename,"undef");
  strcpy(m[0].filename,filename);
  //printf("ispdb: %s\n",filename);
  if(read_molecules(&m[0],'c')==0)
    {
      if(m[0].residues>0)
	{
	  return 1;
	}
    }

  return 0;
}


char *basename(char *fname)
{
  int sl;
    char *s;
    // if (strlen(fname) > STRING_BUFFER-4)
    //  error(NIL,"filename is too long",NIL);
    s = strrchr(fname,'/');
    if (s)
      fname = s+1;
    /* Process suffix */
    return fname;

}

