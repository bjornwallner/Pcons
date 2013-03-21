#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <sys/types.h>
#include <dirent.h>
#include <sys/times.h>
#include <time.h>
#include <limits.h>
#include "molecule.h"
#include "lgscore.h"
//#include "/afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/quality_molecule.h"
#define MAXFILE 2000
#define STRING_BUFFER 1000
void usage();
int SameMethod(char *file1, char *file2);
int ispdb(char *filename);
char *basename(char *fname);

main(int argc,char *argv[])             /* Main routine */
{
  DIR            *dip;
  FILE           *fp;
  struct dirent  *dit;
  char           dir[PATH_MAX]="";
  char           file2[PATH_MAX];
  char           tempfilename[PATH_MAX];
  char           infilename[PATH_MAX]="";
  char           filenames[MAXFILE][100];
  char           filenames_with_path[MAXFILE][PATH_MAX];
  //char           **filenames='\0'; //[MAXFILE][100];
  //char           **filenames_with_path='\0';//[MAXFILE][STRING_BUFFER];
  int            i,j,k;
  int            files=0;
  int            L=4;
  int            step=2;
  double         minsim=121.0;
  double         d0=sqrt(5);
  double         factor=0.5;
  char           temp[PATH_MAX];
  lgscore        LG[1];
  //  double         LG_average[MAXFILE]={0};
  //double         S_average[MAXFILE]={0};
  //double         Sstr[MAXFILE][MAXRES]={{0}}; 
  //  int            number_of_comparisons[MAXFILE]={0};
  double         *LG_average,*S_average;
  double          **Sstr;
  int            *number_of_comparisons;
  int            total_number_of_comparisons=0;
  int            maxlen=0;
  int            userdef_len=-1;
  double         rmsd=0;
  double        sqrt5=sqrt(5);
  char          target[STRING_BUFFER]="T0XXX";
  int           fastmode=0;
  int           memorymode=1;
  int           bytes_freed=0;
  int           target_defined_by_user=0;
  int           caspoutput=0;
  int           lgscoreoutput=0;
  int           maxfile=MAXFILE;
  //molecule      m[MAXFILE];
  dyn_molecule  *dm;
  
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
	  if (strcmp(argv[i],"-d")==0)
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
			  //printf("%s\n",tempfilename);
			  //filenames=(char**)realloc(filenames,sizeof(char)*20); //(files+1));
			  ////filenames[files]=(char*)realloc(filenames[files],sizeof(char)*(strlen(dit->d_name)+1));
			  //filenames[files]=(char*)realloc(filenames[files],sizeof(char)*20000);
			  //printf("%s\n","first");
			  //
			  //filenames_with_path=(char**)realloc(filenames_with_path,sizeof(char)*20); //(files+1));
			  //// filenames_with_path[files]=(char*)realloc(filenames_with_path[files],sizeof(char)*(strlen(tempfilename)+1));
			  //filenames_with_path[files]=(char*)realloc(filenames_with_path[files],sizeof(char)*2000);
			  //printf("%s\n","second");
			  
			  strcpy(filenames[files],dit->d_name);
			  strcpy(filenames_with_path[files],tempfilename);
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
			}
		    }
		}   
	   

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
		      //		      printf("%s %d\n",tempfilename,strlen(tempfilename));
		      if(ispdb(tempfilename))
			{
			  strcpy(filenames_with_path[files],tempfilename);
			  //printf("%s %d\n",filenames_with_path[files],strlen(filenames_with_path[files]));
			  ///			  printf("%s %d\n",filenames_with_path[files],strlen(filenames_with_path[files]));
			  
			  strcpy(filenames[files],basename(tempfilename));
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
		}
	      else
		{
		  printf("Cannot open %s",tempfilename);
		  exit(0);
		}
	      fclose(fp);
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
	   if(strcmp(argv[i],"-casp")==0)
	    {
	      caspoutput=1;
	    }
	   if(strcmp(argv[i],"-lgscore")==0)
	    {
	      lgscoreoutput=1;
	    }
	  i++; 
	}
    
	      
      st_time=times(&st_cpu);
      // Read all CA coordinates to memory
      if(memorymode)
	{
	  // printf("%d",files);
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
	      // printf("%s\n",dm[i].filename);
	      if(read_molecules_dynamic(&dm[i],'c')==0)
	      	{
	      	  if(dm[i].atm[dm[i].residues-1].resnum>maxlen)
		    {
		      maxlen=dm[i].atm[dm[i].residues-1].resnum;
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

     

      LG_average=malloc(sizeof(double)*files);
      S_average=malloc(sizeof(double)*files);
      number_of_comparisons=malloc(sizeof(int)*files);
      Sstr=(double**)malloc(sizeof(double*)*files);
      // printf("%d %d\n",sizeof(double*),sizeof(double));
      //exit(1);
      
      if(userdef_len>maxlen)
	{
	  maxlen=userdef_len;
	}

      //  printf("%d\n",maxlen);
      //exit(1);
      for(i=0;i<files;i++)
	{
	  LG_average[i]=0;
	  S_average[i]=0;
	  number_of_comparisons[i]=0;
	  Sstr[i]=(double*)malloc(sizeof(double)*maxlen);
	  //printf("%d\n",i);
	  for(j=0;j<maxlen;j++)
	    {
	      //printf("%d %d\n",i,j);
	      Sstr[i][j]=0;
	    }
	}
      
      //      printf("#Maximum number of comparisons: %d (probably lower if there are models from same method)\n",(int)files*(files-1)/2);
      for(i=0;i<files;i++)
	{
	  //printf("%d\n",total_number_of_comparisons);
	  for(j=i+1;j<files;j++)
	    {
	      //printf("%s %s %s %s %d %d same:%d\n",dm[i].method,dm[j].method,filenames[i],filenames[j],i,j,!SameMethod(filenames[i],filenames[j]));
	      
	      if((dm[i].method==NULL || dm[j].method==NULL) && !SameMethod(filenames[i],filenames[j]) ||
		 (dm[i].method!=NULL && dm[j].method!=NULL) && strcmp(dm[i].method,dm[j].method)!=0)
		{
		  //printf("IN!\n");
		  //printf("%s %s %d %d\n",dm[i].method,dm[j].method,strcmp(dm[i].method,dm[j].method),SameMethod(filenames[i],filenames[j]));
		  // printf("%d Comparing %s %s %s %s\n",total_number_of_comparisons,filenames[i],filenames[j],dm[i].filename,dm[j].filename);
		  if(memorymode)
		    {
		      LGscore_res_pt(&dm[i],&dm[j],LG,d0,minsim,L,factor,step);
		    }
		  else
		    {
		      LGscore_res(filenames_with_path[i],filenames_with_path[j],LG,d0,minsim,L,factor,step);
		    }

		  //		  printf("%lf\n",LG[0].LGscore);
		  
		  LG_average[i]+=LG[0].LGscore;
		  LG_average[j]+=LG[0].LGscore;
		  S_average[i]+=LG[0].Ssum;
		  S_average[j]+=LG[0].Ssum;

		  number_of_comparisons[i]++;
		  number_of_comparisons[j]++;

		  


		  total_number_of_comparisons++;


		  //if(strcmp(filenames[i],"

		  //if(LG[0].residues>maxlen)
		  //  {
		  //    maxlen=LG[0].residues;
		  //  }
		  // In the case of CASP I can assume that the resnum is the correct index
		  for(k=0;k<LG[0].residues;k++)
		    {
		      // printf("%d %lf\n",LG[0].resnum[k],LG[0].S[k]);
		      Sstr[i][LG[0].resnum[k]-1]+=LG[0].S[k];
		      Sstr[j][LG[0].resnum[k]-1]+=LG[0].S[k];
		    }
		}
	//else
	//	{
	//	  printf("Does not compare %s %s\n",filenames[i],filenames[j]);
	//	}
	    }
	}
      //sqrt(1/${$average_similarity[$i]}[$j]-1)*$sqrt5;
      
      en_time=times(&en_cpu);

     
      
      




      if(files>1)
	{
	  if(userdef_len >= 0 && userdef_len<maxlen)
	    {
	      fprintf(stderr,"WARNING the specified target length is shorter than the longest model (%d vs %d)\n",userdef_len,maxlen);
	      fprintf(stderr,"By removing the -L switch WARNING the specified target length is shorter than the longest model (%d vs %d)\n",userdef_len,maxlen);
	    }


	  printf("PFRMAT QA\n");
	  printf("TARGET %s\n",target);
	  //	  printf("AUTHOR 1645-7339-3688\n");
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
	  else
	    {
	      printf("REMARK Global quality estimate is average S-score\n");
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
	  printf("METHOD * Calculates the model quality based on structral consensus                    *\n");
	  printf("METHOD *                                                                              *\n");
	  printf("METHOD * Reference: Björn Wallner and Arne Elofsson, Protein Sci., 2006 15(4):900-913 *\n");
	  printf("METHOD * For comments, please email to: bjorn@sbc.su.se                               *\n");
	  printf("METHOD *                                                                              *\n");
	  printf("METHOD *                                                                              *\n");
	  printf("METHOD * Global: Average pairwise S-score divided by target length                    *\n");
	  printf("METHOD * Local:  (Pcons-local) The average pairwise S_i=1/(1+(rms_i^2/5) all against  *\n");
	  printf("METHOD *         all then rescaled back to rms.                                       *\n");
	  printf("METHOD ********************************************************************************\n");
	  printf("MODEL 1\n");
	  printf("QMODE 2\n");
	  

	  if(userdef_len>0)
	    {
	      maxlen=userdef_len;
	    }
	  for(i=0;i<files;i++)
	    {
	      LG_average[i]/=number_of_comparisons[i];
	      S_average[i]/=number_of_comparisons[i];
	      S_average[i]/=maxlen;

	      if(lgscoreoutput)
		{
		  printf("%s %5.3lf ",filenames[i],LG_average[i]);
		}
	      else
		{
		  printf("%s %5.3lf ",filenames[i],S_average[i]);
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
	      printf("\n");
	      
	    }
	  printf("END\n");
	}
      else
	{
	  printf("Too few pdbfiles (%d) found\n",files); 
	}

      //Freeing up some allocated memory
     
      bytes_freed=0;
      for(i=0;i<files;i++)
	 {
	   //printf("%d\n",i); 
	  free_dyn_molecule(&dm[i]);
	  //printf("%lf \n",Sstr[i][0]);
	  free(Sstr[i]);
	 }
      
      free(Sstr);
      free(LG_average);
      free(S_average);
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
  fprintf(stderr,"* Calculates structral consensus for all models in a specified directory       *\n");
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
  fprintf(stderr,"\t\t-casp <output local rmsd as local quality (default is average local S-score)>\n");
  fprintf(stderr,"\t\t-lgscore <output average LGscore as global quality measure (default is average S-score)>\n");
  fprintf(stderr,"\t\t-L <target sequence length. default: longest sequence in the set>\n");
  fprintf(stderr,"\t\t-t <target id, default: will look for T0XXX in the beginning of filenames >\n");
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




  if(file1[len1-3]=='T' && file1[len1-2]=='S' && file2[len2-3]=='T' && file2[len2-2]=='S')
    {
      //CASP naming convension ending with TS{rank}
      
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
      
      if(len1-ident == 1 && (diffpos==len1 || diffpos==len1-4 && strcmp(&file1[len1-2],"pdb")))
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
  else
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
      

      // printf("%s %s %d\n",method1,method2,strcmp(method1,method2));
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


int ispdb(char *filename)
{
  molecule m[1];
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

