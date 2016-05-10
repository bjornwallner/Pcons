#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "molecule.h"
//#include "/afs/pdc.kth.se/home/b/bjornw/bjorn/modules/numrecipies/nr.h"
//#include "/afs/pdc.kth.se/home/b/bjornw/bjorn/modules/numrecipies/nrutil.h"
//#include "/afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/nets/nets.h"
#include "nr.h"
#include "nrutil.h"
#include "nets.h"


/* changed input to (molecule *m,char atomflag)
   atomflag == a -> read all atoms (except H)
   atomflag == c -> CA atoms
   atomflag == b -> backbone CA,C,N,O atoms
*/

int read_molecules_dynamic(dyn_molecule *m,char atomflag,int *ignore_res)	/* Reads in molecules to be superimposed */
{
  int	i,j,k,atoms,residues;	/* Counter variables */
  char	buff[3000];	/* Input string */
  char	temp[3000];	/* Throwaway string */
  char	line_flag[3000];	/* PDB file line mode */
  char  desc_flag[3000];
  char	residue[4];	/* PDB atom info for output */
  char	name[4];
  int   resnum;
  char  resname[9];
  char  old_resname[9]="noname";
  char alt_loc_check[2]=" ";
  char x_temp[11];
  char y_temp[11];
  char z_temp[11];
  //char	number[8];
  int number;
  char chain[2];
  char inscode[2]=" ";
  char alt_loc[2]=" ";
  double	x,y,z;		/* Temporary coordinates values */
  double bfactor=9.99;
  FILE *fp;
  char temp_number[11];
  char temp_resnum[11];
  int ss_flag=0;

  i=0; /* It is only one molecule to be read this was done for the moment instead of changing all "i" to 0 */
  

  //for (i=0;i<2;i++)
  //{
  //#ifdef  NOZLIB*/
  //fp=fopen(m[i].filename,"r");	/* Does file exist? */
  //#else
  //fp=gzopen(m[i].filename,"r");	/* Does file exist? */
  //#endif
  //printf("%s\n",m[0].filename);
  //strcpy(m[0].method,"undef");
  m[0].rank=-1;
  m[0].score=-9999;

  fp=fopen(m[0].filename,"r");	/* Does file exist? */
  m[0].sequence='\0';
  m[0].ss='\0';
  m[0].method='\0';
  m[0].atm='\0';
  m[0].CA_ref='\0';
  m[0].residues=0;
  m[0].atoms=0;
  
  if (fp!=NULL)	/* If yes, read in coordinates */
    {
      /* Initialize things */
      //m[0].xcen=m[0].ycen=m[0].zcen=0;
      atoms=0;
      residues=0;
      //#ifdef NOZLIB
      //while(fgets(buff,255,fp)!=NULL)
      //#else
      //while(gzgets(fp,buff,255)!=NULL)
      //#endif
      while(fgets(buff,1000,fp)!=NULL)
	{
	  //sscanf(buff,"%-6s %4d  %-3s%1s%3s %1s%5s    %7.3lf %7.3lf %7.3lf",line_flag,&number,name,residue,chain,&resnum,inscode,&x,&y,&z);
	  //sscanf(buff,"%s %d %s %s %s",line_flag,&number,name,residue,resname);
	  strcpy(line_flag,"undef");
	  strcpy(desc_flag,"undef");
	  strcpy(temp,"undef");
	  sscanf(buff,"%s",line_flag);
	  // printf("%s %s %s",m[0].filename,line_flag,&buff);
	  if(strcmp("REMARK",line_flag)==0)
	    {
	      sscanf(buff,"%s %s %s",line_flag,desc_flag,temp);
	      if(strcmp("SS",desc_flag)==0)
		{
		  m[0].ss=malloc(sizeof(char)*(strlen(temp)+1));
		  strcpy(m[0].ss,temp);
		  ss_flag=1;
		  //printf("%s\n",m[0].ss);
		}
	      if(strcmp("METHOD",desc_flag)==0)
		{
		  m[0].method=malloc(sizeof(char)*(strlen(temp)+1));
		  strcpy(m[0].method,temp);
		  //printf("%s\n",m[0].method);
		}
	      if(strcmp("SCORE",desc_flag)==0)
		{
		  m[0].score=atof(temp);
		  //strcpy(m[0].method,temp);
		  //printf("%lf\n",m[0].score);
		}
    	    }
	  if(strcmp("MODEL",line_flag)==0)
	    {
	      sscanf(buff,"%s %s",line_flag,temp);
	      m[0].rank=atoi(temp);
	    }
	  
	  if(strcmp("ATOM",line_flag)==0 && (buff[12] != 'H' || buff[13] != 'H' ))	/* Is it an ATOM entry? */
	    {
	      //printf("Heja: %s %s",m[0].filename,&buff);
	      strncpy_NULL(temp_number,&buff[6],5);
	      strncpy_NULL(name,&buff[13],3);
	      if(atomflag == 'a' || 
		 (atomflag == 'c' && strcmp("CA ",name) == 0) || 
		 (atomflag == 'b' && (strcmp("CA ",name) == 0 || strcmp("C  ",name) == 0 || strcmp("O  ",name) == 0 || strcmp("N  ",name) ==0)))
		{
		  //( printf("%s\n",name);
		  //printf("%s %s",m[0].filename,&buff);
		  strncpy_NULL(alt_loc,&buff[16],1);
		  strncpy_NULL(residue,&buff[17],3);
		  strncpy_NULL(chain,&buff[21],1);
		  strncpy_NULL(temp_resnum,&buff[22],4);
		  strncpy_NULL(resname,&buff[22],5);
		  strncpy_NULL(x_temp,&buff[30],8);
		  strncpy_NULL(y_temp,&buff[38],8);
		  strncpy_NULL(z_temp,&buff[46],8);
		  number=atoi(temp_number);
		  resnum=atoi(temp_resnum);
		  x=atof(x_temp);
		  y=atof(y_temp);
		  z=atof(z_temp);
	      
	      //printf("test: %s %d %s %s %s %s %d %s %lf %lf %lf\n",line_flag,number,name,alt_loc,residue,chain,resnum,resname,x,y,z);
	      //if (strcmp("N",name)==0)	   /*  Is it an N atom => new residue? */
		  //printf("%s %d %d\n",m[0].filename,resnum,ignore_res[resnum]);
		  //		  printf("%d %d\n",residues,ignore_res[0]); 
		  if(ignore_res[resnum]==0) 
		    {
		    if(strcmp(old_resname,resname)!=0)
		      {
			m[0].sequence=realloc(m[0].sequence,sizeof(char)*(residues+1));
			//m[0].sequence=malloc(sizeof(char)*2000);
			m[0].sequence[residues]=aa321(residue);
			residues++;
			strcpy(alt_loc_check,alt_loc);
		      }
		    //sscanf(&buff[22],"%d %lf %lf %lf",&resnum,&x,&y,&z);
		    //printf("test %d %s %s %d %lf %lf %lf\n",number,name,residue,resnum,x,y,z);
		    if(strcmp(alt_loc_check,alt_loc)==0)
		      {
			//printf("test: %s %s",m[0].filename,&buff);
			m[0].atm=realloc(m[0].atm,sizeof(atm)*(atoms+1));
			m[0].atm[atoms].x=x;
			m[0].atm[atoms].y=y;
			m[0].atm[atoms].z=z;
			m[0].atm[atoms].resnum=resnum;
			m[0].atm[atoms].number=number; //atoi(number);
			m[0].atm[atoms].rescount=residues;
			m[0].atm[atoms].selected=TRUE;
			strcpy(m[0].atm[atoms].name,name);
			strcpy(m[0].atm[atoms].residue,residue);
			strcpy(m[0].atm[atoms].resname,resname);
			if(strcmp("CA ",name) == 0)
			  {
			    //printf("%s %s %s %s %d\n",m[0].filename,resname,residue,alt_loc,residues);
			    m[0].CA_ref=realloc(m[0].CA_ref,sizeof(int)*(residues+1));
			    m[0].CA_ref[residues-1]=atoms;
			  }
			//if(strcmp(old_resname,resname)!=0)
			//	{
			//	  m[0].sequence[residues-1]=aa321(residue);
			//	}
			
			m[0].xcen+=x;
			m[0].ycen+=y;
			m[0].zcen+=z;
			atoms++;
			strcpy(old_resname,resname);
		      }
		  }
		}
	    }
	}
      if(atoms>0)  {
	m[0].sequence=realloc(m[0].sequence,sizeof(char)*(residues+1));
	m[0].sequence[residues]='\0';
	m[0].atoms=atoms;
	m[0].residues=residues;
	m[0].xcen=m[0].xcen/atoms;
	m[0].ycen=m[0].ycen/atoms;
	m[0].zcen=m[0].zcen/atoms;
	if(ss_flag == 0)
	  {
	    //fprintf(stderr,"No secondary structure information in file!\n");
	  }
      }
      fclose(fp);
      //#ifdef NOZLIB
      //  fclose(fp);
      //#else
      //	  gzclose(fp);
      //#endif
	  //	  if (atoms!=m[0].atoms)		/* Are file sizes indentical? */
	  //{
	  //  printf("Inconsistent number of atoms in file %s\n",m[0].filename);
	  //  return(1);
	  //}
       
    }
  else
    {
      printf("Couldn't open file \"%s\"\n",m[0].filename);
      exit(1);
    }
  return(0);
}
int free_dyn_molecule(dyn_molecule *m)
{
  free(m[0].sequence);
  free(m[0].ss);
  free(m[0].method);
  free(m[0].CA_ref);
  free(m[0].atm);
}


int read_molecules(molecule *m,char atomflag)	/* Reads in molecule */
{
  int	i,j,k,atoms,residues;	/* Counter variables */
  char	buff[3000];	/* Input string */
  char	temp[3000];	/* Throwaway string */
  char	line_flag[3000];	/* PDB file line mode */
  char  desc_flag[3000];
  char	residue[4];	/* PDB atom info for output */
  char	name[4];
  int   resnum;
  char  resname[9];
  char  old_resname[9]="noname";
  char alt_loc_check[2]=" ";
  char x_temp[11];
  char y_temp[11];
  char z_temp[11];
  //char	number[8];
  int number;
  char chain[2];
  char inscode[2]=" ";
  char alt_loc[2]=" ";
  double	x,y,z;		/* Temporary coordinates values */
  FILE *fp;
  double bfactor=9.99;
  char temp_number[11];
  char temp_resnum[11];
  char temp_bfactor[11];
  int ss_flag=0;
  
  i=0; /* It is only one molecule to be read this was done for the moment instead of changing all "i" to 0 */
  

  //for (i=0;i<2;i++)
  //{
  //#ifdef  NOZLIB*/
  //fp=fopen(m[i].filename,"r");	/* Does file exist? */
  //#else
  //fp=gzopen(m[i].filename,"r");	/* Does file exist? */
  //#endif
  //printf("%s\n",m[0].filename);
  strcpy(m[0].method,"undef");
  m[0].rank=-1;
  m[0].score=-9999;

  fp=fopen(m[0].filename,"r");	/* Does file exist? */
  if (fp!=NULL)	/* If yes, read in coordinates */
    {
      /* Initialize things */
      //m[0].xcen=m[0].ycen=m[0].zcen=0;
      atoms=0;
      residues=0;
      //#ifdef NOZLIB
      //while(fgets(buff,255,fp)!=NULL)
      //#else
      //while(gzgets(fp,buff,255)!=NULL)
      //#endif
      //printf("file: %s\n",m[0].filename);
      while(fgets(buff,1000,fp)!=NULL)
	{
	  //sscanf(buff,"%-6s %4d  %-3s%1s%3s %1s%5s    %7.3lf %7.3lf %7.3lf",line_flag,&number,name,residue,chain,&resnum,inscode,&x,&y,&z);
	  //sscanf(buff,"%s %d %s %s %s",line_flag,&number,name,residue,resname);
	  strcpy(line_flag,"undef");
	  strcpy(desc_flag,"undef");
	  sscanf(buff,"%s",line_flag);
	  
	  if(strcmp("REMARK",line_flag)==0)
	    {
	      //printf("%s %s %s",m[0].filename,line_flag,&buff);
	      sscanf(buff,"%s %s %s",line_flag,desc_flag,temp);
	      //printf("%s %s %s\n",line_flag,desc_flag,temp);
	      if(strcmp("SS",desc_flag)==0)
		{
		  strcpy(m[0].ss,temp);
		  ss_flag=1;
		  //printf("%s\n",m[0].ss);
		}
	      if(strcmp("METHOD",desc_flag)==0)
		{
		  strcpy(m[0].method,temp);
		  //printf("%s\n",m[0].method);
		}
	      if(strcmp("SCORE",desc_flag)==0)
		{
		  m[0].score=atof(temp);
		  //strcpy(m[0].method,temp);
		  //printf("%lf\n",m[0].score);
		}
 	    }


	  if(strcmp("MODEL",line_flag)==0)
	    {
	      sscanf(buff,"%s %s",line_flag,temp);
	      m[0].rank=atoi(temp);
	    }
	  
	  if(strcmp("ATOM",line_flag)==0 && (buff[12] != 'H' || buff[13] != 'H' ))	/* Is it an ATOM entry? */
	    {
	      //printf("Heja: %s %s",m[0].filename,&buff);
	      strncpy_NULL(temp_number,&buff[6],5);
	      strncpy_NULL(name,&buff[13],3);
	      if(atomflag == 'a' || 
		 (atomflag == 'c' && strcmp("CA ",name) == 0) || 
		 (atomflag == 'b' && (strcmp("CA ",name) == 0 || strcmp("C  ",name) == 0 || strcmp("O  ",name) == 0 || strcmp("N  ",name) ==0)))
		{
		  //( printf("%s\n",name);
		  strncpy_NULL(alt_loc,&buff[16],1);
		  strncpy_NULL(residue,&buff[17],3);
		  strncpy_NULL(chain,&buff[21],1);
		  strncpy_NULL(temp_resnum,&buff[22],4);
		  strncpy_NULL(resname,&buff[22],5);
		  strncpy_NULL(x_temp,&buff[30],8);
		  strncpy_NULL(y_temp,&buff[38],8);
		  strncpy_NULL(z_temp,&buff[46],8);
		  
		  number=atoi(temp_number);
		  resnum=atoi(temp_resnum);
		  if(strlen(buff)>=66)
		    {

		      strncpy_NULL(temp_bfactor,&buff[60],6);
		      bfactor=atof(temp_bfactor);
		      //    printf("%lf\n",bfactor);
		    }

		  x=atof(x_temp);
		  y=atof(y_temp);
		  z=atof(z_temp);
	      
	      //printf("test: %s %d %s %s %s %s %d %s %lf %lf %lf\n",line_flag,number,name,alt_loc,residue,chain,resnum,resname,x,y,z);
	      //if (strcmp("N",name)==0)	   /*  Is it an N atom => new residue? */
	      //printf("%s %s\n",old_resname,resname);
		  if(strcmp(old_resname,resname)!=0)
		    {
		      m[0].sequence[residues]=aa321(residue);
		      residues++;
		      strcpy(alt_loc_check,alt_loc);
		      //printf("%s %s\n",resname,residue);
		    }
		  //sscanf(&buff[22],"%d %lf %lf %lf",&resnum,&x,&y,&z);
		  //printf("test %d %s %s %d %lf %lf %lf\n",number,name,residue,resnum,x,y,z);
		  if(strcmp(alt_loc_check,alt_loc)==0)
		    {
		      m[0].atm[atoms].x=x;
		      m[0].atm[atoms].y=y;
		      m[0].atm[atoms].z=z;
		      m[0].atm[atoms].bfactor=bfactor;
		      m[0].atm[atoms].resnum=resnum;
		      m[0].atm[atoms].number=number; //atoi(number);
		      m[0].atm[atoms].rescount=residues;
		      m[0].atm[atoms].selected=TRUE;
		      m[0].atm[atoms].deleted=FALSE;
		      strcpy(m[0].atm[atoms].name,name);
		      strcpy(m[0].atm[atoms].residue,residue);
		      strcpy(m[0].atm[atoms].resname,resname);
		      strcpy(m[0].atm[atoms].chain,chain);
		      
		      if(strcmp("CA ",name) == 0)
			{
			  //printf("%s %s %s %s %d\n",m[0].filename,resname,residue,alt_loc,residues);
			  m[0].CA_ref[residues-1]=atoms;
			}
		      //if(strcmp(old_resname,resname)!=0)
		      //	{
		      //	  m[0].sequence[residues-1]=aa321(residue);
		      //	}
		  
		      m[0].xcen+=x;
		      m[0].ycen+=y;
		      m[0].zcen+=z;
		      atoms++;
		      strcpy(old_resname,resname);
		    }
		}
	    }
	}
      m[0].sequence[residues]='\0';
      m[0].atoms=atoms;
      m[0].residues=residues;
      m[0].xcen=m[0].xcen/atoms;
      m[0].ycen=m[0].ycen/atoms;
      m[0].zcen=m[0].zcen/atoms;
      if(ss_flag == 0)
	{
	  //fprintf(stderr,"No secondary structure information in file!\n");
	}

      fclose(fp);
      //#ifdef NOZLIB
      //  fclose(fp);
      //#else
      //	  gzclose(fp);
      //#endif
	  //	  if (atoms!=m[0].atoms)		/* Are file sizes indentical? */
	  //{
	  //  printf("Inconsistent number of atoms in file %s\n",m[0].filename);
	  //  return(1);
	  //}
     
    }
  else
    {
      printf("Couldn't open file \"%s\"\n",m[0].filename);
      return(1); //exit(1);
    }
  return(0);
}

int read_molecules_ca(molecule *m)	/* Reads in molecules to be superimposed */
{
  int	i,j,k,atoms,residues;	/* Counter variables */
  char	buff[1200];	/* Input string */
  char	temp[1000];	/* Throwaway string */
  char	line_flag[11];	/* PDB file line mode */
  char  desc_flag[20];
  char	residue[4];	/* PDB atom info for output */
  char	name[4];
  int   resnum;
  char  resname[9];
  char  old_resname[9]="noname";
  char alt_loc_check[2]=" ";
  char x_temp[11];
  char y_temp[11];
  char z_temp[11];
  //char	number[8];
  int number;
  char chain[2];
  char inscode[2]=" ";
  char alt_loc[2]=" ";
  double	x,y,z;		/* Temporary coordinates values */
  FILE *fp;
  char temp_number[11];
  char temp_resnum[11];
  
  i=0; /* It is only one molecule to be read this was done for the moment instead of changing all "i" to 0 */
  

  //for (i=0;i<2;i++)
  //{
  //#ifdef  NOZLIB*/
  //fp=fopen(m[i].filename,"r");	/* Does file exist? */
  //#else
  //fp=gzopen(m[i].filename,"r");	/* Does file exist? */
  //#endif
  //printf("%s\n",m[0].filename);
  strcpy(m[0].method,"undef");
  m[0].rank=-1;
  m[0].score=-9999;
  fp=fopen(m[0].filename,"r");	/* Does file exist? */
  if (fp!=NULL)	/* If yes, read in coordinates */
    {
      /* Initialize things */
      //m[0].xcen=m[0].ycen=m[0].zcen=0;
      atoms=0;
      residues=0;
      //#ifdef NOZLIB
      //while(fgets(buff,255,fp)!=NULL)
      //#else
      //while(gzgets(fp,buff,255)!=NULL)
      //#endif
      while(fgets(buff,1000,fp)!=NULL)
	{
	  //sscanf(buff,"%-6s %4d  %-3s%1s%3s %1s%5s    %7.3lf %7.3lf %7.3lf",line_flag,&number,name,residue,chain,&resnum,inscode,&x,&y,&z);
	  //sscanf(buff,"%s %d %s %s %s",line_flag,&number,name,residue,resname);
	  sscanf(buff,"%s",line_flag);
	  if(strcmp("REMARK",line_flag)==0)
	    {
	      //printf("%s",buff);
	      sscanf(buff,"%s %s %s",line_flag,desc_flag,temp);
	      if(strcmp("SS",desc_flag)==0)
		{
		  //printf("%s\n%d\n\n",buff,strlen(buff));
		  strcpy(m[0].ss,temp);
		  //printf("%s\n%d\n",m[0].ss,strlen(m[0].ss));
		}
	      if(strcmp("METHOD",desc_flag)==0)
		{
		  strcpy(m[0].method,temp);
		  //printf("%s\n",m[0].method);
		}
	      if(strcmp("SCORE",desc_flag)==0)
		{
		  m[0].score=atof(temp);
		  //strcpy(m[0].method,temp);
		  //printf("%lf\n",m[0].score);
		}
    
	    }
	  if(strcmp("MODEL",line_flag)==0)
	    {
	      sscanf(buff,"%s %s",line_flag,temp);
	      m[0].rank=atoi(temp);
	    }

	  if (strcmp("ATOM",line_flag)==0 && buff[13] != 'H')	/* Is it an ATOM entry? */
	    { 
	      strncpy_NULL(temp_number,&buff[6],5);
	      strncpy_NULL(name,&buff[13],3);
	      //printf("%s",&buff[6]);
	      if(strcmp("CA ",name) == 0)
		{
	      //printf("%s",&buff[6]);
		  strncpy_NULL(alt_loc,&buff[16],1);
		  strncpy_NULL(residue,&buff[17],3);
		  strncpy_NULL(chain,&buff[21],1);
		  strncpy_NULL(temp_resnum,&buff[22],4);
		  strncpy_NULL(resname,&buff[22],5);
		  strncpy_NULL(x_temp,&buff[30],8);
		  strncpy_NULL(y_temp,&buff[38],8);
		  strncpy_NULL(z_temp,&buff[46],8);
		  
		  number=atoi(temp_number);
		  resnum=atoi(temp_resnum);
		  x=atof(x_temp);
		  y=atof(y_temp);
		  z=atof(z_temp);
		
	      
	      //printf("test: %s %d %s %s %s %s %d %s %lf %lf %lf\n",line_flag,number,name,alt_loc,residue,chain,resnum,resname,x,y,z);
	      //if (strcmp("N",name)==0)	   /*  Is it an N atom => new residue? */
	      //printf("%s %s\n",old_resname,resname);
		  if(strcmp(old_resname,resname)!=0)
		  {
		    m[0].sequence[residues]=aa321(residue);
		    residues++;
		    strcpy(alt_loc_check,alt_loc);
		    //printf("%s %s\n",resname,residue);
		  }
	      //sscanf(&buff[22],"%d %lf %lf %lf",&resnum,&x,&y,&z);
		  //printf("test %d %s %s %d %lf %lf %lf\n",number,name,residue,resnum,x,y,z);
		  if(strcmp(alt_loc_check,alt_loc)==0)
		    {
		      m[0].atm[atoms].x=x;
		      m[0].atm[atoms].y=y;
		      m[0].atm[atoms].z=z;
		      m[0].atm[atoms].resnum=resnum;
		      m[0].atm[atoms].number=number; //atoi(number);
		      m[0].atm[atoms].rescount=residues;
		      m[0].atm[atoms].selected=TRUE;
		      strcpy(m[0].atm[atoms].name,name);
		      strcpy(m[0].atm[atoms].residue,residue);
		      //if(strcmp(old_resname,resname)!=0)
		      //	{
		      //	  m[0].sequence[residues-1]=aa321(residue);
		      //	}
		      
		      m[0].xcen+=x;
		      m[0].ycen+=y;
		      m[0].zcen+=z;
		      atoms++;
		      strcpy(old_resname,resname);
		    }
		}
	    }
	}
      m[0].sequence[residues]='\0';
      m[0].atoms=atoms;
      m[0].residues=residues;
      m[0].xcen=m[0].xcen/atoms;
      m[0].ycen=m[0].ycen/atoms;
      m[0].zcen=m[0].zcen/atoms;
      fclose(fp);
      //#ifdef NOZLIB
      //  fclose(fp);
      //#else
      //	  gzclose(fp);
      //#endif
	  //	  if (atoms!=m[0].atoms)		/* Are file sizes indentical? */
	  //{
	  //  printf("Inconsistent number of atoms in file %s\n",m[0].filename);
	  //  return(1);
	  //}
       
    }
  else
    {
      printf("Couldn't open file \"%s\"\n",m[0].filename);
      exit(1);
    }
  return(0);
}


int read_molecules_backbone(molecule *m)	/* Reads in molecules to be superimposed */
{
  int	i,j,k,atoms,residues;	/* Counter variables */
  char	buff[1200];	/* Input string */
  char	temp[1000];	/* Throwaway string */
  char	line_flag[11];	/* PDB file line mode */
  char  desc_flag[20];
  char	residue[4];	/* PDB atom info for output */
  char	name[4];
  int   resnum;
  char  resname[9];
  char  old_resname[9]="noname";
  char x_temp[11];
  char y_temp[11];
  char z_temp[11];
  //char	number[8];
  int number;
  char chain[2];
  char inscode[2]=" ";
  char alt_loc[2]=" ";
  double	x,y,z;		/* Temporary coordinates values */
  FILE *fp;
  char temp_number[11];
  char temp_resnum[11];
  i=0; /* It is only one molecule to be read this was done for the moment instead of changing all "i" to 0 */
  //printf("%s\n",m[0].filename);
  strcpy(m[0].method,"undef");
  m[0].rank=-1;
  m[0].score=-9999;

  fp=fopen(m[0].filename,"r");	/* Does file exist? */
  if (fp!=NULL)	/* If yes, read in coordinates */
    {
      /* Initialize things */
       atoms=0;
      residues=0;
       while(fgets(buff,1000,fp)!=NULL)
	{
	  //sscanf(buff,"%-6s %4d  %-3s%1s%3s %1s%5s    %7.3lf %7.3lf %7.3lf",line_flag,&number,name,residue,chain,&resnum,inscode,&x,&y,&z);
	  //sscanf(buff,"%s %d %s %s %s",line_flag,&number,name,residue,resname);
	  sscanf(buff,"%s",line_flag);
	  if(strcmp("REMARK",line_flag)==0)
	    {
	      sscanf(buff,"%s %s %s",line_flag,desc_flag,temp);
	      if(strcmp("SS",desc_flag)==0)
		{
		  //strcpy(m[0].ss,temp);
		  strncpy_NULL(m[0].ss,temp,strlen(temp));
		  //  printf("%s\n",m[0].ss);
		}
	      if(strcmp("METHOD",desc_flag)==0)
		{
		  //strcpy(m[0].method,temp);
		  strncpy_NULL(m[0].method,temp,strlen(temp));
		  //printf("%s\n",m[0].method);
		}
	      if(strcmp("SCORE",desc_flag)==0)
		{
		  m[0].score=atof(temp);
		  //strcpy(m[0].method,temp);
		  //printf("%lf\n",m[0].score);
		}
    
	    }
	  if(strcmp("MODEL",line_flag)==0)
	    {
	      sscanf(buff,"%s %s",line_flag,temp);
	      m[0].rank=atoi(temp);
	    }
	  if(strcmp("ATOM",line_flag)==0)	/* Is it an ATOM entry? */
	    {
	      //printf("%s",&buff[6]);
	      strncpy_NULL(temp_number,&buff[6],5);
	      strncpy_NULL(name,&buff[13],3);
	      //printf("%s",&buff[6]);
	      if(strcmp("CA ",name) == 0 ||
		 strcmp("C  ",name) == 0 || 
		 strcmp("O  ",name) == 0 ||
		 strcmp("N  ",name) ==0)
		{
		  //printf("%d %s",2,&buff[6]);
		  strncpy_NULL(alt_loc,&buff[16],1);
		  strncpy_NULL(residue,&buff[17],3);
		  strncpy_NULL(chain,&buff[21],1);
		  strncpy_NULL(temp_resnum,&buff[22],4);
		  strncpy_NULL(resname,&buff[22],5);
		  strncpy_NULL(x_temp,&buff[30],8);
		  strncpy_NULL(y_temp,&buff[38],8);
		  strncpy_NULL(z_temp,&buff[46],8);


		  //(strcmp("CA ",name) !=0)
		  //
		  number=atoi(temp_number);
		  resnum=atoi(temp_resnum);
		  x=atof(x_temp);
		  y=atof(y_temp);
		  z=atof(z_temp);
		  
		  //printf("test: %s %d %s %s %s %s %d %s %lf %lf %lf\n",line_flag,number,name,alt_loc,residue,chain,resnum,resname,x,y,z);
		  //if (strcmp("N",name)==0)	   /*  Is it an N atom => new residue? */
		  //printf("%s %s\n",old_resname,resname);
		  if(strcmp(old_resname,resname)!=0)
		    {
		      m[0].sequence[residues]=aa321(residue);
		      m[0].res_ref[residues]=atoms;
		      residues++;
		      //printf("%s %s\n",resname,residue);
		    }
		      //printf("%s %s\n",resname,residue);
		      //  }
		  //sscanf(&buff[22],"%d %lf %lf %lf",&resnum,&x,&y,&z);
		      //printf("test %d %s %s %d %lf %lf %lf\n",number,name,residue,resnum,x,y,z);
		  m[0].atm[atoms].x=x;
		  m[0].atm[atoms].y=y;
		  m[0].atm[atoms].z=z;
		  m[0].atm[atoms].resnum=resnum;
		  //printf("%d\n",m[0].atm[atoms].resnum);
		  m[0].atm[atoms].number=number; //atoi(number);
		  m[0].atm[atoms].rescount=residues;
		  m[0].atm[atoms].selected=TRUE;
		  strcpy(m[0].atm[atoms].name,name);
		  strcpy(m[0].atm[atoms].residue,residue);
		  if(strcmp("CA ",name) == 0)
		    {
		      //printf("%s %s %s %d\n",resname,residue,alt_loc,residues);
		      m[0].CA_ref[residues-1]=atoms;
		    }

		  //  if(strcmp(old_resname,resname)!=0)
		  // {
		  //    m[0].sequence[residues-1]=aa321(residue);
		  //  }

		  m[0].xcen+=x;
		  m[0].ycen+=y;
		  m[0].zcen+=z;
		  atoms++;
		  strcpy(old_resname,resname);
		  //}
		}
	    }
	}
       m[0].sequence[residues]='\0';
       m[0].atoms=atoms;
       m[0].residues=residues;
       m[0].xcen=m[0].xcen/atoms;
       m[0].ycen=m[0].ycen/atoms;
       m[0].zcen=m[0].zcen/atoms;


      fclose(fp);
      //#ifdef NOZLIB
      //  fclose(fp);
      //#else
      //	  gzclose(fp);
      //#endif
	  //	  if (atoms!=m[0].atoms)		/* Are file sizes indentical? */
	  //{
	  //  printf("Inconsistent number of atoms in file %s\n",m[0].filename);
	  //  return(1);
	  //}
       
    }
  else
    {
      printf("Couldn't open file \"%s\"\n",m[0].filename);
      exit(1);
    }
  return(0);
}

int read_to_molecule(molecule *m,char** atom_vec,char** residue_vec,int* number_vec,rvec* coord_vec,size_t n)
{
  int	i,atoms,residues;	/* Counter variables */
  int   resnum;
  int   old_resnum=-991999;
  atoms=0;
  residues=0;
  
  for(i=0;i<n;i++)
    {
      //printf("Before if: %d %d %s %s %d %lf %lf %lf\n",n,i,atom_vec[i],residue_vec[i],number_vec[i],10*coord_vec[i][0],10*coord_vec[i][1],10*coord_vec[i][2]);
      if(strcmp("CA",atom_vec[i]) == 0 ||
	 strcmp("C",atom_vec[i]) == 0 || 
	 strcmp("O",atom_vec[i]) == 0 ||
	 strcmp("N",atom_vec[i]) == 0)
	{
	  	  
	  strcpy(m[0].atm[atoms].residue,residue_vec[i]);
	  m[0].atm[atoms].residue[3]='\0';//Garantee to be 3 characters long, gromacs add some strange things.
	  resnum=number_vec[i];
	  if(old_resnum!=resnum)
	    {
	      m[0].sequence[residues]=aa321(m[0].atm[atoms].residue);
	      m[0].res_ref[residues]=atoms;
	      //printf("%d %d %s %c\n",old_resnum,resnum,residue_vec[i],m[0].sequence[residues]);
	      residues++;
	    }
	  m[0].atm[atoms].x=10*coord_vec[i][0];
	  m[0].atm[atoms].y=10*coord_vec[i][1];
	  m[0].atm[atoms].z=10*coord_vec[i][2];
	  m[0].atm[atoms].resnum=resnum;
	  m[0].atm[atoms].rescount=residues;
	  m[0].atm[atoms].selected=TRUE;
	  strcpy(m[0].atm[atoms].name,atom_vec[i]);
	  if(strcmp("CA",atom_vec[i]) == 0)
	    {
	      //printf("CA: %s %d %d %lf %lf %lf\n",m[0].atm[atoms].residue,m[0].atm[atoms].resnum,m[0].atm[atoms].rescount,m[0].atm[atoms].x,m[0].atm[atoms].y,m[0].atm[atoms].z);
	      m[0].CA_ref[residues-1]=atoms;
	    }
	  //printf("%s %s %d %d %lf %lf %lf\n",m[0].atm[atoms].name,m[0].atm[atoms].residue,m[0].atm[atoms].resnum,m[0].atm[atoms].rescount,m[0].atm[atoms].x,m[0].atm[atoms].y,m[0].atm[atoms].z);
	  m[0].xcen+=coord_vec[i][0];
	  m[0].ycen+=coord_vec[i][1];
	  m[0].zcen+=coord_vec[i][2];
	  atoms++;
	  old_resnum=resnum;
	}

    }
  //printf("After\n");
   m[0].sequence[residues]='\0';
   m[0].atoms=atoms;
   m[0].residues=residues;
   m[0].xcen=m[0].xcen/atoms;
   m[0].ycen=m[0].ycen/atoms;
   m[0].zcen=m[0].zcen/atoms; 
   //exit(0);
   return(0);
}

void strncpy_NULL(char *dest, char *src, size_t n)
{
  strncpy(dest, src, n);
  dest[n]='\0';
}

int get_atomtype3(char *name, char *res) 
{
  /* C         - 0
   * N         - 1
   * O         - 2
   * other     - 3
   */
  int type=get_atomtype(name,res);
  
  if(type == 0 || type == 3 || type == 4 || type == 5 || type == 6)
    return 0;
  if(type == 1 || type == 7 || type == 8)
    return 1;
  if(type == 2 || type == 9 || type == 10 || type == 11)
    return 2;
  if(type == 12 || type == 13)
    return 3;

}




int get_atomtype(char *name, char *res)        /* Function which takes a name and residue and return type */
{
  /* Types  - Return number
   *
   * C         - 0
   * N         - 1
   * O         - 2
   * CA        - 3
   * CH3       - 4
   * CH/CH2    - 5
   * C(OO-)    - 6
   * NH        - 7
   * NH2       - 8
   * (C)OO-    - 9
   * =O        - 10
   * OH        - 11
   * S         - 12    
   * OXT       - 13
   */
  //printf("%s %s\n",name,res);
  if(strcmp("C  ",name)==0)
    return 0;
  if(strcmp("N  ",name)==0)
    return 1;
  if(strcmp("O  ",name)==0)
    return 2;
  if(strcmp("CA ",name)==0)
    return 3;
  if((strcmp("ILE",res)==0 && (strcmp("CD1",name)==0 || strcmp("CG2",name)==0)) ||
     (strcmp("LEU",res)==0 && (strcmp("CD1",name)==0 || strcmp("CD2",name)==0)) ||
     (strcmp("MET",res)==0 && strcmp("CE ",name)==0) ||
     (strcmp("THR",res)==0 && strcmp("CG2",name)==0) ||
     (strcmp("ALA",res)==0 && strcmp("CB ",name)==0) ||
     (strcmp("VAL",res)==0 && (strcmp("CG1",name)==0 || strcmp("CG2",name)==0)))
    return 4; 
  if((strcmp("ASP",res)==0 && strcmp("CG ",name)==0) ||
     (strcmp("GLU",res)==0 && strcmp("CG ",name)==0))
    return 6;
  if(strcmp("NE ",name)==0 || strcmp("ND1",name)==0 || strcmp("NE1",name)==0)
    return 7;
  if(strcmp("NH1",name)==0 || strcmp("NH2",name)==0 || strcmp("ND2",name)==0
     || strcmp("NE2",name)==0 || strcmp("NZ ",name)==0)
    return 8;
  if((strcmp("ASP",res)==0 && (strcmp("OD2",name)==0 || strcmp("OD1",name)==0)) ||
     (strcmp("GLU",res)==0 && (strcmp("OE1",name)==0 || strcmp("OE2",name)==0)))
    return 9;
  if((strcmp("ASN",res)==0 && strcmp("OD1",name)==0) ||
     (strcmp("GLN",res)==0 && strcmp("OE1",name)==0))
    return 10; 
  if(strcmp("OG ",name)==0 || strcmp("OH ",name)==0 || strcmp("OG1",name)==0)
    return 11;
  if(strcmp("SG ",name)==0 || strcmp("SD ",name)==0)
    return 12;
  if(strcmp("OXT",name)==0)
    return 13;
  return 5;
}

void print_type(int type_no, FILE *fp)
{
  if(type_no==0)
    fprintf(fp,"C     ");
  if(type_no==1)
    fprintf(fp,"N     ");
  if(type_no==2)
    fprintf(fp,"O     ");
  if(type_no==3)
    fprintf(fp,"CA    ");
  if(type_no==4)
    fprintf(fp,"CH3   ");
  if(type_no==5)
    fprintf(fp,"CH/CH2");
  if(type_no==6)
    fprintf(fp,"C(OO-)");
  if(type_no==7)
    fprintf(fp,"NH    ");
  if(type_no==8)
    fprintf(fp,"NH2   ");
  if(type_no==9)
    fprintf(fp,"(C)OO-");
  if(type_no==10)
    fprintf(fp,"=0    ");
  if(type_no==11)
    fprintf(fp,"OH    ");
  if(type_no==12)
    fprintf(fp,"S     ");
  if(type_no==13)
    fprintf(fp,"OXT   ");
}
int get_res6(char* res)
{
     if(strcmp("ARG",res)==0 || strcmp("LYS",res)==0)
       return 0;
     if(strcmp("ASP",res)==0 || strcmp("GLU",res)==0)
       return 1;
     if(strcmp("HIS",res)==0 || strcmp("PHE",res)==0 || strcmp("TRP",res)==0 || strcmp("TYR",res)==0)
       return 2;
     if(strcmp("ASN",res)==0 || strcmp("GLN",res)==0 || strcmp("SER",res)==0 || strcmp("THR",res)==0)
       return 3;
     if(strcmp("ALA",res)==0 || strcmp("ILE",res)==0 || strcmp("LEU",res)==0 || strcmp("MET",res)==0 || strcmp("VAL",res)==0 || strcmp("CYS",res)==0)
       return 4;
     if(strcmp("GLY",res)==0 || strcmp("PRO",res)==0)
       return 5;
     if(strcmp("R",res)==0 || strcmp("K",res)==0)
       return 0;
     if(strcmp("D",res)==0 || strcmp("E",res)==0)
       return 1;
     if(strcmp("H",res)==0 || strcmp("F",res)==0 || strcmp("W",res)==0 || strcmp("Y",res)==0)
       return 2;
     if(strcmp("N",res)==0 || strcmp("Q",res)==0 || strcmp("S",res)==0 || strcmp("T",res)==0)
       return 3;
     if(strcmp("A",res)==0 || strcmp("I",res)==0 || strcmp("L",res)==0 || strcmp("M",res)==0 || strcmp("V",res)==0 || strcmp("C",res)==0)
       return 4;
     if(strcmp("G",res)==0 || strcmp("P",res)==0)
       return 5;
}
int get_res6_no_pointer(char res)
{
  if(res=='R' || res=='K')
    return 0;
  if(res=='D' || res=='E')
    return 1;
  if(res=='H' || res=='F' || res=='W' || res=='Y')
    return 2;
  if(res=='N' || res=='Q' || res=='S' || res=='T')
    return 3;
  if(res=='A' || res=='I' || res=='L' || res=='M' || res=='V' || res=='C')
    return 4;
  if(res=='G' || res=='P')
    return 5;
}




int get_res(char* res)
{
  if(strcmp("ALA",res)==0)
    return 0;
  if(strcmp("ARG",res)==0)
    return 1;
  if(strcmp("ASN",res)==0)
    return 2;
  if(strcmp("ASP",res)==0)
    return 3;
  if(strcmp("CYS",res)==0)
    return 4;
  if(strcmp("GLN",res)==0)
    return 5;
  if(strcmp("GLU",res)==0)
    return 6;
  if(strcmp("GLY",res)==0)
    return 7;
  if(strcmp("HIS",res)==0)
    return 8;
  if(strcmp("ILE",res)==0)
    return 9;
  if(strcmp("LEU",res)==0)
    return 10;
  if(strcmp("LYS",res)==0)
    return 11;
  if(strcmp("MET",res)==0)
    return 12;
  if(strcmp("PHE",res)==0)
    return 13;
  if(strcmp("PRO",res)==0)
    return 14;
  if(strcmp("SER",res)==0)
    return 15;
  if(strcmp("THR",res)==0)
    return 16;
  if(strcmp("TRP",res)==0)
    return 17;
  if(strcmp("TYR",res)==0)
    return 18;
  if(strcmp("VAL",res)==0)
    return 19;
}
void print_res(int res,FILE *fp)
{
  if(res==0)
    fprintf(fp,"ALA");
  if(res == 1)
    fprintf(fp,"ARG");
  if(res == 2)
    fprintf(fp,"ASN");
  if(res == 3)
    fprintf(fp,"ASP");
  if(res == 4)
    fprintf(fp,"CYS");
  if(res == 5)
    fprintf(fp,"GLN");
  if(res == 6)
    fprintf(fp,"GLU");
  if(res == 7)
    fprintf(fp,"GLY");
  if(res == 8)
    fprintf(fp,"HIS");
  if(res == 9)
    fprintf(fp,"ILE");
  if(res == 10)
    fprintf(fp,"LEU");
  if(res == 11)
    fprintf(fp,"LYS");
  if(res == 12)
    fprintf(fp,"MET");
  if(res == 13)
    fprintf(fp,"PHE");
  if(res == 14)
    fprintf(fp,"PRO");
  if(res == 15)
    fprintf(fp,"SER");
  if(res == 16)
    fprintf(fp,"THR");
  if(res == 17)
    fprintf(fp,"TRP");
  if(res == 18)
    fprintf(fp,"TYR");
  if(res == 19)
    fprintf(fp,"VAL");

}

double crd(molecule *m,int atomno1, int atomno2)       /* atomnoX is the first atom of the residue */
{
  int i,j;
  double dist,lowest_dist;
  lowest_dist=999999;
  for(i=atomno1;m[0].atm[i].rescount == m[0].atm[atomno1].rescount;i++)  
    {
      if(strcmp("C  ",m[0].atm[i].name)!=0 && 
	 strcmp("O  ",m[0].atm[i].name)!=0 &&
	 strcmp("N  ",m[0].atm[i].name)!=0)
	{
	  for(j=atomno2;m[0].atm[j].rescount == m[0].atm[atomno2].rescount;j++)
	    {
	      if(strcmp("C  ",m[0].atm[j].name)!=0 && 
		 strcmp("O  ",m[0].atm[j].name)!=0 &&
		 strcmp("N  ",m[0].atm[j].name)!=0)
		{
		  dist=distance(m,i,j);
		  //printf("%f ",dist);
		  if(dist<lowest_dist)
		    lowest_dist=dist;
		  //	  printf("%s %s %f\n",m[0].atm[j].residue,m[0].atm[j].name,lowest_dist);
		}
	    }
	}
    }
  //printf("%s %s %d %d\n",m[0].atm[i].residue,m[0].atm[i].name,m[0].atm[i].rescount,m[0].atm[atomno1].rescount);
  //printf("%s %s %d %d\n",m[0].atm[j].residue,m[0].atm[j].name,m[0].atm[j].rescount,m[0].atm[atomno2].rescount);
  return lowest_dist;
}





double distance(molecule *m,int atomno1,int atomno2)
{
  //printf("%s %s\n",m[0].atm[atomno1-1].name,m[0].atm[atomno2-1].name);

  return (m[0].atm[atomno1].x-m[0].atm[atomno2].x)*(m[0].atm[atomno1].x-m[0].atm[atomno2].x)+(m[0].atm[atomno1].y-m[0].atm[atomno2].y)*(m[0].atm[atomno1].y-m[0].atm[atomno2].y)+(m[0].atm[atomno1].z-m[0].atm[atomno2].z)*(m[0].atm[atomno1].z-m[0].atm[atomno2].z);
}

//void jacobi(float **a, int n, float d[], float **v, int *nrot);

double fatness(molecule *m)
{
  int i,j,nrot;
  float  trans_x,trans_y,trans_z,trans_x2,trans_y2,trans_z2;
  float a,b,c,max,min;
  float  V[6]={};
  float  I[3][3]={};
  float *d,**v,**e;
  int NP=3;

  for(i=0;i<m[0].atoms;i++)
    {
      trans_x=(float)m[0].atm[i].x-m[0].xcen;
      trans_y=(float)m[0].atm[i].y-m[0].ycen;
      trans_z=(float)m[0].atm[i].z-m[0].zcen;

      trans_x2=trans_x*trans_x;
      trans_y2=trans_y*trans_y;
      trans_z2=trans_z*trans_z;
      
      V[0]=V[0]+trans_y2+trans_z2;
      V[1]=V[1]+trans_x*trans_y;
      V[2]=V[2]+trans_x*trans_z;

      
      V[3]=V[3]+trans_x2+trans_z2;
      V[4]=V[4]+trans_y*trans_z;
      
      V[5]=V[5]+trans_x2+trans_y2;
    }
  I[0][0]=V[0];
  I[1][0]=-V[1];
  I[2][0]=-V[2];
  I[0][1]=-V[1];
  I[1][1]=V[3];
  I[2][1]=-V[4];
  I[0][2]=-V[2];
  I[1][2]=-V[4];
  I[2][2]=V[5];

  
  
  d=vector(1,NP);
  v=matrix_nr(1,NP,1,NP);
  e=convert_matrix(&I[0][0],1,NP,1,NP);
  //printf("%5.3f %5.3f %5.3f\n",e[1][1],e[1][2],e[1][3]);
  //printf("%5.3f %5.3f %5.3f\n",e[2][1],e[2][2],e[2][3]);
  //printf("%5.3f %5.3f %5.3f\n",e[3][1],e[3][2],e[3][3]);
#ifdef GMXLIB
  jacobi_single(e,NP,d,v,&nrot);
#else
  jacobi(e,NP,d,v,&nrot);
#endif
  //eigsrt(d,v,NP);
  a=sqrt(2.5*(d[2]+d[3]-d[1])/m[0].atoms);
  b=sqrt(2.5*(d[1]+d[3]-d[2])/m[0].atoms);
  c=sqrt(2.5*(d[1]+d[2]-d[3])/m[0].atoms);

  //printf("%5.3f %5.3f %5.3f\n=======\n",a,b,c); 
  min=a;
  if(b<min)
    min=b;
  if(c<min)
    min=c;

  max=a;
  if(b>max)
    max=b;
  if(c>max)
    max=c;
  
  //printf("%5.3f %5.3f %5.3f\n=======\n",a,b,c); 
  //printf("%5.3f",max/min); 
  
  

  
  //printf("unsorted eigenvectors:\n");
  //for (i=1;i<=NP;i++) {
  //  printf("eigenvalue %3d = %12.6f\n",i,d[i]);
  //  printf("eigenvector:\n");
  //  for (j=1;j<=NP;j++) {
  //	printf("%12.6f",v[j][i]);
  //	if ((j % 5) == 0) printf("\n");
  //  }
  //  printf("\n");
  //}
  //printf("\n****** Sorting Eigenvectors ******\n\n");
  //eigsrt(d,v,NP);
  //printf("sorted eigenvectors:\n");
  //for (i=1;i<=NP;i++) {
  //  printf("eigenvalue %3d = %12.6f\n",i,d[i]);
  //  printf("eigenvector:\n");
  //  for (j=1;j<=NP;j++) {
  //	printf("%12.6f",v[j][i]);
  //	if ((j % 5) == 0) printf("\n");
  //  }
  //  printf("\n");
  //}
  free_convert_matrix(e,1,NP,1,NP);
  free_matrix(v,1,NP,1,NP);
  free_vector(d,1,NP);
  
  return max/min;
}
double fatness2(molecule *m)
{
  int i,j,nrot;
  float  trans_x,trans_y,trans_z,trans_x2,trans_y2,trans_z2;
  float a,b,c,max,min;
  float  V[6]={};
  float  I[3][3]={};
  float *d,**v,**e;
  int NP=3;

  for(i=0;i<m[0].atoms;i++)
    {
      trans_x=(float)m[0].atm[i].x-m[0].xcen;
      trans_y=(float)m[0].atm[i].y-m[0].ycen;
      trans_z=(float)m[0].atm[i].z-m[0].zcen;

      trans_x2=trans_x*trans_x;
      trans_y2=trans_y*trans_y;
      trans_z2=trans_z*trans_z;
      
      V[0]=V[0]+trans_y2+trans_z2;
      V[1]=V[1]+trans_x*trans_y;
      V[2]=V[2]+trans_x*trans_z;

      
      V[3]=V[3]+trans_x2+trans_z2;
      V[4]=V[4]+trans_y*trans_z;
      
      V[5]=V[5]+trans_x2+trans_y2;
    }
  I[0][0]=V[0];
  I[1][0]=-V[1];
  I[2][0]=-V[2];
  I[0][1]=-V[1];
  I[1][1]=V[3];
  I[2][1]=-V[4];
  I[0][2]=-V[2];
  I[1][2]=-V[4];
  I[2][2]=V[5];

  
  
  d=vector(1,NP);
  v=matrix_nr(1,NP,1,NP);
  e=convert_matrix(&I[0][0],1,NP,1,NP);
  //printf("%5.3f %5.3f %5.3f\n",e[1][1],e[1][2],e[1][3]);
  //printf("%5.3f %5.3f %5.3f\n",e[2][1],e[2][2],e[2][3]);
  //printf("%5.3f %5.3f %5.3f\n",e[3][1],e[3][2],e[3][3]);
#ifdef GMXLIB
  jacobi_single(e,NP,d,v,&nrot);
#else
  jacobi(e,NP,d,v,&nrot);
#endif
  //eigsrt(d,v,NP);
  a=sqrt(2.5*(d[2]+d[3]-d[1])/m[0].atoms);
  b=sqrt(2.5*(d[1]+d[3]-d[2])/m[0].atoms);
  c=sqrt(2.5*(d[1]+d[2]-d[3])/m[0].atoms);

  //printf("%5.3f %5.3f %5.3f\n=======\n",a,b,c); 
  min=a;
  if(b<min)
    min=b;
  if(c<min)
    min=c;

  max=a;
  if(b>max)
    max=b;
  if(c>max)
    max=c;
  
  //printf("%5.3f %5.3f %5.3f\n=======\n",a,b,c); 
  //printf("%5.3f",max/min); 
  
  

  
  //printf("unsorted eigenvectors:\n");
  //for (i=1;i<=NP;i++) {
  //  printf("eigenvalue %3d = %12.6f\n",i,d[i]);
  //  printf("eigenvector:\n");
  //  for (j=1;j<=NP;j++) {
  //	printf("%12.6f",v[j][i]);
  //	if ((j % 5) == 0) printf("\n");
  //  }
  //  printf("\n");
  //}
  //printf("\n****** Sorting Eigenvectors ******\n\n");
  //eigsrt(d,v,NP);
  //printf("sorted eigenvectors:\n");
  //for (i=1;i<=NP;i++) {
  // printf("eigenvalue %3d = %12.6lf\n",i,d[i]);
  // printf("eigenvector:\n");
  // for (j=1;j<=NP;j++) {
  //   printf("%12.6f",v[j][i]);
  //   if ((j % 5) == 0) printf("\n");
  // }
  // printf("\n");
  //}
  free_convert_matrix(e,1,NP,1,NP);
  free_matrix(v,1,NP,1,NP);
  free_vector(d,1,NP);
  
  return max;
}


char aa321(const char* res)
{
  //char check_res[4];
   //check_res[0]=res[0];
  //check_res[1]=res[1];
  //check_res[2]=res[2];
  //res[3]='\0';
  //printf("%s %s\n",res,res);
  
  if(strcmp("ALA",res)==0)
    return 'A';
  if(strcmp("ARG",res)==0)
    return 'R';
  if(strcmp("ASN",res)==0)
    return 'N';
  if(strcmp("ASP",res)==0)
    return 'D';
  if(strcmp("CYS",res)==0)
    return 'C';
  if(strcmp("GLN",res)==0)
    return 'Q';
  if(strcmp("GLU",res)==0)
    return 'E';
  if(strcmp("GLY",res)==0)
    return 'G';
  if(strcmp("HIS",res)==0)
    return 'H';
  if(strcmp("ILE",res)==0)
    return 'I';
  if(strcmp("LEU",res)==0)
    return 'L';
  if(strcmp("LYS",res)==0)
    return 'K';
  if(strcmp("MET",res)==0)
    return 'M';
  if(strcmp("PHE",res)==0)
    return 'F';
  if(strcmp("PRO",res)==0)
    return 'P';
  if(strcmp("SER",res)==0)
    return 'S';
  if(strcmp("THR",res)==0)
    return 'T';
  if(strcmp("TRP",res)==0)
    return 'W';
  if(strcmp("TYR",res)==0)
    return 'Y';
  if(strcmp("VAL",res)==0)
    return 'V';
  return 'X';
}

char*  assign_ss(molecule *m,float cutoff,float angle)
{
  char  *ss;
  int   *hb;
  int   j,i;
  int   number_of_res;
  int   *alphatemp,*betatemp1,*betatemp2,*alphatest,*betatest;
  int hbond_check=0;
  int hbond_error=0;
  
  
   
  //printf("Assigning Secondary structure\n");
  //printf("Residues: %d\n",m[0].residues);
  ss=malloc(m[0].residues*sizeof(char)+1);
  alphatemp=malloc(m[0].residues*sizeof(int));
  betatemp1=malloc(m[0].residues*sizeof(int));
  betatemp2=malloc(m[0].residues*sizeof(int));
  alphatest=malloc(m[0].residues*sizeof(int));
  betatest=malloc(m[0].residues*sizeof(int));
  
  hb=malloc(m[0].residues*m[0].residues*sizeof(int));
  //printf("Leffe\n");
  
  for(i=0;i<m[0].residues;i++)
    {
      alphatemp[i]=0;
      betatemp1[i]=0;
      betatemp2[i]=0;
      alphatest[i]=0;
      betatest[i]=0;
    }
  for(i=0;i<m[0].residues*m[0].residues;i++)
    {
      hb[i]=0;
    }
  
  for(i=0;i<m[0].residues;i++)
    {
      //for(j=i+1;j<m[0].residues;j++)
      for(j=0;j<m[0].residues;j++)
	{
	  hbond_check=hbond(m,m[0].res_ref[i],m[0].res_ref[j],cutoff,angle);
	  if(hbond_check==1)
	    {
	      hb[calc_index(m[0].residues,i,j)]=1;
	      //hb[calc_index(m[0].residues,j,i)]=1;
	      //      printf("Hbond found> %d %d\n",i,j);
	    }
	  else
	    {
	      if(hbond_check==-1)
		hbond_error++;
	// if(hbond_error>m[0].residues*m[0].residues)
	//	{
	//	  //fprintf(stderr,"hej %s\n",m[0].filename,hbond_error);
	//	  fprintf(stderr,"%s Unable to assign ss to structure... %d %d %5.3f\n",m[0].filename,hbond_error,m[0].residues,hbond_error/m[0].residues);
	//	  // fprintf(stderr,"%s Unable to assign ss to structure... %d \n",m[0].filename,hbond_error);
	//	  free(hb);
	//	  free(alphatemp);
	//	  free(betatemp1);
	//	  free(betatemp2);
	//	  free(alphatest);
	//	  free(betatest);
	//	  ss[0]='\0';
	//	  return ss;
	//	}
	    }
	    
	}
    }
  if(hbond_error==m[0].residues*m[0].residues)
    fprintf(stderr,"file: %s contains no hydrogen bonds (might be a CA-model)\n",m[0].filename);
  for(i=1;i<m[0].residues-2;i++)
    {
      //printf("%d %d %d %d %d %d %d\n",i,i+3,i-1,i+2,m[0].residues,hbond(m,m[0].res_ref[i],m[0].res_ref[i+3]),hbond(m,m[0].res_ref[i-1],m[0].res_ref[i+2]));
      //  if(hbond(m,m[0].res_ref[i],m[0].res_ref[i+3]) &&
      //	 hbond(m,m[0].res_ref[i-1],m[0].res_ref[i+2]))
      //	{
      //if(hb[i*m[0].residues+i+3] &&
      //	 hb[(i-1)*m[0].residues+i+2])

      if(hb[calc_index(m[0].residues,i,i+3)] &&
	 hb[calc_index(m[0].residues,i-1,i+2)])
	{
	  for(j=i;j<=i+3;j++)
	    {
	      alphatemp[j]=3;
	      //printf("%d %d\n",i,j);
	    }
	}
    }
  for(i=1;i<m[0].residues-3;i++)
    {
      //printf("%d %d %d %d %d\n",i,i+4,i-1,i+3,m[0].residues);
      //if(hbond(m,m[0].res_ref[i],m[0].res_ref[i+4]) &&
      //	 hbond(m,m[0].res_ref[i-1],m[0].res_ref[i+3]))
      //	{
      //if(hb[i*m[0].residues+i+4] &&
      //	 hb[(i-1)*m[0].residues+i+3])
      if(hb[calc_index(m[0].residues,i,i+4)] &&
	 hb[calc_index(m[0].residues,i-1,i+3)])
	{
	  for(j=i;j<=i+3;j++)
	    {
	      alphatemp[j]=4;
	      //     printf("%d %d\n",i,j);
	    }
	}
    }

 
  for(i=0;i<m[0].residues;i++)
    {
      for(j=0;j<m[0].residues;j++)
	 { 
	   if(abs(j-i)>2)
	     {
	       //printf("%d %d %d\n",i,j,m[0].residues);
	       if((hb[calc_index(m[0].residues,max(i-1,0),j)] && hb[calc_index(m[0].residues,j,min(i+1,m[0].residues-1))]) ||
		  (hb[calc_index(m[0].residues,max(j-1,0),i)] && hb[calc_index(m[0].residues,i,min(j+1,m[0].residues-1))]))
		 {
		   
		   betatemp1[max(i-1,0)]=1;
		   betatemp1[i]=1;
		   betatemp1[min(m[0].residues-1,i+1)]=1;
		   betatemp1[max(j-1,0)]=1;
		   betatemp1[j]=1;
		   betatemp1[min(m[0].residues-1,j+1)]=1;
		   
		 }
	       if((hb[calc_index(m[0].residues,i,j)] && hb[calc_index(m[0].residues,j,i)])||
		  (hb[calc_index(m[0].residues,max(i-1,0),min(j+1,m[0].residues-1))] &&
		   hb[calc_index(m[0].residues,max(j-1,0),min(i+1,m[0].residues-1))]))
		 {
		   //printf("BETA\n");
		   betatemp2[max(i-1,0)]=1;
		   betatemp2[i]=1;
		   betatemp2[min(m[0].residues-1,i+1)]=1;
		   betatemp2[max(j-1,0)]=1;
		   betatemp2[j]=1;
		   betatemp2[min(m[0].residues-1,j+1)]=1;
		 }
	     }
	 }
    }

  for(i=0;i<m[0].residues;i++)
    {
      if(alphatemp[i]>1)
	{
	  alphatest[i]=1;
	  //printf("INFO> Alpha-helix %d\n",i);
	}
      if((betatemp1[max(0,i-2)] && betatemp1[i]) ||
	 (betatemp2[max(0,i-2)] && betatemp2[i]))
	{
	  betatest[i]=1;
          betatest[max(0,i-1)]=1;
          betatest[max(0,i-2)]=1;
	  //printf("BETA\n");
	}
      else
	{
	  if((betatemp1[min(m[0].residues-1,i+2)] && betatemp1[i]) ||
	     (betatemp2[min(m[0].residues-1,i+2)] && betatemp2[i]))
	    {
	      betatest[i]=1;
	      betatest[min(m[0].residues-1,i+1)]=1;
	      betatest[min(m[0].residues-1,i+2)]=1;
	      //printf("BETA\n");
	    }
	  

	}
    }
  for(i=0;i<m[0].residues;i++)
    {
      //printf("%s ",m[0].atm[m[0].CA_ref[i]].residue);
      if(alphatest[i]==1)
	{
	  //printf("H ");
	  ss[i]='H';
	}
      else
	{
	  if(betatest[i]==1)
	    {
	      //printf("E ");
	      ss[i]='E';
	    }
	  else
	    {
	      //if(betatest[i]!=1 && alphatest[i]!=1)
	      //	{
	      ss[i]='C';
		  //printf("C ");
	    }
	}
      //printf("\n");
    }
  ss[i]='\0';

  free(hb);
  free(alphatemp);
  free(betatemp1);
  free(betatemp2);
  free(alphatest);
  free(betatest);
  //printf("%s\n%s\n",m[0].sequence,ss);


  return ss;
}

int calc_index(int size,int row,int col)
{
  return size*row+col;
}



int hbond(molecule *m,int atomno1, int atomno2,float cutoff, float angle)     /* atomnoX is the first atom of the residue */
{
  int N1,C1,O1,CA1,N2,C2,O2,CA2,i,j,chksum=0;

  double nvecN1,nvecN2,nvecO1,nvecO2;
  double angle1,angle2;
  my_vector vecN1,vecN2,vecO1,vecO2;
  int debug=0;
  
  //if(strcmp(m[0].filename,"/afs/pdc.kth.se/home/b/bjornw/bjorn/modelling/pcons5/models/LB4/ss_added/1e9sE.1b6g.fugue.7.bb")==0)
  //  {
  //    debug=1;
  //  }

  //Assign indices to the different atomtypes
  if(debug)
    printf("Assigning\n");
  for(i=atomno1;m[0].atm[i].rescount == m[0].atm[atomno1].rescount;i++)  
    {
      if(strcmp("C  ",m[0].atm[i].name)==0 || strcmp("C",m[0].atm[i].name)==0){
	C1=i;
	if(debug)
	  printf("1: %s %s\n",m[0].atm[atomno1].residue,m[0].atm[i].name);
	chksum++;
      }
      if(strcmp("O  ",m[0].atm[i].name)==0 || strcmp("O",m[0].atm[i].name)==0){
	O1=i;
	if(debug)
	  printf("1: %s %s\n",m[0].atm[atomno1].residue,m[0].atm[i].name);
	chksum++;
      }
      if(strcmp("N  ",m[0].atm[i].name)==0 || strcmp("N",m[0].atm[i].name)==0){
	N1=i;
	if(debug)
	  printf("1: %s %s\n",m[0].atm[atomno1].residue,m[0].atm[i].name);
	chksum++;
      }
      if(strcmp("CA ",m[0].atm[i].name)==0 || strcmp("CA",m[0].atm[i].name)==0){
	CA1=i;
	if(debug)
	  printf("1: %s %s\n",m[0].atm[atomno1].residue,m[0].atm[i].name);
	chksum++;
      }
    }
  for(i=atomno2;m[0].atm[i].rescount == m[0].atm[atomno2].rescount;i++)  
    {
      if(strcmp("C  ",m[0].atm[i].name)==0 || strcmp("C",m[0].atm[i].name)==0){
	C2=i;
	if(debug)
	  printf("2: %s %s\n",m[0].atm[atomno2].residue,m[0].atm[i].name);
	chksum++;
      }
      if(strcmp("O  ",m[0].atm[i].name)==0 || strcmp("O",m[0].atm[i].name)==0){
	O2=i;
	if(debug)
	  printf("2: %s %s\n",m[0].atm[atomno2].residue,m[0].atm[i].name);
	chksum++;
      }
      if(strcmp("N  ",m[0].atm[i].name)==0 || strcmp("N",m[0].atm[i].name)==0){
	N2=i;
	if(debug)
	  printf("2: %s %s\n",m[0].atm[atomno2].residue,m[0].atm[i].name);
	chksum++;
      }
      if(strcmp("CA ",m[0].atm[i].name)==0 || strcmp("CA",m[0].atm[i].name)==0){
	CA2=i;
	if(debug)
	  printf("2: %s %s\n",m[0].atm[atomno2].residue,m[0].atm[i].name);
	chksum++;
      }
    }
  //In case of "Alternative location" for atoms the last will be chosen.
  if(debug)
    {
      printf("End Assigning\n");
      printf("%d %d %d %d %d %d %d %d\n",N1,CA1,C1,O1,N2,CA2,C2,O2);
    }

  
  //  vecN1.x=m[0].atm[N1].x-m[0].atm[O2].x;
  //  vecN1.y=m[0].atm[N1].y-m[0].atm[O2].y;
  //  vecN1.z=m[0].atm[N1].z-m[0].atm[O2].z;
  if(chksum<8)
    {
      
      //if(chksum!=2)
      //{
      //	  fprintf(stderr,"Unable to assign all atomtypes in function hbond, %d, %s %d - %s %d %s\n",chksum,m[0].atm[atomno1].residue,m[0].atm[atomno1].resnum,m[0].atm[atomno2].residue,m[0].atm[atomno2].resnum,m[0].filename);
      //	}
      // if(chksum == 2)
      //	{
      //	  fprintf(stderr,"Might be CA model...\n");
      //	}
      return -1;

    }

  vecN2.x=m[0].atm[N2].x-m[0].atm[O1].x;
  vecN2.y=m[0].atm[N2].y-m[0].atm[O1].y;
  vecN2.z=m[0].atm[N2].z-m[0].atm[O1].z;
 

  // Check bond length to be sure..... 1.22-1.25
  vecO1.x=m[0].atm[O1].x-m[0].atm[C1].x;
  vecO1.y=m[0].atm[O1].y-m[0].atm[C1].y;
  vecO1.z=m[0].atm[O1].z-m[0].atm[C1].z;
     
  //vecO2.x=m[0].atm[O2].x-m[0].atm[C2].x;
  //vecO2.y=m[0].atm[O2].y-m[0].atm[C2].y;
  //vecO2.z=m[0].atm[O2].z-m[0].atm[C2].z;
  

  //Direction from C=O to H-N CO1->N2
  // printf("%d %s %s %lf\n",atomno1,m[0].atm[atomno1].name,m[0].atm[atomno1].residue,CO2.x*CO2.x+CO2.y*CO2.y+CO2.z*CO2.z);
  
  //nvecN1=sqrt(vecN1.x*vecN1.x+vecN1.y*vecN1.y+vecN1.z*vecN1.z);
  nvecN2=sqrt(vecN2.x*vecN2.x+vecN2.y*vecN2.y+vecN2.z*vecN2.z);
  //printf("%lf %lf\n",nvecN1,nvecN2);
 
  if(nvecN2<cutoff)
    {
      nvecO1=sqrt(vecO1.x*vecO1.x+vecO1.y*vecO1.y+vecO1.z*vecO1.z);
      angle2=acos((vecN2.x*vecO1.x+vecN2.y*vecO1.y+vecN2.z*vecO1.z)/(nvecO1*nvecN2));

      //      printf("2 %3d %3s %3s %4d %4d %3s %3s %4d %lf %lf %lf\n",atomno1,m[0].atm[atomno1].name,m[0].atm[atomno1].residue,m[0].atm[atomno1].resnum,atomno2,m[0].atm[atomno2].name,m[0].atm[atomno2].residue,m[0].atm[atomno2].resnum,nvecO1,nvecN2,angle2*180/PI);
      if(angle2<angle)  //completary angle
	{
	  return 1;
	  //printf("2 %3d %3s %3s %4d %4d %3s %3s %4d %lf %lf %lf\n",atomno1,m[0].atm[atomno1].name,m[0].atm[atomno1].residue,m[0].atm[atomno1].resnum,atomno2,m[0].atm[atomno2].name,m[0].atm[atomno2].residue,m[0].atm[atomno2].resnum,nvecO2,nvecN2,angle2*180/PI);
	  //return 1;
	}
    }
  
//  if(nvecN1<3.5)   
//    {
//	nvecO2=sqrt(vecO2.x*vecO2.x+vecO2.y*vecO2.y+vecO2.z*vecO2.z);
//	angle1=acos((vecN1.x*vecO2.x+vecN1.y*vecO2.y+vecN1.z*vecO2.z)/(nvecO2*nvecN1));
//	//printf("1 %3d %3s %3s %4d %4d %3s %3s %4d %lf %lf %lf\n",atomno1,m[0].atm[atomno1].name,m[0].atm[atomno1].residue,m[0].atm[atomno1].resnum,atomno2,m[0].atm[atomno2].name,m[0].atm[atomno2].residue,m[0].atm[atomno2].resnum,nvecO2,nvecN1,angle1*180/PI);
//	if(angle1<1.57) //completary angle
//	  {
//	    return 1;
//	    //printf("1 %3d %3s %3s %4d %4d %3s %3s %4d %lf %lf %lf\n",atomno1,m[0].atm[atomno1].name,m[0].atm[atomno1].residue,m[0].atm[atomno1].resnum,atomno2,m[0].atm[atomno2].name,m[0].atm[atomno2].residue,m[0].atm[atomno2].resnum,nvecO1,nvecN1,angle1*180/PI);
//	    // return 1; 
//	  }
//	}
  //printf("%4d %4d %4d %4d %4d %4d %4d %4d %4d %4d\n",atomno1,N1,C1,O1,CA1,atomno2,N2,C2,O2,CA2);
  //printf("%d %s %s %d %s %s\n",atomno1,m[0].atm[atomno1].name,m[0].atm[atomno1].residue,atomno2,m[0].atm[atomno2].name,m[0].atm[atomno2].residue);
  return 0;

}


char*  read_psipred(char* filename)
{
  char  *ss;
  char  *temp;
  char	buff[512];	/* Input string */
  char	line_flag[11];	
  char  str[200];
  FILE *fp;
 
  ss=malloc(2000*sizeof(char));
  strcpy(ss,"");
  
  //printf("len temp: %d\n",strlen(temp));
  //printf("len ss: %d\n",strlen(ss));
  fp=fopen(filename,"r");	/* Does file exist? */
  if (fp!=NULL)
    {
      while(fgets(buff,255,fp)!=NULL)
	{
	  sscanf(buff,"%s %s",line_flag,str);
	  if(strcmp(line_flag,"Pred:")==0)
	    {
	      //printf("%d\n",strlen(ss));
	      //printf("%s",buff);
	      //printf("%s\n",str);
	      strcat(ss,str);

	    }
	  // printf("%s\n",line_flag);
	}
      //   printf("%s\n",ss);
      return ss;
    }
  else
    {
      return 0;
    }
}

char*  read_psipred2(char* filename,double* coil,double* helix,double* sheet)
{
  char  *ss;
  char	buff[512];	/* Input string */
  char	line_flag[11];	
  char  str[200];
  FILE *fp;
  int temp;
  char seq;
  int i=0;
  int j=0;
  char coil_temp[10];
  char helix_temp[10];
  char sheet_temp[10];
  ss=malloc(MAXRES*sizeof(char));
  helix[i]=0.1;
  coil[i]=0.1;
  sheet[i]=0.3;
  //printf("len temp: %d\n",strlen(temp));
  //printf("len ss: %d\n",strlen(ss));
  fp=fopen(filename,"r");	/* Does file exist? */
  if (fp!=NULL)
    {
      while(fgets(buff,255,fp)!=NULL)
	{
	  if(strlen(buff)==31)
	    {
	      //printf("%d %s ",strlen(buff),buff);
	      ss[i]=buff[7];
	      strncpy_NULL(coil_temp,&buff[11],5);
	      strncpy_NULL(helix_temp,&buff[18],5);
	      strncpy_NULL(sheet_temp,&buff[25],5);
	      coil[i]=atof(coil_temp);
	      helix[i]=atof(helix_temp);
	      sheet[i]=atof(sheet_temp);
	      i++;
	    }
	}
    }
  else
    {
      //return 0;
      fprintf(stderr,"Could not open %s\n",filename);
      exit(1);
    }
  //for(j=0;j<i;j++)
  //  {
  //    printf("%c %lf %lf %lf\n",ss[j],coil[j],helix[j],sheet[j]);
  //  }
  ss[i]='\0';
  return ss;
}

int aa(char aa)
{
  if(aa=='A')
    return 0;
  if(aa=='R')
    return 1;
  if(aa=='N')
    return 2;
  if(aa=='D')
    return 3;
  if(aa=='C')
    return 4;
  if(aa=='Q')
    return 5;
  if(aa=='E')
    return 6;
  if(aa=='G')
    return 7;
  if(aa=='H')
    return 8;
  if(aa=='I')
    return 9;
  if(aa=='L')
    return 10;
  if(aa=='K')
    return 11;
  if(aa=='M')
    return 12;
  if(aa=='F')
    return 13;
  if(aa=='P')
    return 14;
  if(aa=='S')
    return 15;
  if(aa=='T')
    return 16;
  if(aa=='W')
    return 17;
  if(aa=='Y')
    return 18;
  if(aa=='V')
    return 19;
}

char profile_index_to_aa(int i)
{
  if(i==0)
    return 'A';
  if(i==1)
    return 'R';
  if(i==2)
    return 'N';
  if(i==3)
    return 'D';
  if(i==4)
    return 'C';
  if(i==5)
    return 'Q';
  if(i==6)
    return 'E';
  if(i==7)
    return 'G';
  if(i==8)
    return 'H';
  if(i==9)
    return 'I';
  if(i==10)
    return 'L';
  if(i==11)
    return 'K';
  if(i==12)
    return 'M';
  if(i==13)
    return 'F';
  if(i==14)
    return 'P';
  if(i==15)
    return 'S';
  if(i==16)
    return 'T';
  if(i==17)
    return 'W';
  if(i==18)
    return 'Y';
  if(i==19)
    return 'V';
}



void read_profile(char* filename,double (*prof)[22],char* seq)
{
  int i=0;
  int j=0;
  int check=0;
  int seq_i=0;
  int prof_col=0;
  char	buff[512];
  char  temp[15];
  int   temp2=0;
  
  FILE *fp;
  //printf("TEST: %d\n", prof[1][0]);
  //prof[seq_i][prof_col]=0;
  fp=fopen(filename,"r");	/* Does file exist? */
  if (fp!=NULL)
    {
      while(fgets(buff,255,fp)!=NULL)
	{
	  if(strlen(buff)==162 || strlen(buff)==165)
	    {
	      //      printf("%s",buff);
	      prof_col=0;
	      strncpy_NULL(temp,&buff[6],1);
	      strcat(seq,temp);
	      //fractions
	      for(i=71;i<150;i=i+4)
		{
		  
		  strncpy_NULL(temp,&buff[i],3);
		  //printf("%3d ",atoi(temp));
		  //temp2=atoi(temp);
		  //printf("%d %d: %d %d\n",seq_i,prof_col,i,temp2);
		  prof[seq_i][prof_col]=(double)atoi(temp)/100; //(double)temp2/100;
		  //printf("%d %d: %d\n",seq_i,prof_col,i);
		  prof_col++;
		}
	      //the last two columns
	      if(strlen(buff)==162)
		{
		  strncpy_NULL(temp,&buff[152],4);
		  //temp2=atof(temp)*100;
		  prof[seq_i][20]=atof(temp); //temp2; //atof(temp);
		  //printf("%4.2lf ",atof(temp));
		  strncpy_NULL(temp,&buff[157],4);
		  //printf("%4.2lf\n",atof(temp));
		  // temp2=atof(temp)*100;
		  prof[seq_i][21]=atof(temp); //temp2; //atof(temp);
		}
	      seq_i++;
	    }
	  
	}
      // exit(1);
      seq[seq_i]='\0';
      
      //printf("SEQ: %s\n",seq);
      //  printf("SEQ: %s\n",seq);
    }
  else
    {
      //return 0;
      fprintf(stderr,"Could not open %s\n",filename);
      exit(1);
    }

}


char*  read_stride(char* filename)
{
  char  *ss;
  char  *temp;
  char	buff[512];	/* Input string */
  char	line_flag[11];
  char  seq[2000];
  char  str[200];
  int   i=0;
  FILE *fp;
  
  ss=malloc(2000*sizeof(char));
  strcpy(ss,"");
  
  //printf("len temp: %d\n",strlen(temp));
  //printf("len ss: %d\n",strlen(ss));
  fp=fopen(filename,"r");	/* Does file exist? */
  if (fp!=NULL)
    {
      while(fgets(buff,255,fp)!=NULL)
	{
	  sscanf(buff,"%s",line_flag);
	  if(strcmp(line_flag,"SEQ")==0)
	    {
	      //printf("%s\n",buff);
	      strncpy(str,&buff[10],50);
	      

	      //strncpy(ss[i],&buff[24],1);

	      //printf("%d\n",strlen(ss));
	      //printf("%s",buff);
	      //printf("%s\n",str);
	      strcat(seq,str);

	    }
	  if(strcmp(line_flag,"STR")==0)
	    {
	      // printf("%s\n",buff);
	      strncpy(str,&buff[10],50);
	      strcat(ss,str);
	    }
	  if(strcmp(line_flag,"LOC")==0)
	    break;
	}
	  
      strcat(seq,"\0");
      for(i=0;i<strlen(seq);i++){
	if(seq[i] != ' ')
	  {
	    //	    printf("%c ",ss[i]);
	    if(ss[i] == 'H' || ss[i] == 'G' || ss[i] == 'I')
	      {
		ss[i]='H';
	      }
	    else{
	      if(ss[i]=='E')
		{
		  ss[i]='E';
		}
	      else
		{
		  ss[i]='C';
		}
	    }
	    //printf("%c %c\n",seq[i],ss[i]);
	  }
	//printf("%c\n",seq[i]);
	else
	  {
	    break;
	  }
      }
      
      //printf("%d %d\n",i,strlen(seq));
      //printf("%s",seq);
      ss[i]='\0';
      return ss;
    }
  else
    {
      fprintf(stderr,"Could not open %s\n",filename);
      //return NULL;
      exit(1);
      //return 0;
      
    }
}

double* calculate_parameters(char* pdbfile,char* psipred)
{
  molecule      m[1];
  int           i,j,k;
  int           atoms=0;
  int           residues=0;
  int           atom_i=0;
  int           atom_j=0;
  int           exposed_cut;
  int           buried_cut;
  int           total_exposed=0;
  int           total_buried=0;
  int           ss_correct=0;
  int           tot_res_contacts=0;
  int           temp=0;
  int           res_contacts[6][6]={{0}};
  int           contacts[10000]={0};
  int           exposed_residues[6]={0};
  int           buried_residues[6]={0};
  double        dist=0;
  double        mean=0;
  double        fat=0;
  double        cutoff=196;
  double        hcut=3.6;
  double        hangle=1.2;
  double        Q3=0;
  double        frac_res_contacts[21]={0};   // 6*7/2=21
  double        frac_exposed_residues[6]={0};
  double        frac_buried_residues[6]={0};
  double        *parameters;
  //double        *parameters;
  FILE          *fp;
  char          *ss;
  char	        psipredfiles[2000];
  
  
  //printf("%s %s\n",pdbfile,psipred);
  strcpy(m[0].filename,pdbfile);
  parameters=(double*)malloc(35*sizeof(double));

  //TETRA if(read_molecules_backbone(m)==0)
  if(read_molecules(m,'b')==0)
    {  
      //printf("%d\n",m[0].residues);
      for(i=0;i<m[0].residues;i++)  
	{
	  for(j=i+1;j<m[0].residues;j++)
	    {
	      atom_i=m[0].CA_ref[i];
	      atom_j=m[0].CA_ref[j];
	      dist=distance(m,atom_i,atom_j);
	      // printf("%lf\n ",dist);
	      if(dist<100)
		{
		  //printf("%lf\n ",dist);
		  contacts[m[0].atm[atom_i].rescount-1]++;
		  contacts[m[0].atm[atom_j].rescount-1]++;
		}
	      if(abs(m[0].atm[atom_i].rescount-m[0].atm[atom_j].rescount)>5 &&
		 dist<cutoff)
		{
		  res_contacts[get_res6(m[0].atm[atom_i].residue)][get_res6(m[0].atm[atom_j].residue)]++;
		  res_contacts[get_res6(m[0].atm[atom_j].residue)][get_res6(m[0].atm[atom_i].residue)]++;
		  tot_res_contacts=tot_res_contacts+2;
		}
	    }
	}
      for(i=0;i<m[0].residues;i++)
	{
	  if(contacts[i]<16)  //exposed
	    {
	      exposed_residues[get_res6(m[0].atm[m[0].CA_ref[i]].residue)]++;
	      total_exposed++;
			      
	    }
	  if(contacts[i]>20)  //buried
	    {
	      
	      buried_residues[get_res6(m[0].atm[m[0].CA_ref[i]].residue)]++;
	      total_buried++;
	    }
	}
      //printf("Exposed Buried\n");
      for(i=0;i<6;i++)
	{
	  //printf("%d/%d %d/%d\n",exposed_residues[i],total_exposed,buried_residues[i],total_buried);
	  if(total_exposed!=0)
	    frac_exposed_residues[i]=(double)exposed_residues[i]/total_exposed;
	  if(total_buried!=0)
	    frac_buried_residues[i]=(double)buried_residues[i]/total_buried;
	  //printf("%5.4f %5.4f\n",frac_exposed_residues[i],frac_buried_residues[i]);
	}
      //printf("\n");
      for(i=0;i<6;i++)
	{
	  for(j=0;j<6;j++)
	    {
	      //  printf("%d ",res_contacts[i][j]);
	      if(i!=j)
		{
		  res_contacts[i][j]=2*res_contacts[i][j];
		}
	    }
	  //printf("\n");
	}
      for(k=0,i=0;i<6;i++)
	{ 
	  for(j=i;j<6;j++,k++) 
	    {
	      if(tot_res_contacts != 0)
		{
		  frac_res_contacts[k]=(double)res_contacts[i][j]/tot_res_contacts;
		  //printf("%f %d %d %d\n",frac_res_contacts[k],res_contacts[i][j],tot_res_contacts,&tot_res_contacts);
		}
	    }
	}
    }
  //printf("\n");
  fat=fatness(m);
  
  //FIX THAT SS=0 MEANS ERROR
  ss=assign_ss(m,hcut,hangle);
  // if(ss!=0)
  //  {
      if(strlen(psipred)==strlen(ss))
	{
	  for(i=0;i<strlen(psipred);i++)
	    {
	      //printf("%d %d\n",psipred[i],ss[i]);
	      //printf("%d\n",(psipred[i]==ss[i]));
	      if(psipred[i]==ss[i])
		{
		  ss_correct++;
		}
	    }
	  //printf("Q3: %d %5.3lf\n",ss_correct,(double)ss_correct/strlen(psipred));
	  Q3=(double)ss_correct/strlen(psipred);
	  
	  for(j=0,i=0;i<21;i++,j++)
	    {
	      parameters[j]=frac_res_contacts[i];
	      //printf("%f ",frac_res_contacts[i]);
	      
	    }
	  for(i=0;i<6;i++,j++)
	    {
	      parameters[j]=frac_exposed_residues[i];
	      //printf("%f ",frac_exposed_residues[i]);
	      
	    }
	  for(i=0;i<6;i++,j++)
	    {
	      parameters[j]=frac_buried_residues[i];
	      //  printf("%f ",frac_buried_residues[i]);
	    }
	  parameters[j]=fat;
	  parameters[j+1]=Q3;
	  //printf("%f %f %d",fat,Q3,j);
	  //printf("\n");
	  //for(i=0;i<35;i++)
	  //	    printf("%lf ",parameters[i]);
	  //printf("\n");
	}
      else
	{
	  printf("Different strlen!!!!fgfdgfdg %lu %lu\n%s\%s\n",strlen(psipred),strlen(ss),psipred,ss);
	  
	}
    
  free(ss);
  return parameters;
}
 
double* ProQCA(char* pdbfile) //,char* psipred)
{
  molecule      m[1];
  network       net_lg[5];
  network       net_mx[5];
  int           i,j,k;
  int           atoms=0;
  int           residues=0;
  int           atom_i=0;
  int           atom_j=0;
  int           exposed_cut;
  int           buried_cut;
  int           total_exposed=0;
  int           total_buried=0;
  int           ss_correct=0;
  int           tot_res_contacts=0;
  int           temp=0;
  int           res_contacts[6][6]={{0}};
  int           contacts[10000]={0};
  int           exposed_residues[6]={0};
  int           buried_residues[6]={0};
  double        dist=0;
  double        mean=0;
  double        fat=0;
  double        cutoff=196;
  double        hcut=3.6;
  double        hangle=1.2;
  double        Q3=0;
  double        pred_lg=0;
  double        pred_mx=0;
  double        frac_res_contacts[21]={0};   // 6*7/2=21
  double        frac_exposed_residues[6]={0};
  double        frac_buried_residues[6]={0};
  double        parameters[35];
  //double        *parameters;
  double        *quality; //(LG,MX)
  FILE          *fp;
  char          *ss;
  //char	        psipredfiles[2000];
  char          lg1[200],lg2[200],lg3[200],lg4[200],lg5[200],mx1[200],mx2[200],mx3[200],mx4[200],mx5[200];

  
  //printf("%s %s\n",pdbfile,psipred);
  strcpy(m[0].filename,pdbfile);
  quality=(double*)malloc(2*sizeof(double));
  quality[0]=0;
  quality[1]=0;
  
  if(getenv("PROQDIR")==NULL)
    {
      fprintf(stderr,"ENV variable PROQDIR needs to be set properly\n");
      exit(1);
    }
  strcpy(lg1,getenv("PROQDIR"));
  strcpy(lg2,lg1);
  strcpy(lg3,lg1);
  strcpy(lg4,lg1);
  strcpy(lg5,lg1);
  strcpy(mx1,lg1);
  strcpy(mx2,lg1);
  strcpy(mx3,lg1);
  strcpy(mx4,lg1);
  strcpy(mx5,lg1);

  read_net((char*)strcat(lg1,"lg1"),&net_lg[0]);
  read_net((char*)strcat(lg2,"lg2"),&net_lg[1]);
  read_net((char*)strcat(lg3,"lg3"),&net_lg[2]);
  read_net((char*)strcat(lg4,"lg4"),&net_lg[3]);
  read_net((char*)strcat(lg5,"lg5"),&net_lg[4]);

  read_net((char*)strcat(mx1,"mx1"),&net_mx[0]);
  read_net((char*)strcat(mx2,"mx2"),&net_mx[1]);
  read_net((char*)strcat(mx3,"mx3"),&net_mx[2]);
  read_net((char*)strcat(mx4,"mx4"),&net_mx[3]);
  read_net((char*)strcat(mx5,"mx5"),&net_mx[4]);


  //TETRA if(read_molecules_backbone(m)==0)
  if(read_molecules(m,'b')==0)
    {  
      //printf("%d\n",m[0].residues);
      if(strlen(m[0].sequence)<19)
	{
	  return quality;
	}
      if(m[0].residues==m[0].atoms)
	{
	  quality[0]=-2;
	  quality[1]=-2;
	  return quality;
	}
      for(i=0;i<m[0].residues;i++)  
	{
	  for(j=i+1;j<m[0].residues;j++)
	    {
	      atom_i=m[0].CA_ref[i];
	      atom_j=m[0].CA_ref[j];
	      dist=distance(m,atom_i,atom_j);
	      // printf("%lf\n ",dist);
	      if(dist<100)
		{
		  //printf("%lf\n ",dist);
		  contacts[m[0].atm[atom_i].rescount-1]++;
		  contacts[m[0].atm[atom_j].rescount-1]++;
		}
	      if(abs(m[0].atm[atom_i].rescount-m[0].atm[atom_j].rescount)>5 &&
		 dist<cutoff)
		{
		  res_contacts[get_res6(m[0].atm[atom_i].residue)][get_res6(m[0].atm[atom_j].residue)]++;
		  res_contacts[get_res6(m[0].atm[atom_j].residue)][get_res6(m[0].atm[atom_i].residue)]++;
		  tot_res_contacts=tot_res_contacts+2;
		}
	    }
	}
      for(i=0;i<m[0].residues;i++)
	{
	  if(contacts[i]<16)  //exposed
	    {
	      exposed_residues[get_res6(m[0].atm[m[0].CA_ref[i]].residue)]++;
	      total_exposed++;
			      
	    }
	  if(contacts[i]>20)  //buried
	    {
	      
	      buried_residues[get_res6(m[0].atm[m[0].CA_ref[i]].residue)]++;
	      total_buried++;
	    }
	}
      //printf("Exposed Buried\n");
      for(i=0;i<6;i++)
	{
	  //printf("%d/%d %d/%d\n",exposed_residues[i],total_exposed,buried_residues[i],total_buried);
	  if(total_exposed!=0)
	    frac_exposed_residues[i]=(double)exposed_residues[i]/total_exposed;
	  if(total_buried!=0)
	    frac_buried_residues[i]=(double)buried_residues[i]/total_buried;
	  //printf("%5.4f %5.4f\n",frac_exposed_residues[i],frac_buried_residues[i]);
	}
      //printf("\n");
      for(i=0;i<6;i++)
	{
	  for(j=0;j<6;j++)
	    {
	      //  printf("%d ",res_contacts[i][j]);
	      if(i!=j)
		{
		  res_contacts[i][j]=2*res_contacts[i][j];
		}
	    }
	  //printf("\n");
	}
      for(k=0,i=0;i<6;i++)
	{ 
	  for(j=i;j<6;j++,k++) 
	    {
	      if(tot_res_contacts != 0)
		{
		  frac_res_contacts[k]=(double)res_contacts[i][j]/tot_res_contacts;
		  //printf("%f %d %d %d\n",frac_res_contacts[k],res_contacts[i][j],tot_res_contacts,&tot_res_contacts);
		}
	    }
	}
    }
  //printf("\n");
  fat=fatness(m);
  //printf("SS\n");
  ss=assign_ss(m,hcut,hangle);
  //printf("AFTER SS\n");
  //printf("%s\nSTR:  %s\n",m[0].filename,m[0].ss);
  if(strlen(m[0].ss)==strlen(ss))
    {
      //  printf("PRED: %s\n",ss);
      for(i=0;i<strlen(m[0].ss);i++)
	{
	  //printf("%d %d\n",m[0].ss[i],ss[i]);
	  //printf("%d\n",(m[0].ss[i]==ss[i]));
	  if(m[0].ss[i]==ss[i])
	    {
	      ss_correct++;
	    }
	}
      //printf("Q3: %d %5.3lf\n",ss_correct,(double)ss_correct/strlen(m[0].ss));
      Q3=(double)ss_correct/strlen(m[0].ss);
      //Q3=; //0.75;
      for(j=0,i=0;i<21;i++,j++)
	{
	  parameters[j]=frac_res_contacts[i];
	  //printf("%f ",frac_res_contacts[i]);
	  
	}
      for(i=0;i<6;i++,j++)
	{
	  parameters[j]=frac_exposed_residues[i];
	  //printf("%f ",frac_exposed_residues[i]);
	  
	}
      for(i=0;i<6;i++,j++)
	{
	  parameters[j]=frac_buried_residues[i];
	  //  printf("%f ",frac_buried_residues[i]);
	}
      //if(fat > 8)
      //	{
      //	  fprintf(stderr,"%s Fatness (%lf) > 8 assigned fatness 8 to avoid strange results\n",m[0].filename,fat);
      //	  fat=8;
      //	}
      parameters[j]=fat;
      parameters[j+1]=Q3;
      //printf("%f %f %d",fat,Q3,j);
      //printf("\n");
      //for(i=0;i<35;i++)
      //	printf("%lf ",parameters[i]);
      //printf("\n");
      
      
    }
  else
    {
      if(strlen(m[0].ss)!=strlen(ss))
	fprintf(stderr,"%s Different strlen\nPDB-SS %s (%lu)\nASSIGN-SS %s (%lu)\nSEQ %s\n",m[0].filename,m[0].ss,strlen(m[0].ss),ss,strlen(ss),m[0].sequence);
      //if(fat > 8)
      //	fprintf(stderr,"%s Fatness (%lf) > 8 assign quality 0 to avoid strange results\n",m[0].filename,fat);
      free(ss);
      quality[0]=0;
      quality[1]=0;
      return quality;
    }
  free(ss);
  if(fat > 8)
    {
      fprintf(stderr,"%s Fatness (%lf) > 8 assigned quality 0 to avoid strange results\n",m[0].filename,fat);
      quality[0]=0;
      quality[1]=0;
    }
  else
    {
      for(i=0;i<5;i++)
	quality[0]+=netfwd(parameters,&net_lg[i]);
      for(i=0;i<5;i++)
	quality[1]+=netfwd(parameters,&net_mx[i]);
      quality[0]=quality[0]/5;
      quality[1]=quality[1]/5;
    }
  //printf("%-30s %6.4lf %6.4lf\n",pdbfile,quality[0]/5,quality[1]/5);
  return quality;
}
 
void ProQ(char** atom_vec,char** residue_vec,int* number_vec,rvec* coord_vec,size_t n,char* psipred,double* quality)
{
  molecule      m[1];
  network       net_lg[5];
  network       net_mx[5];
  int           i,j,k;
  int           atoms=0;
  int           residues=0;
  int           atom_i=0;
  int           atom_j=0;
  int           exposed_cut;
  int           buried_cut;
  int           total_exposed=0;
  int           total_buried=0;
  int           ss_correct=0;
  int           tot_res_contacts=0;
  int           temp=0;
  int           res_contacts[6][6]={{0}};
  int           contacts[10000]={0};
  int           exposed_residues[6]={0};
  int           buried_residues[6]={0};
  double        dist=0;
  double        mean=0;
  double        fat=0;
  double        cutoff=196;
  double        hcut=3.6;
  double        hangle=1.2;
  double        Q3=0;
  double        pred_lg=0;
  double        pred_mx=0;
  double        frac_res_contacts[21]={0};   // 6*7/2=21
  double        frac_exposed_residues[6]={0};
  double        frac_buried_residues[6]={0};
  double        parameters[35];
  //double        *parameters;
  //double        *quality; //(LG,MX)
  FILE          *fp;
  char          *ss;
  char	        psipredfiles[2000];
  char          lg1[200],lg2[200],lg3[200],lg4[200],lg5[200],mx1[200],mx2[200],mx3[200],mx4[200],mx5[200];


  //printf("In ProQ\n");
  //printf("%s %s\n",pdbfile,psipred);
  //strcpy(m[0].filename,pdbfile);
  //quality=(double*)malloc(2*sizeof(double));
  quality[0]=0;
  quality[1]=0;
  

  strcpy(lg1,getenv("PROQDIR"));
  strcpy(lg2,lg1);
  strcpy(lg3,lg1);
  strcpy(lg4,lg1);
  strcpy(lg5,lg1);
  strcpy(mx1,lg1);
  strcpy(mx2,lg1);
  strcpy(mx3,lg1);
  strcpy(mx4,lg1);
  strcpy(mx5,lg1);

  read_net((char*)strcat(lg1,"lg1"),&net_lg[0]);
  read_net((char*)strcat(lg2,"lg2"),&net_lg[1]);
  read_net((char*)strcat(lg3,"lg3"),&net_lg[2]);
  read_net((char*)strcat(lg4,"lg4"),&net_lg[3]);
  read_net((char*)strcat(lg5,"lg5"),&net_lg[4]);

  read_net((char*)strcat(mx1,"mx1"),&net_mx[0]);
  read_net((char*)strcat(mx2,"mx2"),&net_mx[1]);
  read_net((char*)strcat(mx3,"mx3"),&net_mx[2]);
  read_net((char*)strcat(mx4,"mx4"),&net_mx[3]);
  read_net((char*)strcat(mx5,"mx5"),&net_mx[4]);
  //printf("In ProQ\n");
  //printf("%lf %lf %lf\n",coord_vec[0][0],coord_vec[0][1],coord_vec[0][2]);
  if(read_to_molecule(m,atom_vec,residue_vec,number_vec,coord_vec,n)==0 && strlen(m[0].sequence)>18)
    {  
      //printf("%d\n",m[0].residues);
      //if(strlen(m[0].sequence)<19)
      //	{
      //	  return quality;
      //	}
      for(i=0;i<m[0].residues;i++)  
	{
	  for(j=i+1;j<m[0].residues;j++)
	    {
	      atom_i=m[0].CA_ref[i];
	      atom_j=m[0].CA_ref[j];
	      dist=distance(m,atom_i,atom_j);
	      //printf("%d %d %lf\n",atom_i,atom_j,dist);
	      if(dist<100)
		{
		  //printf("%lf\n ",dist);
		  contacts[m[0].atm[atom_i].rescount-1]++;
		  contacts[m[0].atm[atom_j].rescount-1]++;
		}
	      if(abs(m[0].atm[atom_i].rescount-m[0].atm[atom_j].rescount)>5 &&
		 dist<cutoff)
		{
		  res_contacts[get_res6(m[0].atm[atom_i].residue)][get_res6(m[0].atm[atom_j].residue)]++;
		  res_contacts[get_res6(m[0].atm[atom_j].residue)][get_res6(m[0].atm[atom_i].residue)]++;
		  tot_res_contacts=tot_res_contacts+2;
		}
	    }
	}
      for(i=0;i<m[0].residues;i++)
	{
	  //printf("%c contacts: %d\n",m[0].sequence[i],contacts[i]);
	  if(contacts[i]<16)  //exposed
	    {
	      exposed_residues[get_res6(m[0].atm[m[0].CA_ref[i]].residue)]++;
	      total_exposed++;
	    }
	  if(contacts[i]>20)  //buried
	    {
	      
	      buried_residues[get_res6(m[0].atm[m[0].CA_ref[i]].residue)]++;
	      total_buried++;
	    }
	}
      //printf("Exposed Buried\n");
      for(i=0;i<6;i++)
	{
	  if(total_exposed!=0)
	    frac_exposed_residues[i]=(double)exposed_residues[i]/total_exposed;
	  if(total_buried!=0)
	    frac_buried_residues[i]=(double)buried_residues[i]/total_buried;
	  //printf("%d/%d=%5.4f %d/%d=%5.4f\n",exposed_residues[i],total_exposed,frac_exposed_residues[i],buried_residues[i],total_buried,frac_buried_residues[i]);
	}
      //printf("\n");
      for(i=0;i<6;i++)
	{
	  for(j=0;j<6;j++)
	    {
	      //  printf("%d ",res_contacts[i][j]);
	      if(i!=j)
		{
		  res_contacts[i][j]=2*res_contacts[i][j];
		}
	    }
	  //printf("\n");
	}
      for(k=0,i=0;i<6;i++)
	{ 
	  for(j=i;j<6;j++,k++) 
	    {
	      if(tot_res_contacts != 0)
		{
		  frac_res_contacts[k]=(double)res_contacts[i][j]/tot_res_contacts;
		  //printf("%f %d %d %d\n",frac_res_contacts[k],res_contacts[i][j],tot_res_contacts,&tot_res_contacts);
		}
	    }
	}
    }
  //printf("\n");
  fat=fatness(m);
  ss=assign_ss(m,hcut,hangle);
  //  printf("PRED: %s\n  SS: %s\n",psipred,ss);
  if(strlen(psipred)==strlen(ss))
    {
      for(i=0;i<strlen(psipred);i++)
      {
	//printf("%d %d\n",psipred[i],ss[i]);
	//printf("%d\n",(psipred[i]==ss[i]));
	if(psipred[i]==ss[i])
	  {
	    ss_correct++;
	  }
      }
      //printf("Q3: %d %5.3lf\n",ss_correct,(double)ss_correct/strlen(psipred));
      Q3=(double)ss_correct/strlen(psipred);
      //if(Q3<0.4)
      //	{
      //	  Q3=0;
      //	}
      //Q3=; //0.75;
      for(j=0,i=0;i<21;i++,j++)
	{
	  parameters[j]=frac_res_contacts[i];
	  //printf("%f ",frac_res_contacts[i]);
	  
	}
      for(i=0;i<6;i++,j++)
	{
	  parameters[j]=frac_exposed_residues[i];
	  //printf("%f ",frac_exposed_residues[i]);
	  
	}
      for(i=0;i<6;i++,j++)
	{
	  parameters[j]=frac_buried_residues[i];
	  //  printf("%f ",frac_buried_residues[i]);
	}
      if(fat > 8)
	{
	  fprintf(stderr,"%s Fatness (%lf) > 8 assigned fatness 8 to avoid strange results\n",m[0].filename,fat);
	  fat=8;
	}
      parameters[j]=fat;
      
      parameters[j+1]=Q3;
      //printf("%f %f %d",fat,Q3,j);
      //printf("\n");
      //      for(i=0;i<35;i++)
      //        printf("%lf ",parameters[i]);
      // printf("\n");


    }
     else
    {
      fprintf(stderr,"Different strlen %lu %lu\nPRED: %s\n  SS: %s\n SEQ: %s\n",strlen(psipred),strlen(ss),psipred,ss,m[0].sequence);
      free(ss);
      quality[0]=-1;
      quality[1]=-1;
      //return quality;
    }
  free(ss);
  for(i=0;i<5;i++)
    quality[0]+=netfwd(parameters,&net_lg[i]);
  for(i=0;i<5;i++)
    quality[1]+=netfwd(parameters,&net_mx[i]);
  quality[0]=quality[0]/5;
  quality[1]=quality[1]/5;
  
  //printf("LG: %6.4lf MX: %6.4lf\n",quality[0],quality[1]);
  //printf("%-30s %6.4lf %6.4lf\n",pdbfile,quality[0]/5,quality[1]/5);
  //return quality;
}

