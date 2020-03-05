/*
  Program LGSCORE   written by Scott M. Le Grand
  Even more modified and completely changed in 1999 by Arne Elofsson 
  Modified by Arne Elofsson in 1994
  Copyright 1992, Scott M. Le Grand and Arne Elofsson (1999)
  The program is available under the Gnu public license (see www.fsf.org)

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License. With the
  exception that if you use this program in any scientific work you have
  to explicitly state that you have used LGSCORE and cite the relevant
  publication (dependent on what you have used LGSCORE for). My
  publicationlist can be found at http://www.sbc.su.se/~arne/papers/.
  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details. You should have received a
  copy of the GNU General Public License along with this program (in 
  the file gpl.txt); if not, write to the Free Software Foundation, 
  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-13
*/


#include <stdio.h>
#include <math.h> 
#include <stdlib.h>
#include <string.h>
//#include "/afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/molecule.h"
//#include "/afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/quality_measure/lgscore.h"
#include "molecule.h"
#include "lgscore.h"


void Dscore_res(char* file1,char* file2,lgscore *LG, double d0, double minsim,int L,double factor) //,int step)
{

  
  molecule	m[3];		/* Molecules to be compared */
  int		out_flag;		
  int		b_flag;		
  int		align_flag;		
  int		select1_flag;		
  int		select2_flag;		
  int		verbose_flag;		/* Verbose output flag */
  int             fraction_flag;            /* Output flag */
  int             short_flag;
  int             arglist;
  char            outfile[800];            /* Output name */
  double		rms;			/* final RMS error */
  int 		ma[2];
  /*int             temp;*/


  //molecule      m0,m1;
  int	files=0;	/* Number of files parsed */
  int	i,j,k,a,b;	/* Counter variable */
  //  int	L=4;	// Minimum length of fragment for FindMaxSub matching
  int    minatoms=19; /* Variable for min no of atoms for P-value*/
  double	s[3][3];	/* Final transformation matrix */
  double          maxrmsd,minrmsd;
  double          fraction=1;
  double          numfrac=0;
  double          maxpvalue=1,maxpvalue1=1.;
  double          score,maxrms,maxrms1,maxscore,maxscore1;
  double          pvalue;
  double          logDscore=0;
  double          temp=0,temp2=0;
  
  int             maxatoms,maxpos,minpos,atoms;
  int             length=0;
  double          bestsizescore[MAXATMS];
  int             maxatoms1;
  int             all_read=1;
  double          cutoffrmsd=0.;
  //  double          minsim=25.0;
  //  double          cutoffrms=3.;
  double          cutoff;
  double ss;//Added by Fang
  //double          *Sstr;
  //double          *Sstr_return;
  int             count;
  int             maxcount=10;
  char            filename[400][800];
  char            target[400][10];
  char            template[400][10];
  char            method[400][20];
  int             rank[400];
  char            filelist[800];
  char            bb[3];
  int             len_filelist;
  int             dotcount=0;
  int             lastchar=' ';
  FILE            *fp;  
  double          Ssum=0;
  double          TMsum=0;
  double          d0_TM=5; 
  

  /* Initialize */
  files=2;
  verbose_flag=FALSE;
  out_flag=FALSE;
  select2_flag=FALSE;
  select1_flag=FALSE;
  short_flag=FALSE;
  arglist=FALSE;
  fraction_flag=FALSE; //TRUE; //FALSE;
  b_flag=FALSE;
  align_flag=FALSE;
  

  /* Parse command line for PDB files and options */
  i=1;
  
  
  
  //  for(a=0;a<files;a++)
  // {
  //   for(b=a+1;b<files;b++)
  //	{
  //printf("%d %d\n",a,b);
  //if(strcmp(method[a],method[b]))
  //{
  maxpvalue=1;
  maxpvalue1=1.;
  strcpy(m[0].filename,file1);
  strcpy(m[1].filename,file2);
  //strcpy(m[2].filename,file1);
  //printf("%-40s %-40s\t",file1,file2);
  //TETRA if(read_molecules_ca(&m[0])==0 && read_molecules_ca(&m[1])==0 && read_molecules_ca(&m[2])==0)
  if(read_molecules(&m[0],'c')==0 && read_molecules(&m[1],'c')==0)
    {
      length=max(m[0].atoms,m[1].atoms);
      check_molecules(&m[0],&m[1]); // this deletes some atoms
      //printf("testing %d %d\n",m[0].atoms,m[1].atoms);
      //atoms=m[0].atoms;
      
      for (i=0;i<m[0].atoms;i++){bestsizescore[i]=0;}
      for (i=0;i<m[0].atoms;i++){bestsizescore[i]=0;
	for (j=0;j<m[0].atoms;j++){
	  
	}
      }
      score=
      }
    }
  else
    {
      printf("Unable to do what you told me... My fault... I'm stupied\n");
    }
  //}
  //}
      
  //}

  
  
  
  for (i=0;i<m[0].atoms;i++){m[0].atm[i].selected=maxselect1_0[i];}
  for (i=0;i<m[1].atoms;i++){m[1].atm[i].selected=maxselect1_1[i];}
  center_molecule(&m[0]);
  center_molecule(&m[1]);
  rms=superimpose_molecules(&m[0],&m[1],s,0.0001);  //stricter error cut for the last superimpose.

    
  //Sstr=(double*)malloc(sizeof(double)*(m[2].atoms+1));
  //Sstr_return=(double*)malloc(sizeof(double)*(m[2].atm[m[2].atoms-1].resnum+1));
  for(i=0;i<m[2].atm[m[2].atoms-1].resnum+1;i++)
    {
      //printf("%d %d %d %d %d\n",i,m[2].atoms,m[2].atm[0].resnum,m[2].atm[m[2].atoms-1].resnum,m[2].atoms+m[2].atm[0].resnum);
      D[0].S[i]=0;
      D[0].TM[i]=0;
    }
  D[0].S[i-1]=-1;
  D[0].TM[i-1]=-1;
  //exit(1);
  
  //  d0=25;
  //Sstr[0]=m[2].atoms+1;
  D[0].residues=0;
  for (i=0,j=0;i<m[2].atoms;i++)
    {
      
      if(m[2].atm[i].resnum == m[0].atm[j].resnum) // && m[0].atm[j].selected)
	{
	  if(isnan(m[0].atm[j].rms))
	    {
	      D[0].S[i]=0.0;
	      D[0].TM[i]=0.0;
	      // if(m[0].atm[j].selected)
	      //  temp2+=1/(1+m[0].atm[j].rms/d0);
	      m[2].atm[i].rms=998001;
	      m[0].atm[j].rms=998001;
	    }
	  else
	    {
	      //printf("%lf\n",d0);
	      D[0].S[i]=1/(1+m[0].atm[j].rms/(d0*d0));
	      D[0].TM[i]=1/(1+m[0].atm[j].rms/(d0_TM*d0_TM));
	      // if(m[0].atm[j].selected)
	      //  temp2+=1/(1+m[0].atm[j].rms/d0);
	      m[2].atm[i].rms=m[0].atm[j].rms;
	    }
	  j++;
	  
	}
      else
	{
	  D[0].S[i]=0;
	  D[0].TM[i]=0;
	  m[2].atm[i].rms=998001;
	  D[0].S[i]=1/(1+998001/(d0*d0));
	  D[0].TM[i]=1/(1+998001/(d0_TM*d0_TM));
	}
      //Sstr[i]=1/(1+m[2].atm[i].rms/5);
      //if(m[2].atm[i].selected)
      //	temp2+=1/(1+m[2].atm[i].rms/5);
      //Sstr_return[m[2].atm[i].resnum-1]=Sstr[i];
      Ssum=Ssum+D[0].S[i];
      TMsum=TMsum+D[0].TM[i];
      
      D[0].residue[i]=aa321(m[2].atm[i].residue);
      D[0].resnum[i]=m[2].atm[i].resnum;
      D[0].residues++;
      // printf("%d %c %lf %lf\n",m[2].atm[i].resnum,aa321(m[2].atm[i].residue),D[0].S[i],D[0].TM[i]); //,m[0].atm[j-1].rms);
       //Sstr_return[m[2].atm[i].resnum]
    }

  
  //printf("SCORE:\t%d\t%d\t%.1f\t%.1f\t%e\t%e\t%7.5lf\n",length,maxatoms1,maxrms1,maxscore1,maxpvalue1,ss,D[0].Dscore);
 

  // free(Sstr);
  //return Sstr_return;
  //exit(0);
  
}
  
  


void Dscore_res_pt(dyn_molecule *m1,dyn_molecule *m2,lgscore *D, double d0, double minsim,int L,double factor,int step)
{
  
  //molecule	m[3];		/* Molecules to be compared */
    molecule *m;
    int		out_flag;
    //    int superimpose_all;
    int		b_flag;		
    int		align_flag;		
    int		select1_flag;		
    int		select2_flag;		
    int		verbose_flag;		/* Verbose output flag */
    int             fraction_flag;            /* Output flag */
    int             short_flag;
    int             arglist;
    char            outfile[800];            /* Output name */
    double		rms;			/* final RMS error */
    int 		ma[2];
    int	files=0;	/* Number of files parsed */
    int	i,j,k,l,a,b;	/* Counter variable */
    int    minatoms=19; /* Variable for min no of atoms for P-value*/
    double	s[3][3];	/* Final transformation matrix */
    double	best_s[3][3];	/* Final transformation matrix */
    double          maxrmsd,minrmsd;
    double          fraction=1;
    double          numfrac=0;
    double          maxpvalue=1,maxpvalue1=1.;
    double          score,maxrms,maxrms1,maxscore,maxscore1;
    double          pvalue;
    double          logDscore=0;
    double          temp=0,temp2=0;
    int             maxatoms,maxpos,minpos,atoms;
    int             length=0;
    int             maxsub_15=0;
    int             maxsub_20=0;
    int             maxsub_25=0;
    int             maxsub_30=0;
    int             maxsub_35=0;
    int             maxselect1_0[MAXATMS],maxselect1_1[MAXATMS];
    int             maxselect2_0[MAXATMS],maxselect2_1[MAXATMS];
    double          bestsizescore[MAXATMS];
    int             maxatoms1;
    int             all_read=1;
    double          cutoffrmsd=0.;
    double          cutoff;
    double ss;//Added by Fang
    int             count;
    int             maxcount=10;
    char            filename[400][800];
    char            target[400][10];
    char            template[400][10];
    char            method[400][20];
    int             rank[400];
    char            filelist[800];
    char            bb[3];
    int             len_filelist;
    int             dotcount=0;
    int             lastchar=' ';
    FILE            *fp;  
    double          Ssum=0;
    double          TMsum=0;
    double          d0_TM=5; 
    /* Initialize */
    files=2;
    b_flag=FALSE;
    verbose_flag=FALSE;
    out_flag=FALSE;
    select2_flag=FALSE;
    select1_flag=FALSE;
    short_flag=FALSE;
    arglist=FALSE;
    fraction_flag=FALSE; //TRUE; //FALSE;
    align_flag=FALSE;

    /* Parse command line for PDB files and options */
    i=1;
  
    
    //  for(a=0;a<files;a++)
    // {
    //   for(b=a+1;b<files;b++)
    //	{
    //printf("%d %d\n",a,b);
    //if(strcmp(method[a],method[b]))
    //{
    maxpvalue=1;
  
    
    //printf("Hello1\n");
    //return;
    m=malloc(3*sizeof(molecule));
    // BW: this is stupied... but the old code deletes everything that 
    // is not the same in both, so before that is rewritten it has to 
    // be copied to keep the original...
    copymolecule(&m[0],m1); //copies only the CA...
    copymolecule(&m[1],m2);
    copymolecule(&m[2],m1);

   
  
  
    //printf("Hello1\n");
    //return;
 
  
    //printf("Hello2\n");
    //if(arglist)
    //	{
    //	  molecule m[2];
    ///	  copymolecule(&m[0],&m_list[a]);
    //  copymolecule(&m[1],&m_list[b]);
    //printf("%-30s(%d) %-30s(%d) %-30s(%d)\n",m[0].filename,m[0].atoms,m[1].filename,m[1].atoms,m[2].filename,m[2].atoms);
    //	}
    length=max(m[0].atoms,m[1].atoms);
    d0_TM=(double)1.24*pow(length-15,0.3333333333)-1.8;
    //printf("length: %d\nd0_TM: %lf\ntemp: %lf",length,d0_TM,pow(length-15,0.333333333));
    // exit(1);
    check_molecules(&m[0],&m[1]); // this deletes some atoms
    //printf("testing %d %d\n",m[0].atoms,m[1].atoms);
    //atoms=m[0].atoms;
  
      
      for (i=0;i<m[0].atoms;i++){bestsizescore[i]=0;}
      for (i=0;i<m[0].atoms-L;i=i+step){
	//printf("%d %d",i,L);
	for (j=0;j<m[0].atoms;j++){ // deselect all atoms
	  m[0].atm[j].selected=FALSE;
	  m[1].atm[j].selected=FALSE;
	}
	for (j=i;j<i+L;j++){ // select some
	  //printf("%d %d %d\n",j,L,i);
	  m[0].atm[j].selected=TRUE;
	  m[1].atm[j].selected=TRUE;
	}
	j=L;
	center_molecule(&m[0]);
	center_molecule(&m[1]);
	reset_rms(&m[0]);
	reset_rms(&m[1]);
	reset_rms(&m[2]);
	//printf("calling superimpose 1\n");
	rms=superimpose_molecules(&m[0],&m[1],s,0.1);	   
	count=0;
	while (j<m[0].atoms && rms<(double)(factor*(j+225)/225)){
	  count++;
	  // printf("%d %d\n",j,count);
	  minrmsd=9999.0;
	  for (k=0;k<m[0].atoms;k++) {
	    if(!m[0].atm[k].selected)
	      {
		if(m[0].atm[k].rms < minrmsd)
		  {
		    minpos=k;
		    minrmsd=m[0].atm[k].rms;
		  }
		//else 
		if(m[0].atm[k].rms < minsim)
		  {
		    //printf("Selecting %d %f\n",k,m[0].atm[k].rms);
		    m[0].atm[k].selected=TRUE;
		    m[1].atm[k].selected=TRUE;
		    j++;
		  }
	      }
	  }
	  j++;
	  //printf("Selecting %d %f %d\n",minpos,minrmsd,j);
	  m[0].atm[minpos].selected=TRUE;
	  m[1].atm[minpos].selected=TRUE;
	  
	  //	     printf("testing2a %f\n",m[1].atm[i].x);
	  center_molecule(&m[0]);
	  center_molecule(&m[1]);
	  //printf("calling superimpose 2\t %d %d\n",i,j);
	  rms=superimpose_molecules(&m[0],&m[1],s,0.1); //Not so strict error cut of here only 0.1.
#ifdef Sscore
	  score=Levitt_Gerstein(&m[0],&m[1],d0*d0);
#else
	  score=Levitt_Gerstein(&m[0],&m[1],5);
#endif
	  //pvalue=D_pvalue(j,score);
	  pvalue=D_pvalueF(j,score);
	  
	  //printf("TEST: %e\n", D_pvalueF(227,1532.3));
	  if (b_flag)
	    {
	      if (score > bestsizescore[j]){bestsizescore[j]=score;}
	    }
	  if (fraction_flag)
	    {
	      if ((double)j/atoms <= fraction )
		{
		  printf("FRAC1:\t%.2f\t%d\t%.1f\t%.1f\t%e\n",
			 (double)j/atoms,j,rms,score,
			 pvalue);
		  fraction-=numfrac;
		  //printf("TEST %f %f\n",fraction,numfrac);
		} 
	    }
#ifdef Sscore
	  if ((score >= maxscore) && (j>minatoms)) 
#else
	  if ((pvalue <= maxpvalue) && (j>minatoms)) 
#endif
	    {
	      if (fraction_flag)
		{
		  printf("GOOD1:\t%.2f\t%d\t%.1f\t%.1f\t%e\n",
			 (double)j/atoms,j,rms,score,
			 pvalue);
		}
	      //printf ("TEST:\t%d\t%f\t%f\t%f\n",m[0].atoms,pvalue,rms,score);
	      maxpvalue=pvalue;
	      maxatoms=j;
	      maxrms=rms;
	      maxscore=score;
	      for (k=0;k<m[0].atoms;k++)
		{
		  maxselect1_0[k]=m[0].atm[k].selected;
		}
	      for (k=0;k<m[1].atoms;k++)
		{
		  maxselect1_1[k]=m[1].atm[k].selected;
		}

	      for(k=0;k<3;k++)
		{
		  for(l=0;l<3;l++)
		    {
		      //		      printf("%8.5lf ",s[k][l]);
		      best_s[k][l]=s[k][l];
		    }
		  // printf("\n");
		}
	    }
	}
	//	     printf("testing5 %f\n",m[1].atm[i].x);
      }
      //for(k=0;k<3;k++)
      //	{
      //	  for(l=0;l<3;l++)
      //	    {
      //	      //		      printf("%8.5lf ",s[k][l]);
      //	      printf("%8.5lf ",best_s[k][l]);
      //	    }
      //	  printf("\n");
      //	}
      ///printf("%d %d\n",count,j);
      // printf("BEST1:\t%.2f\t%d\t%.1f\t%.1f\t%e\n",
      //	  (double)maxatoms/(double)atoms,maxatoms,maxrms,maxscore,maxpvalue);
      maxatoms1=maxatoms;
      maxrms1=maxrms;
      maxscore1=maxscore;
      maxpvalue1=maxpvalue;
      //(double)maxatoms/(double)m[0].atoms,maxatoms,maxrms,maxscore,maxpvalue,maxpvalue1);
      if (maxatoms1<minatoms){maxatoms1=maxrms1=maxscore1=0;maxpvalue1=1;}
      //           ss=(log10(maxpvalue)+0.6108)/ma1;
	      
      //lg[a][b]=(log10(1/maxpvalue1));
      //lg[b][a]=lg[a][b];
      if(short_flag)
	{
	  printf("%7.5lf\n",(log10(1/maxpvalue1)));
	}
      else
	{
	  ss=-1000*(log10(maxpvalue1)+0.6108)/m[0].atoms; //Why m[0]????? dependent on order
	  //printf("SCORE:\t%d\t%d\t%.1f\t%.1f\t%e\t%e\t%7.5lf\n",
	  // length,maxatoms1,maxrms1,maxscore1,maxpvalue1,ss,(log10(1/maxpvalue1)));
	  logDscore=(log10(1/maxpvalue1));
	}
      if(align_flag){
	for (i=0;i<m[0].atoms;i++){m[0].atm[i].selected=maxselect1_0[i];}
	for (i=0;i<m[1].atoms;i++){m[1].atm[i].selected=maxselect1_1[i];}
	center_molecule(&m[0]);
	center_molecule(&m[1]);
	//printf("calling superimpose 3\n");
	rms=superimpose_molecules(&m[0],&m[1],s,0.0001);
	
      }
      //MAXSUB:   \t1.5\t2.0\t2.5\t3.0\t3.5\n
      if (b_flag){
	for (i=0;i<m[0].atoms;i++){
	  //printf("SCORE-TEST\t%d\t%f\t%e\n",i,bestsizescore[i],D_pvalue(i,bestsizescore[i]));
	  printf("SCORE-TEST\t%d\t%f\t%e\n",i,bestsizescore[i],D_pvalueF(i,bestsizescore[i]));
	}
      }
      if (out_flag){
	for (k=0;k<m[0].atoms;k++){m[0].atm[k].selected=maxselect1_0[k];
	  //	    printf("TEST %d %d\n",k,m[0].atm[k].selected);
	}
	for (k=0;k<m[1].atoms;k++){m[1].atm[k].selected=maxselect1_1[k];}
	center_molecule(&m[0]);
	center_molecule(&m[1]);
	rms=superimpose_molecules(&m[0],&m[1],s,0.0001);	   
	//printf("TEST-out: %s\n",outfile);
      }
  
  
  
      for (i=0;i<m[0].atoms;i++){m[0].atm[i].selected=maxselect1_0[i];}
      for (i=0;i<m[1].atoms;i++){m[1].atm[i].selected=maxselect1_1[i];}

    

    center_molecule(&m[0]);
    center_molecule(&m[1]);
    rms=superimpose_molecules(&m[0],&m[1],s,0.0001);

  //for (i=0;i<m[0].atoms;i++)
  //  {
  //   printf("RMS: %d %d %lf\n",i+1,m[0].atm[i].resnum,m[0].atm[i].rms);
  //  }
  //Sstr=(double*)malloc(sizeof(double)*(m[2].atoms+1));
  //Sstr_return=(double*)malloc(sizeof(double)*(m[2].atm[m[2].atoms-1].resnum+1));
  for(i=0;i<m[2].atm[m[2].atoms-1].resnum+1;i++)
    {
      //printf("%d %d %d %d %d\n",i,m[2].atoms,m[2].atm[0].resnum,m[2].atm[m[2].atoms-1].resnum,m[2].atoms+m[2].atm[0].resnum);
      D[0].S[i]=0;
      D[0].TM[i]=0;
    }
  D[0].S[i-1]=-1;
  D[0].TM[i-1]=-1;
  //exit(1);
  
  

  //  d0=25;
  //Sstr[0]=m[2].atoms+1;
  D[0].residues=0;
  
  for(i=0;i<length;i++)
    {
      //printf("%d %d %d %d %d\n",i,m[2].atoms,m[2].atm[0].resnum,m[2].atm[m[2].atoms-1].resnum,m[2].atoms+m[2].atm[0].resnum);
      D[0].S[i]=0;
      D[0].TM[i]=0;
    }
 

  for (i=0,j=0;i<m[2].atoms && j<m[0].atoms;i++)
    {
      
      if(m[2].atm[i].resnum == m[0].atm[j].resnum) // && m[0].atm[j].selected)
	{
	  // printf("YES! %d %d %d %lf\n",i,m[2].atm[i].resnum,m[0].atm[j].resnum,m[0].atm[j].rms);
	  if(isnan(m[0].atm[j].rms))
	    {
	      D[0].S[i]=0.0;
	      D[0].TM[i]=0.0;
	      // if(m[0].atm[j].selected)
	      //  temp2+=1/(1+m[0].atm[j].rms/d0);
	      m[2].atm[i].rms=999999999999; //998001;998001;
	      m[0].atm[j].rms=999999999999; //998001;998001;
	    }
	  else
	    {
	      //printf("%lf\n",d0);
	      D[0].S[i]=1/(1+m[0].atm[j].rms/(d0*d0));
	      D[0].TM[i]=1/(1+m[0].atm[j].rms/(d0_TM*d0_TM));
	      // if(m[0].atm[j].selected)
	      //  temp2+=1/(1+m[0].atm[j].rms/d0);
	      m[2].atm[i].rms=m[0].atm[j].rms;
	    }
	  j++;
	  
	}
      else
	{
	  D[0].S[i]=0;
	  D[0].TM[i]=0;
	  m[2].atm[i].rms=999999999999; //998001;
	  //	  D[0].S[i]=1/(1+998001/(d0*d0));
	  //	  D[0].TM[i]=1/(1+998001/(d0_TM*d0_TM));
	}
      //Sstr[i]=1/(1+m[2].atm[i].rms/5);
      //if(m[2].atm[i].selected)
      //	temp2+=1/(1+m[2].atm[i].rms/5);
      //Sstr_return[m[2].atm[i].resnum-1]=Sstr[i];
      Ssum=Ssum+D[0].S[i];
      TMsum=TMsum+D[0].TM[i];
      
      D[0].residue[i]=aa321(m[2].atm[i].residue);
      D[0].resnum[i]=m[2].atm[i].resnum;
      D[0].residues++;
      // printf("S: %d %d %c %lf %lf %d\n",i,m[2].atm[i].resnum,aa321(m[2].atm[i].residue),D[0].S[i],D[0].TM[i],length); //,m[0].atm[j-1].rms);
       //Sstr_return[m[2].atm[i].resnum]
    }
  D[0].Ssum=Ssum;
  D[0].TMsum=TMsum;
  // D[0].residues=m[2].residues;
  // printf("%d \n",D[0].residues);
  //exit(1);
  D[0].Dscore=log10(1/maxpvalue1);
  
  //printf("Ssum: %8.5lf over %d residues gives a mean of %8.5lf\n",Ssum,m[2].residues,Ssum/m[2].residues);
  //printf("TMsum: %8.5lf over %d residues gives a mean of %8.5lf\n",TMsum,m[2].residues,TMsum/m[2].residues);
  //printf("MAT: %lf %lf %lf\nMAT: %lf %lf %lf\nMAT: %lf %lf %lf\n",s[0][0],s[0][1],s[0][2],s[1][0],s[1][1],s[1][2],s[2][0],s[2][1],s[2][2]);
  //printf("%d\n",sizeof(Sstr));
  D[0].length=length;
  D[0].maxatoms=maxatoms1;
  D[0].maxrms=maxrms1;
  D[0].maxscore=maxscore1;
  D[0].maxpvalue=maxpvalue1;

  //printf("SCORE:\t%d\t%d\t%.1f\t%.1f\t%e\t%e\t%7.5lf\n",length,maxatoms1,maxrms1,maxscore1,maxpvalue1,ss,D[0].Dscore);
 

  // free(Sstr);
  //return Sstr_return;
  //exit(0);
  free(m);
}

