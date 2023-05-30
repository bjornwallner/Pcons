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

double LGscore(char *file1,char *file2,double minsim,int L, double factor)
{

  
  molecule	m[2];		/* Molecules to be compared */
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
  double          logLGscore=0;
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
  //  double          cutoffrms=3.;
  // double          minsim=25.0;
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
 
  double     x,y,z;
  FILE            *fp;  

  /* Initialize */
  files=2;
  verbose_flag=FALSE;
  out_flag=FALSE;
  select2_flag=FALSE;
  select1_flag=FALSE;
  short_flag=FALSE;
  arglist=FALSE;
  fraction_flag=FALSE;
  align_flag=TRUE;

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
  //printf("%-40s %-40s\t",file1,file2);
  //TETRA if(read_molecules_ca(&m[0])==0 && read_molecules_ca(&m[1])==0)
  if(read_molecules(&m[0],'c')==0 && read_molecules(&m[1],'c')==0)
    {
      length=max(m[0].atoms,m[1].atoms);
      check_molecules(&m[0],&m[1]); // this deletes some atoms
      atoms=m[0].atoms;
      for (i=0;i<m[0].atoms;i++){bestsizescore[i]=0;}
      for (i=0;i<m[0].atoms-L;i=i+1){
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
	//printf("calling superimpose 1\n");
	rms=superimpose_molecules(&m[0],&m[1],s,0.0001);	   
	count=0;
	while (j<m[0].atoms && rms<(double)(factor*(j+225)/225)){
	  count++;
	  //printf("%d %d\n",j,count);
	  minrmsd=9999.0;
	  for (k=0;k<m[0].atoms;k++) {
	    if(!m[0].atm[k].selected)
	      {
		if(m[0].atm[k].rms < minrmsd)
		  {
		    minpos=k;
		    minrmsd=m[0].atm[k].rms;
		  }
		else if(m[0].atm[k].rms < minsim)
		  {
		    //printf("Selecting %d %f\n",k,m[0].atm[k].rms);
		    m[0].atm[k].selected=TRUE;
		    m[1].atm[k].selected=TRUE;
		    j++;
		  }
	      }
	  }
	  
	  if(!m[0].atm[minpos].selected){
	    j++;
	  //printf("Selecting %d %f %d\n",minpos,minrmsd,j);
	    m[0].atm[minpos].selected=TRUE;
	    m[1].atm[minpos].selected=TRUE;
	  }
	  //j++;
	  ////printf("Selecting %d %f %d\n",minpos,minrmsd,j);
	  //m[0].atm[minpos].selected=TRUE;
	  //m[1].atm[minpos].selected=TRUE;
		  
	  //	     printf("testing2a %f\n",m[1].atm[i].x);
	  center_molecule(&m[0]);
	  center_molecule(&m[1]);
	  //printf("calling superimpose 2\t %d %d\n",i,j);
	  rms=superimpose_molecules(&m[0],&m[1],s,0.0001);
	  score=Levitt_Gerstein(&m[0],&m[1],5);
	  //pvalue=LG_pvalue(j,score);
	  pvalue=LG_pvalueF(j,score);
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
	  if ((pvalue <= maxpvalue) && (j>minatoms))
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
	    }
	}
	//	     printf("testing5 %f\n",m[1].atm[i].x);
      }
      //printf("%d %d\n",count,j);
      //	 printf("BEST1:\t%.2f\t%d\t%.1f\t%.1f\t%e\n",
      // (double)maxatoms/(double)atoms,maxatoms,maxrms,maxscore,maxpvalue);
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
	  ss=-1000*(log10(maxpvalue1)+0.6108)/m[0].atoms; //m[0].atoms;
	  //	  printf("%lf %d\n",ss,m[0].atoms);
	  //Why ma[0]????? dependent on order
	  //printf("SCORE:\t%d\t%d\t%.1f\t%.1f\t%e\t%e\t%7.5lf\n",
	  // length,maxatoms1,maxrms1,maxscore1,maxpvalue1,ss,(log10(1/maxpvalue1)));
	  logLGscore=(log10(1/maxpvalue1));
	}
      if (align_flag){
	//for (i=0;i<m[0].atoms;i++){m[0].atm[i].selected=maxselect1_0[i];}
	//for (i=0;i<m[1].atoms;i++){m[1].atm[i].selected=maxselect1_1[i];}
	//center_molecule(&m[0]);
	//center_molecule(&m[1]);
	//printf("calling superimpose 3\n");
	//rms=superimpose_molecules(&m[0],&m[1],s);
	//for(i=0;i<3;i++)
	//  {
	//    for(j=0;j<3;j++)
	//      {
	//		printf("%8.5lf ",s[i][j]);
	//     }
	//    printf("\n");
	//  }
	read_molecules(&m[0],'c');
	read_molecules(&m[1],'c');
	//fprintf(stderr,"LAST\n");
	center_molecule(&m[0]);
	center_molecule(&m[1]);
	for (i=0;i<m[0].atoms;i++){m[0].atm[i].selected=maxselect1_0[i];}
	for (i=0;i<m[1].atoms;i++){m[1].atm[i].selected=maxselect1_1[i];}
	rms=superimpose_molecules(&m[0],&m[1],s,0.0001);
	
	//for(i=0;i<3;i++)
	//  {
	//    for(j=0;j<3;j++)
	//      {
	//	printf("%8.5lf ",s[i][j]);
	//      }
	//    printf("\n");
	//  }
	//exit(1);
	for (i=0;i<m[1].atoms;i++)
	  {

	    //dist=m2->atm[i].x*m2->atm[i].x+m2->atm[i].y*m2->atm[i].y+m2->atm[i].z*m2->atm[i].z;
	    //rintf(stderr,"TEST1-m2\t %d\t%6.3f*%6.3f %6.3f*%6.3f %6.3f*%6.3f\n",i,s[0][0],m[1].atm[i].x,s[0][1],m[1].atm[i].y,s[0][2],m[1].atm[i].z);
	    //printf(stderr,"TEST2-m2\t %d\t%6.3f*%6.3f %6.3f*%6.3f %6.3f*%6.3f\n",i,s[1][0],m[1].atm[i].x,s[1][1],m[1].atm[i].y,s[1][2],m[1].atm[i].z);
	    //fprintf(stderr,"TEST3-m2\t %d\t%6.3f*%6.3f %6.3f*%6.3f %6.3f*%6.3f\n",i,s[2][0],m[1].atm[i].x,s[2][1],m[1].atm[i].y,s[2][2],m[1].atm[i].z);
	    //	    fprintf(stderr,"TEST1-m2\t %d\t%6.3f %6.3f %6.3f\n",i,m[1].atm[i].x,m[1].atm[i].y,m[1].atm[i].z);
	    x=s[0][0]*m[1].atm[i].x+s[0][1]*m[1].atm[i].y+s[0][2]*m[1].atm[i].z;
	    y=s[1][0]*m[1].atm[i].x+s[1][1]*m[1].atm[i].y+s[1][2]*m[1].atm[i].z;
	    z=s[2][0]*m[1].atm[i].x+s[2][1]*m[1].atm[i].y+s[2][2]*m[1].atm[i].z;
	    m[1].atm[i].x=x;
	    m[1].atm[i].y=y;
	    m[1].atm[i].z=z;
	    //	    fprintf(stderr,"TEST2-m2\t %d\t%6.3f %6.3f %6.3f\n",i,m[1].atm[i].x,m[1].atm[i].y,m[1].atm[i].z);
	    
	  }
	//fp=fopen("foo10.pdb","w");
	//if (fp!=NULL)
	//  {
	//    for (i=0;i<m[1].atoms;i++)
	//      {
	//	fprintf(fp,"ATOM   %-4d  %-4s%-3s  %-4d     %7.3lf %7.3lf %7.3lf  1.0   %-1d.00     \n",
	//		m[1].atm[i].number,m[1].atm[i].name,
	//		m[1].atm[i].residue,m[1].atm[i].resnum,
	//		m[1].atm[i].x,m[1].atm[i].y,m[1].atm[i].z,
	//		m[0].atm[i].selected);
	//      }
	//    fprintf(fp,"END \n");
	//    /*	    printf("test %d  %d \n",m[1].atoms,m[0].atoms);*/
	//    for (i=0;i<m[0].atoms;i++)
	//      {
	//	fprintf(fp,"ATOM   %-4d  %-4s%-3s  %-4d     %7.3lf %7.3lf %7.3lf  1.0   %-1d.50     \n",
	//		m[0].atm[i].number,m[0].atm[i].name,
	//		m[0].atm[i].residue,m[0].atm[i].resnum,
	//		m[0].atm[i].x,m[0].atm[i].y,m[0].atm[i].z,
	//		m[0].atm[i].selected);
	//      }
	//    fprintf(fp,"END \n");
	//  }
	//fclose(fp);
	
      }
	//MAXSUB:   \t1.5\t2.0\t2.5\t3.0\t3.5\n
      if (b_flag){
	for (i=0;i<m[0].atoms;i++){
	  //printf("SCORE-TEST\t%d\t%f\t%e\n",i,bestsizescore[i],LG_pvalue(i,bestsizescore[i]));
	  printf("SCORE-TEST-LG\t%d\t%f\t%e\n",i,bestsizescore[i],LG_pvalueF(i,bestsizescore[i]));
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
    }
  else
    {
      printf("Unable to do what you told me... My fault... I'm stupied\n");
    }
  //}
  //}
      
  //}
  return logLGscore;
}


double superimpose(char *file1,char *file2,char* file3,double minsim,int L, double factor)
{

  
  molecule	m[4];		/* Molecules to be compared */
  molecule	all_atom[2];		/* Molecules to be compared */
  vector        temp_center,best_center;
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
  int	i,j,k,l,a,b;	/* Counter variable */
  //  int	L=4;	// Minimum length of fragment for FindMaxSub matching
  int    minatoms=19; /* Variable for min no of atoms for P-value*/
  double	s[3][3];	/* Final transformation matrix */
  double	rotation_matrix[3][3];	/* Final transformation matrix */
  double	best_rotation_matrix[3][3];	/* Final transformation matrix */
  double        t[3][3];	/* Temporary storage matrix */
  double          maxrmsd,minrmsd;
  double          fraction=1;
  double          numfrac=0;
  double          maxpvalue=1,maxpvalue1=1.;
  double          score,maxrms,maxrms1,maxscore,maxscore1;
  double          pvalue;
  double          logLGscore=0;
  int             maxatoms,maxpos,minpos,minpos1,minpos2,atoms;
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
  //  double          cutoffrms=3.;
  // double          minsim=25.0;
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
  int             pos1,pos2;
  double     x,y,z;
  FILE            *fp;  

  /* Initialize */
  files=2;
  verbose_flag=FALSE;
  out_flag=FALSE;
  select2_flag=FALSE;
  select1_flag=FALSE;
  short_flag=FALSE;
  arglist=FALSE;
  fraction_flag=FALSE;
  align_flag=TRUE;

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
  //printf("%-40s %-40s\t",file1,file2);
  //TETRA if(read_molecules_ca(&m[0])==0 && read_molecules_ca(&m[1])==0)
  if(read_molecules(&m[0],'c')==0 && read_molecules(&m[1],'c')==0)
    {
      //printf("Hello\n");
      //if(arglist)
      //	{
      //	  molecule m[2];
      ///	  copymolecule(&m[0],&m_list[a]);
      //  copymolecule(&m[1],&m_list[b]);
      //	  printf("%-30s %-30s\t",m[0].filename,m[1].filename);
      //	}
      length=max(m[0].atoms,m[1].atoms);
      // printf("testing %d %d\n",m[0].atoms,m[1].atoms);
      check_molecules(&m[0],&m[1]); // this delete some  
      
      copy_molecule(&m[0],&m[2]); //copy m[0] to m[2] 
      copy_molecule(&m[1],&m[3]); //copy m[1] to m[3] 
      for (i=0;i<m[0].atoms;i++){bestsizescore[i]=0;}
      
      for (i=0;i<m[0].atoms-L;i++){
	//	printf("%d %d\n",i,L);
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
	copy_coord(&m[2],&m[0]);
	copy_coord(&m[3],&m[1]);
	temp_center=center_molecule(&m[0]);
	center_molecule(&m[1]);
	//	printf("calling superimpose 1 %d\n",i);
	rms=superimpose_molecules(&m[0],&m[1],s,0.0001);
	//copy_matrix(s,rotation_matrix);
	//printf("MAT: %lf %lf %lf\nMAT: %lf %lf %lf\nMAT: %lf %lf %lf\n\n",s[0][0],s[0][1],s[0][2],s[1][0],s[1][1],s[1][2],s[2][0],s[2][1],s[2][2]);

	count=0;
	while (j<m[0].atoms && rms<(double)(factor*(j+225)/225)){
	  count++;
	  //  printf("1: %d %d %d\n",j,i,count);
	  minrmsd=9999.0;
	  for (k=0;k<m[0].atoms;k++) {
	    if(!m[0].atm[k].selected)
	      {
		if(m[0].atm[k].rms < minrmsd)
		  {
		    minpos=k;
		    minrmsd=m[0].atm[k].rms;
		  }
		if(m[0].atm[k].rms < minsim)
		  {
		    //printf("Selecting %d %f\n",k,m[0].atm[k].rms);
		    m[0].atm[k].selected=TRUE;
		    m[1].atm[k].selected=TRUE;
		    j++;
		  }
	      }
	  }
	  if(!m[0].atm[minpos].selected)
	    {
	      j++;
	      //printf("Selecting %d %f %d\n",minpos,minrmsd,j);
	      m[0].atm[minpos].selected=TRUE;
	      m[1].atm[minpos].selected=TRUE;
	    }
	    //j++;
	  ////printf("Selecting %d %f %d\n",minpos,minrmsd,j);
	  //m[0].atm[minpos].selected=TRUE;
	  //m[1].atm[minpos].selected=TRUE;
		  
	  //	     printf("testing2a %f\n",m[1].atm[i].x);
	  copy_coord(&m[2],&m[0]);
	  copy_coord(&m[3],&m[1]);
	  temp_center=center_molecule(&m[0]);
	  center_molecule(&m[1]);
	  //printf("calling superimpose 2\t %d %d\n",i,j);
	  // printf("CEN: %lf %lf %lf\n",temp_center.x,temp_center.y,temp_center.z);
	  rms=superimpose_molecules(&m[0],&m[1],s,0.0001);

	  //multiply_matrix(rotation_matrix,s,t);
	  //copy_matrix(t,rotation_matrix);
	  score=Levitt_Gerstein(&m[0],&m[1],5);
	  //pvalue=LG_pvalue(j,score);
	  pvalue=LG_pvalueF(j,score);
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
	  if ((pvalue <= maxpvalue) && (j>minatoms))
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

	      copy_matrix(s,best_rotation_matrix);
	      best_center.x=temp_center.x;
	      best_center.y=temp_center.y;
	      best_center.z=temp_center.z;
	    }
	}
	//	     printf("testing5 %f\n",m[1].atm[i].x);
      }

    

      //printf("%d %d\n",count,j);
      //	 printf("BEST1:\t%.2f\t%d\t%.1f\t%.1f\t%e\n",
      //		(double)maxatoms/(double)atoms,maxatoms,maxrms,maxscore,maxpvalue);
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
	  ss=-1000*(log10(maxpvalue1)+0.6108)/m[0].atoms; //m[0].atoms;
	  //	  printf("%lf %d\n",ss,m[0].atoms);
	  //Why ma[0]????? dependent on order
	  //printf("SCORE:\t%d\t%d\t%.1f\t%.1f\t%e\t%e\t%7.5lf\n",
	  // length,maxatoms1,maxrms1,maxscore1,maxpvalue1,ss,(log10(1/maxpvalue1)));
	  logLGscore=(log10(1/maxpvalue1));
	}
      if (align_flag){
	
	//	printf ("s1:\t%f\t%f\t%f\n",best_rotation_matrix[0][0],best_rotation_matrix[0][1],best_rotation_matrix[0][2]);
	//printf ("s2:\t%f\t%f\t%f\n",best_rotation_matrix[1][0],best_rotation_matrix[1][1],best_rotation_matrix[1][2]);
	//	printf ("s3:\t%f\t%f\t%f\n",best_rotation_matrix[2][0],best_rotation_matrix[2][1],best_rotation_matrix[2][2]);

	//	printf("CEN: %lf %lf %lf\n",best_center.x,best_center.y,best_center.z);
	printf("FRAC: %.2f\tATOMS: %d\tRMS: %.2f\tS: %.1f\tP: %e\tLG: %.3f\n%s rotated in the frame of %s and saved as %s\n",
	       (double)maxatoms/(double)m[0].atoms,maxatoms,maxrms,maxscore,maxpvalue,log10(1/maxpvalue),file2,file1,file3);
	strcpy(all_atom[0].filename,m[0].filename);
	strcpy(all_atom[1].filename,m[1].filename);

	read_molecules(&all_atom[0],'a');
	read_molecules(&all_atom[1],'a');


	for (i=0;i<m[0].atoms;i++){m[0].atm[i].selected=maxselect1_0[i];}
	for (i=0;i<m[1].atoms;i++){m[1].atm[i].selected=maxselect1_1[i];}

	
	for (i=0;i<all_atom[0].atoms;i++)
	  {
	     if(atom_selected(all_atom[0].atm[i].name,all_atom[0].atm[i].resname,&m[0]))
	       {
		 all_atom[0].atm[i].selected=TRUE;
		 //printf("1: %-6s%5d %4s%1s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n","ATOM",all_atom[0].atm[i].number,all_atom[0].atm[i].name," ",all_atom[0].atm[i].residue,all_atom[0].atm[i].chain,all_atom[0].atm[i].resnum,all_atom[0].atm[i].x,all_atom[0].atm[i].y,all_atom[0].atm[i].z,1.00,all_atom[0].atm[i].bfactor);
	       }
	     else
	       {
		 all_atom[0].atm[i].selected=FALSE;
	       }
	  }
	//exit(1);
	for (i=0;i<all_atom[1].atoms;i++)
	  {
	     if(atom_selected(all_atom[1].atm[i].name,all_atom[1].atm[i].resname,&m[1]))
	       {
		 all_atom[1].atm[i].selected=TRUE;
		 // printf("2: %-6s%5d %4s%1s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n","ATOM",all_atom[1].atm[i].number,all_atom[1].atm[i].name," ",all_atom[1].atm[i].residue,all_atom[1].atm[i].chain,all_atom[1].atm[i].resnum,all_atom[1].atm[i].x,all_atom[1].atm[i].y,all_atom[1].atm[i].z,1.00,all_atom[1].atm[i].bfactor);
	       }
	     else
	       {
		 all_atom[1].atm[i].selected=FALSE;
	       }
	  }

	//	exit(1);
	//	center_molecule(&all_atom[0]);
	center_molecule(&all_atom[1]);
	
	for (i=0;i<all_atom[1].atoms;i++)
	  {
	    
	    //dist=m2->atm[i].x*m2->atm[i].x+m2->atm[i].y*m2->atm[i].y+m2->atm[i].z*m2->atm[i].z;
	    //rintf(stderr,"TEST1-m2\t %d\t%6.3f*%6.3f %6.3f*%6.3f %6.3f*%6.3f\n",i,s[0][0],m[1].atm[i].x,s[0][1],m[1].atm[i].y,s[0][2],m[1].atm[i].z);
	    //printf(stderr,"TEST2-m2\t %d\t%6.3f*%6.3f %6.3f*%6.3f %6.3f*%6.3f\n",i,s[1][0],m[1].atm[i].x,s[1][1],m[1].atm[i].y,s[1][2],m[1].atm[i].z);
	    //fprintf(stderr,"TEST3-m2\t %d\t%6.3f*%6.3f %6.3f*%6.3f %6.3f*%6.3f\n",i,s[2][0],m[1].atm[i].x,s[2][1],m[1].atm[i].y,s[2][2],m[1].atm[i].z);
	    //	    fprintf(stderr,"TEST1-m2\t %d\t%6.3f %6.3f %6.3f\n",i,m[1].atm[i].x,m[1].atm[i].y,m[1].atm[i].z);
	    x=best_rotation_matrix[0][0]*all_atom[1].atm[i].x+best_rotation_matrix[0][1]*all_atom[1].atm[i].y+best_rotation_matrix[0][2]*all_atom[1].atm[i].z;
	    y=best_rotation_matrix[1][0]*all_atom[1].atm[i].x+best_rotation_matrix[1][1]*all_atom[1].atm[i].y+best_rotation_matrix[1][2]*all_atom[1].atm[i].z;
	    z=best_rotation_matrix[2][0]*all_atom[1].atm[i].x+best_rotation_matrix[2][1]*all_atom[1].atm[i].y+best_rotation_matrix[2][2]*all_atom[1].atm[i].z;
	    x+=best_center.x;
	    y+=best_center.y;
	    z+=best_center.z;
	    
	    all_atom[1].atm[i].x=x;
	    all_atom[1].atm[i].y=y;
	    all_atom[1].atm[i].z=z;
	    //	    fprintf(stderr,"TEST2-m2\t %d\t%6.3f %6.3f %6.3f\n",i,m[1].atm[i].x,m[1].atm[i].y,m[1].atm[i].z);
	    
	  }



	fp=fopen(file3,"w");
	if (fp!=NULL)
	  {
	    for (i=0;i<all_atom[1].atoms;i++)
	      {
		//	printf $filehandle ("%-6s%5d %-4s%1s%3s %1s%5s   %8.3f%8.3f%8.3f%6.2f%6.2f\n", $self->{field}[$i],$self->{atomno}[$i],$self->{atomtype}[$i],$self->{atom_alt_loc}[$i],$self->{res}[$i],$chain,$lookup_hash{$self->{resno}[$i].$self->{insertion_code}[$i]},$self->{x}[$i],$self->{y}[$i],$self->{z}[$i],1,$self->{bfactor}[$i]);
		
		//printf("%d %s\n",i,m[1].atm[i].chain);
		fprintf(fp,"%-6s%5d %4s%1s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n","ATOM",all_atom[1].atm[i].number,all_atom[1].atm[i].name," ",all_atom[1].atm[i].residue,all_atom[1].atm[i].chain,all_atom[1].atm[i].resnum,all_atom[1].atm[i].x,all_atom[1].atm[i].y,all_atom[1].atm[i].z,1.00,all_atom[1].atm[i].bfactor);		
		//fprintf(fp,"%-6s%5d %-4s%1s%3s %1s%5s   %8.3f%8.3f%8.3f%6.2f%6.2f\n","ATOM",
		//	m[1].atm[i].number,m[1].atm[i].name,"",
		//	m[1].atm[i].residue,m[1].atm[i].chain,m[1].atm[i].resnum,
		//	m[1].atm[i].x,m[1].atm[i].y,m[1].atm[i].z,1,m[1].atm[i].bfactor);
		//m[0].atm[i].selected);
	      }
	    fprintf(fp,"TER\nEND\n");
	    ///*	    printf("test %d  %d \n",m[1].atoms,m[0].atoms);*/
	      //for (i=0;i<all_atom[0].atoms;i++)
	      //{
	      //	fprintf(fp,"%-6s%5d %4s%1s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n","ATOM",all_atom[0].atm[i].number,all_atom[0].atm[i].name," ",all_atom[0].atm[i].residue,all_atom[0].atm[i].chain,all_atom[0].atm[i].resnum,all_atom[0].atm[i].x,all_atom[0].atm[i].y,all_atom[0].atm[i].z,1.00,all_atom[0].atm[i].bfactor);
	      //
	      //}
	      //fprintf(fp,"TER \n");
	      //fprintf(fp,"END \n");

	  }
	fclose(fp);

	//check_molecules(&all_atom[0],&all_atom[1]);
	//
	//for (i=0;i<all_atom[0].atoms;i++)
	//  {
	//     if(atom_selected(all_atom[0].atm[i].name,all_atom[0].atm[i].resname,&m[0]))
	//       {
	//	 all_atom[0].atm[i].selected=TRUE;
	//	 //printf("1: %-6s%5d %4s%1s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n","ATOM",all_atom[0].atm[i].number,all_atom[0].atm[i].name," ",all_atom[0].atm[i].residue,all_atom[0].atm[i].chain,all_atom[0].atm[i].resnum,all_atom[0].atm[i].x,all_atom[0].atm[i].y,all_atom[0].atm[i].z,1.00,all_atom[0].atm[i].bfactor);
	//       }
	//     else
	//       {
	//	 all_atom[0].atm[i].selected=FALSE;
	//       }
	//  }
	//exit(1);
	//for (i=0;i<all_atom[1].atoms;i++)
	//  {
	//     if(atom_selected(all_atom[1].atm[i].name,all_atom[1].atm[i].resname,&m[1]))
	//       {
	//	 all_atom[1].atm[i].selected=TRUE;
	//	 // printf("2: %-6s%5d %4s%1s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n","ATOM",all_atom[1].atm[i].number,all_atom[1].atm[i].name," ",all_atom[1].atm[i].residue,all_atom[1].atm[i].chain,all_atom[1].atm[i].resnum,all_atom[1].atm[i].x,all_atom[1].atm[i].y,all_atom[1].atm[i].z,1.00,all_atom[1].atm[i].bfactor);
	//       }
	//     else
	//       {
	//	 all_atom[1].atm[i].selected=FALSE;
	//       }
	//  }

	//rms=0;
	//atoms=0;
	//or (i=0;i<all_atom[0].atoms;i++)
	// {
	//   if(all_atom[0].atm[i].selected)
	//     {
	//	//printf("1: %-6s%5d %4s%1s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n","ATOM",all_atom[0].atm[i].number,all_atom[0].atm[i].name," ",all_atom[0].atm[i].residue,all_atom[0].atm[i].chain,all_atom[0].atm[i].resnum,all_atom[0].atm[i].x,all_atom[0].atm[i].y,all_atom[0].atm[i].z,1.00,all_atom[0].atm[i].bfactor);
	//  	   
	//	//printf("2: %-6s%5d %4s%1s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n","ATOM",all_atom[1].atm[i].number,all_atom[1].atm[i].name," ",all_atom[1].atm[i].residue,all_atom[1].atm[i].chain,all_atom[1].atm[i].resnum,all_atom[1].atm[i].x,all_atom[1].atm[i].y,all_atom[1].atm[i].z,1.00,all_atom[1].atm[i].bfactor);
	//  
	//	//    rms+=(all_atom[0].atm[i].x-all_atom[1].atm[i].x)*(all_atom[0].atm[i].x-all_atom[1].atm[i].x)+(all_atom[0].atm[i].y-all_atom[1].atm[i].y)*(all_atom[0].atm[i].y-all_atom[1].atm[i].y)+(all_atom[0].atm[i].z-all_atom[1].atm[i].z)*(all_atom[0].atm[i].z-all_atom[1].atm[i].z);
	//	//atoms++;
	//     }
	//     //printf("RMSD: %lf\n",rms);
	// }
	//	printf("RMSD: %lf %d\n",sqrt(rms/atoms),atoms);


	
      }
	//MAXSUB:   \t1.5\t2.0\t2.5\t3.0\t3.5\n
      if (b_flag){
	for (i=0;i<m[0].atoms;i++){
	  //printf("SCORE-TEST\t%d\t%f\t%e\n",i,bestsizescore[i],LG_pvalue(i,bestsizescore[i]));
	  printf("SCORE-TEST\t%d\t%f\t%e\n",i,bestsizescore[i],LG_pvalueF(i,bestsizescore[i]));
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
    }
  else
    {
      printf("Unable to do what you told me... My fault... I'm stupied\n");
    }
  //}
  //}
      
  //}
  return logLGscore;
}



void LGscore_res(char* file1,char* file2,lgscore *LG, double d0, double minsim,int L,double factor,int step)
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
  double          logLGscore=0;
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
  strcpy(m[2].filename,file1);
  //printf("%-40s %-40s\t",file1,file2);
  //TETRA if(read_molecules_ca(&m[0])==0 && read_molecules_ca(&m[1])==0 && read_molecules_ca(&m[2])==0)
  if(read_molecules(&m[0],'c')==0 && read_molecules(&m[1],'c')==0 && read_molecules(&m[2],'c')==0)
    {
      //      printf("Hello\n");
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
	  rms=superimpose_molecules(&m[0],&m[1],s,0.1);  //not so strict error cutoff on the superimpose
#ifdef Sscore
	  score=Levitt_Gerstein(&m[0],&m[1],d0*d0);
#else
	  score=Levitt_Gerstein(&m[0],&m[1],5);
#endif
	  //pvalue=LG_pvalue(j,score);
	  pvalue=LG_pvalueF(j,score);
	  
	  //printf("TEST: %e\n", LG_pvalueF(227,1532.3));
	  if (b_flag)
	    {
	      if (score > bestsizescore[j]){
		bestsizescore[j]=score;
	      }
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

	  //	  if ((score >= maxscore) && (j>minatoms)){
	  //LGscore
	  //  if ((pvalue <= maxpvalue) && (j>minatoms))
#ifdef Sscore
	  if ((score >= maxscore) && (j>minatoms))
	    //printf("Sscore rule!\n");
	    //	  if ((pvalue <= maxpvalue) && (j>minatoms)) {  
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
	    }
	}
	//	     printf("testing5 %f\n",m[1].atm[i].x);
      }
      //printf("%d %d\n",count,j);
      //	 printf("BEST1:\t%.2f\t%d\t%.1f\t%.1f\t%e\n",
      // (double)maxatoms/(double)atoms,maxatoms,maxrms,maxscore,maxpvalue);
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
	  logLGscore=(log10(1/maxpvalue1));
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
	  //printf("SCORE-TEST\t%d\t%f\t%e\n",i,bestsizescore[i],LG_pvalue(i,bestsizescore[i]));
	  printf("SCORE-TEST\t%d\t%f\t%e\n",i,bestsizescore[i],LG_pvalueF(i,bestsizescore[i]));
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
      LG[0].S[i]=0;
      LG[0].TM[i]=0;
    }
  LG[0].S[i-1]=-1;
  LG[0].TM[i-1]=-1;
  //exit(1);
  
  //  d0=25;
  //Sstr[0]=m[2].atoms+1;
  LG[0].residues=0;
  for (i=0,j=0;i<m[2].atoms;i++)
    {
      
      if(m[2].atm[i].resnum == m[0].atm[j].resnum) // && m[0].atm[j].selected)
	{
	  if(isnan(m[0].atm[j].rms))
	    {
	      LG[0].S[i]=0.0;
	      LG[0].TM[i]=0.0;
	      // if(m[0].atm[j].selected)
	      //  temp2+=1/(1+m[0].atm[j].rms/d0);
	      m[2].atm[i].rms=998001;
	      m[0].atm[j].rms=998001;
	    }
	  else
	    {
	      //printf("%lf\n",d0);
	      LG[0].S[i]=1/(1+m[0].atm[j].rms/(d0*d0));
	      LG[0].TM[i]=1/(1+m[0].atm[j].rms/(d0_TM*d0_TM));
	      // if(m[0].atm[j].selected)
	      //  temp2+=1/(1+m[0].atm[j].rms/d0);
	      m[2].atm[i].rms=m[0].atm[j].rms;
	    }
	  j++;
	  
	}
      else
	{
	  LG[0].S[i]=0;
	  LG[0].TM[i]=0;
	  m[2].atm[i].rms=998001;
	  LG[0].S[i]=1/(1+998001/(d0*d0));
	  LG[0].TM[i]=1/(1+998001/(d0_TM*d0_TM));
	}
      //Sstr[i]=1/(1+m[2].atm[i].rms/5);
      //if(m[2].atm[i].selected)
      //	temp2+=1/(1+m[2].atm[i].rms/5);
      //Sstr_return[m[2].atm[i].resnum-1]=Sstr[i];
      Ssum=Ssum+LG[0].S[i];
      TMsum=TMsum+LG[0].TM[i];
      
      LG[0].residue[i]=aa321(m[2].atm[i].residue);
      LG[0].resnum[i]=m[2].atm[i].resnum;
      LG[0].residues++;
      // printf("%d %c %lf %lf\n",m[2].atm[i].resnum,aa321(m[2].atm[i].residue),LG[0].S[i],LG[0].TM[i]); //,m[0].atm[j-1].rms);
       //Sstr_return[m[2].atm[i].resnum]
    }
  LG[0].Ssum=Ssum;
  LG[0].TMsum=TMsum;
  // LG[0].residues=m[2].residues;
  // printf("%d \n",LG[0].residues);
  //exit(1);
  LG[0].LGscore=log10(1/maxpvalue1);

  //printf("Ssum: %8.5lf over %d residues gives a mean of %8.5lf\n",Ssum,m[2].residues,Ssum/m[2].residues);
  //printf("TMsum: %8.5lf over %d residues gives a mean of %8.5lf\n",TMsum,m[2].residues,TMsum/m[2].residues);
  //printf("MAT: %lf %lf %lf\nMAT: %lf %lf %lf\nMAT: %lf %lf %lf\n",s[0][0],s[0][1],s[0][2],s[1][0],s[1][1],s[1][2],s[2][0],s[2][1],s[2][2]);
  //printf("%d\n",sizeof(Sstr));
  LG[0].length=length;
  LG[0].maxatoms=maxatoms1;
  LG[0].maxrms=maxrms1;
  LG[0].maxscore=maxscore1;
  LG[0].maxpvalue=maxpvalue1;

  //printf("SCORE:\t%d\t%d\t%.1f\t%.1f\t%e\t%e\t%7.5lf\n",length,maxatoms1,maxrms1,maxscore1,maxpvalue1,ss,LG[0].LGscore);
 

  // free(Sstr);
  //return Sstr_return;
  //exit(0);
  
}
  
  
void rms(dyn_molecule *m1,dyn_molecule *m2,lgscore *LG, double d0)
  {
  
    molecule	m[3];		/* Molecules to be compared */
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
    double          logLGscore=0;
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
    
    copymolecule(&m[0],m1); //copies only the CA...
    copymolecule(&m[1],m2);
    copymolecule(&m[2],m1);
    

  
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
      
    reset_rms(&m[0]);
    reset_rms(&m[1]);
    for (j=0;j<m[0].atoms;j++){ // select all atoms
      m[0].atm[j].selected=TRUE;
      m[1].atm[j].selected=TRUE;
    }
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
      LG[0].S[i]=0;
      LG[0].TM[i]=0;
    }
  LG[0].S[i-1]=-1;
  LG[0].TM[i-1]=-1;
  //exit(1);
  
  

  //  d0=25;
  //Sstr[0]=m[2].atoms+1;
  LG[0].residues=0;
  
  for(i=0;i<length;i++)
    {
      //printf("%d %d %d %d %d\n",i,m[2].atoms,m[2].atm[0].resnum,m[2].atm[m[2].atoms-1].resnum,m[2].atoms+m[2].atm[0].resnum);
      LG[0].S[i]=0;
      LG[0].TM[i]=0;
    }
 

  for (i=0,j=0;i<m[2].atoms && j<m[0].atoms;i++)
    {
      
      if(m[2].atm[i].resnum == m[0].atm[j].resnum) // && m[0].atm[j].selected)
	{
	  // printf("YES! %d %d %d %lf\n",i,m[2].atm[i].resnum,m[0].atm[j].resnum,m[0].atm[j].rms);
	  if(isnan(m[0].atm[j].rms))
	    {
	      LG[0].S[i]=0.0;
	      LG[0].TM[i]=0.0;
	      // if(m[0].atm[j].selected)
	      //  temp2+=1/(1+m[0].atm[j].rms/d0);
	      m[2].atm[i].rms=999999999999; //998001;998001;
	      m[0].atm[j].rms=999999999999; //998001;998001;
	    }
	  else
	    {
	      //printf("%lf\n",d0);
	      LG[0].S[i]=1/(1+m[0].atm[j].rms/(d0*d0));
	      LG[0].TM[i]=1/(1+m[0].atm[j].rms/(d0_TM*d0_TM));
	      // if(m[0].atm[j].selected)
	      //  temp2+=1/(1+m[0].atm[j].rms/d0);
	      m[2].atm[i].rms=m[0].atm[j].rms;
	    }
	  j++;
	  
	}
      else
	{
	  LG[0].S[i]=0;
	  LG[0].TM[i]=0;
	  m[2].atm[i].rms=999999999999; //998001;
	  //	  LG[0].S[i]=1/(1+998001/(d0*d0));
	  //	  LG[0].TM[i]=1/(1+998001/(d0_TM*d0_TM));
	}
      //Sstr[i]=1/(1+m[2].atm[i].rms/5);
      //if(m[2].atm[i].selected)
      //	temp2+=1/(1+m[2].atm[i].rms/5);
      //Sstr_return[m[2].atm[i].resnum-1]=Sstr[i];
      Ssum=Ssum+LG[0].S[i];
      TMsum=TMsum+LG[0].TM[i];
      
      LG[0].residue[i]=aa321(m[2].atm[i].residue);
      LG[0].resnum[i]=m[2].atm[i].resnum;
      LG[0].residues++;
      // printf("S: %d %d %c %lf %lf %d\n",i,m[2].atm[i].resnum,aa321(m[2].atm[i].residue),LG[0].S[i],LG[0].TM[i],length); //,m[0].atm[j-1].rms);
       //Sstr_return[m[2].atm[i].resnum]
    }
  LG[0].Ssum=Ssum;
  LG[0].TMsum=TMsum;
  // LG[0].residues=m[2].residues;
  // printf("%d \n",LG[0].residues);
  //exit(1);
  LG[0].LGscore=log10(1/maxpvalue1);

  //printf("Ssum: %8.5lf over %d residues gives a mean of %8.5lf\n",Ssum,m[2].residues,Ssum/m[2].residues);
  //printf("TMsum: %8.5lf over %d residues gives a mean of %8.5lf\n",TMsum,m[2].residues,TMsum/m[2].residues);
  //printf("MAT: %lf %lf %lf\nMAT: %lf %lf %lf\nMAT: %lf %lf %lf\n",s[0][0],s[0][1],s[0][2],s[1][0],s[1][1],s[1][2],s[2][0],s[2][1],s[2][2]);
  //printf("%d\n",sizeof(Sstr));
  LG[0].length=length;
  LG[0].maxatoms=maxatoms1;
  LG[0].maxrms=maxrms1;
  LG[0].maxscore=maxscore1;
  LG[0].maxpvalue=maxpvalue;

  //printf("SCORE:\t%d\t%d\t%.1f\t%.1f\t%e\t%e\t%7.5lf\n",length,maxatoms1,maxrms1,maxscore1,maxpvalue1,ss,LG[0].LGscore);
 

  // free(Sstr);
  //return Sstr_return;
  //exit(0);
  
}



void LGscore_res_pt(dyn_molecule *m1,dyn_molecule *m2,lgscore *LG, double d0, double minsim,int L,double factor,int step)
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
    double          logLGscore=0;
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
	  //pvalue=LG_pvalue(j,score);
	  pvalue=LG_pvalueF(j,score);
	  
	  //printf("TEST: %e\n", LG_pvalueF(227,1532.3));
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
	  logLGscore=(log10(1/maxpvalue1));
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
	  //printf("SCORE-TEST\t%d\t%f\t%e\n",i,bestsizescore[i],LG_pvalue(i,bestsizescore[i]));
	  printf("SCORE-TEST\t%d\t%f\t%e\n",i,bestsizescore[i],LG_pvalueF(i,bestsizescore[i]));
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
      LG[0].S[i]=0;
      LG[0].TM[i]=0;
    }
  LG[0].S[i-1]=-1;
  LG[0].TM[i-1]=-1;
  //exit(1);
  
  

  //  d0=25;
  //Sstr[0]=m[2].atoms+1;
  LG[0].residues=0;
  
  for(i=0;i<length;i++)
    {
      //printf("%d %d %d %d %d\n",i,m[2].atoms,m[2].atm[0].resnum,m[2].atm[m[2].atoms-1].resnum,m[2].atoms+m[2].atm[0].resnum);
      LG[0].S[i]=0;
      LG[0].TM[i]=0;
    }
 

  for (i=0,j=0;i<m[2].atoms && j<m[0].atoms;i++)
    {
      
      if(m[2].atm[i].resnum == m[0].atm[j].resnum) // && m[0].atm[j].selected)
	{
	  // printf("YES! %d %d %d %lf\n",i,m[2].atm[i].resnum,m[0].atm[j].resnum,m[0].atm[j].rms);
	  if(isnan(m[0].atm[j].rms))
	    {
	      LG[0].S[i]=0.0;
	      LG[0].TM[i]=0.0;
	      // if(m[0].atm[j].selected)
	      //  temp2+=1/(1+m[0].atm[j].rms/d0);
	      m[2].atm[i].rms=999999999999; //998001;998001;
	      m[0].atm[j].rms=999999999999; //998001;998001;
	    }
	  else
	    {
	      //printf("%lf\n",d0);
	      LG[0].S[i]=1/(1+m[0].atm[j].rms/(d0*d0));
	      LG[0].TM[i]=1/(1+m[0].atm[j].rms/(d0_TM*d0_TM));
	      // if(m[0].atm[j].selected)
	      //  temp2+=1/(1+m[0].atm[j].rms/d0);
	      m[2].atm[i].rms=m[0].atm[j].rms;
	    }
	  j++;
	  
	}
      else
	{
	  LG[0].S[i]=0;
	  LG[0].TM[i]=0;
	  m[2].atm[i].rms=999999999999; //998001;
	  //	  LG[0].S[i]=1/(1+998001/(d0*d0));
	  //	  LG[0].TM[i]=1/(1+998001/(d0_TM*d0_TM));
	}
      //Sstr[i]=1/(1+m[2].atm[i].rms/5);
      //if(m[2].atm[i].selected)
      //	temp2+=1/(1+m[2].atm[i].rms/5);
      //Sstr_return[m[2].atm[i].resnum-1]=Sstr[i];
      Ssum=Ssum+LG[0].S[i];
      TMsum=TMsum+LG[0].TM[i];
      
      LG[0].residue[i]=aa321(m[2].atm[i].residue);
      LG[0].resnum[i]=m[2].atm[i].resnum;
      LG[0].residues++;
      // printf("S: %d %d %c %lf %lf %d\n",i,m[2].atm[i].resnum,aa321(m[2].atm[i].residue),LG[0].S[i],LG[0].TM[i],length); //,m[0].atm[j-1].rms);
       //Sstr_return[m[2].atm[i].resnum]
    }
  LG[0].Ssum=Ssum;
  LG[0].TMsum=TMsum;
  // LG[0].residues=m[2].residues;
  // printf("%d \n",LG[0].residues);
  //exit(1);
  LG[0].LGscore=log10(1/maxpvalue1);
  
  //printf("Ssum: %8.5lf over %d residues gives a mean of %8.5lf\n",Ssum,m[2].residues,Ssum/m[2].residues);
  //printf("TMsum: %8.5lf over %d residues gives a mean of %8.5lf\n",TMsum,m[2].residues,TMsum/m[2].residues);
  //printf("MAT: %lf %lf %lf\nMAT: %lf %lf %lf\nMAT: %lf %lf %lf\n",s[0][0],s[0][1],s[0][2],s[1][0],s[1][1],s[1][2],s[2][0],s[2][1],s[2][2]);
  //printf("%d\n",sizeof(Sstr));
  LG[0].length=length;
  LG[0].maxatoms=maxatoms1;
  LG[0].maxrms=maxrms1;
  LG[0].maxscore=maxscore1;
  LG[0].maxpvalue=maxpvalue1;

  //printf("SCORE:\t%d\t%d\t%.1f\t%.1f\t%e\t%e\t%7.5lf\n",length,maxatoms1,maxrms1,maxscore1,maxpvalue1,ss,LG[0].LGscore);
 

  // free(Sstr);
  //return Sstr_return;
  //exit(0);
  free(m);
}



vector center_molecule(molecule *m)	/* Reads in molecules to be superimposed */
     //    molecule *m;
{
  int	i;	/* Counter variables */
  int   natoms;  // Number of selected atoms
  double	xcen,ycen,zcen;		/* Temporary coordinates values */
  vector vec;
  
  natoms=0;
  xcen=ycen=zcen=0;
  for (i=0;i<m->atoms;i++){
    if (m->atm[i].selected){
      xcen+=m->atm[i].x;
      ycen+=m->atm[i].y;
      zcen+=m->atm[i].z;
      natoms++;
    }
  }
  /* Now center molecule */
  xcen/=(double)natoms;
  ycen/=(double)natoms;
  zcen/=(double)natoms;
  vec.x=xcen;
  vec.y=ycen;
  vec.z=zcen;
  
  //  printf("TEST CEN natoms: %d %f %f %f \n",natoms,xcen,ycen,zcen);
  for (i=0;i<m->atoms;i++)
    {
      m->atm[i].x-=xcen;
      m->atm[i].y-=ycen;
      m->atm[i].z-=zcen;
      // printf("TEST natoms: %d %f %f %f \n",i,m->atm[i].x,m->atm[i].y,m->atm[i].z);
    }
  //i--;
  //printf("TEST natoms: %d %f %f %f \n",i,m->atm[i].x,m->atm[i].y,m->atm[i].z);
  return vec;
}
int multiply_matrix(a,b,c)  /* computes C=AB */
     double a[3][3],b[3][3],c[3][3];
{
  int i,j,k;
  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      {
	c[i][j]=0.0;
	for (k=0;k<3;k++)
	  c[i][j]+=a[i][k]*b[k][j];
      }
}
int copy_matrix(f,t) /* copy matrix f into matrix t */
     double f[3][3],t[3][3];
{
  int i,j;
  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      t[i][j]=f[i][j];
}
int transpose_matrix(m) /* Transpose a 3x3 matrix */
     double m[3][3];
{
  double dummy;
  dummy=m[0][1]; m[0][1]=m[1][0]; m[1][0]=dummy;
  dummy=m[0][2]; m[0][2]=m[2][0]; m[2][0]=dummy;
  dummy=m[1][2]; m[1][2]=m[2][1]; m[2][1]=dummy;
}
int delete_atom(molecule *m1,int num)	
{
  int          i,j,k;
  for (i=num;i<m1->atoms;i++)
    {
      j=i+1;
      strcpy(m1->atm[i].name,m1->atm[j].name);
      strcpy(m1->atm[i].residue,m1->atm[j].residue);
      strcpy(m1->atm[i].resname,m1->atm[j].resname);
      m1->atm[i].x=m1->atm[j].x;
      m1->atm[i].y=m1->atm[j].y;
      m1->atm[i].z=m1->atm[j].z;
      /*      m1->atm[i].rms=m1->atm[j].rms; */
      m1->atm[i].number=m1->atm[j].number;
      m1->atm[i].resnum=m1->atm[j].resnum;
      

    }
  m1->atoms--;
}

int check_molecules2(molecule *m1,molecule *m2)	
{
  /* this will delete all atoms that are not in both molecules */
  /* Now we should extract the part of the molecules that exists in both */
  int          i,j,minlength;
  int          found;

  minlength=m1->atoms;
  if (minlength>m1->atoms) minlength=m2->atoms;
  i=0;
  for(i=0;i<m1->atoms;i++)
    {
      for(j=0;j<m2->atoms;j++)
	{
	  

	}

    }


  while (i<=m1->atoms &&  i<=m2->atoms ){
    while ( m1->atm[i].resnum < m2->atm[i].resnum && i< m1->atoms ){
      //  printf("Deleting m1: %d %d %d\n",i,m1->atm[i].resnum,m2->atm[i].resnum);
      delete_atom(m1,i);	
    }
    while (m1->atm[i].resnum > m2->atm[i].resnum &&  i< m2->atoms ){
      //printf("Deleting m2: %d %d %d\n",i,m1->atm[i].resnum,m2->atm[i].resnum);
      delete_atom(m2,i);	
    }
    i++;
  }
  while (m1->atoms < m2->atoms){    delete_atom(m2,m2->atoms)  ;}
  while (m1->atoms > m2->atoms){    delete_atom(m1,m1->atoms)  ;}
}
int check_molecules(m1,m2)	
    molecule 	*m1,*m2;		
{
  /* this will delete all atoms that are not in both molecules */
  /* Now we should extract the part of the molecules that exists in both */
  int          i,j,minlength;
  int          found;
  int          current_resnum;
  char         current_atom[4];
  minlength=m1->atoms;
  if (minlength>m2->atoms) minlength=m2->atoms;
  i=0;
  while(i<m1->atoms)
    {
      if(!atom_exist(m1->atm[i].name,m1->atm[i].resname,m2))
	{
	  //delete_atom(m1,i);
	  //printf("NO m1: %d %s %s %d \"%s\"\n",
	  //	 i,
	  //	 m1->atm[i].name,
	  //	 m1->atm[i].residue,
	  //	 m1->atm[i].resnum,
	  //	 m1->atm[i].resname);
		 
	  delete_atom(m1,i);
	}
      else
	{
	  i++;
	  //printf("NO\n");
	  //delete_atom(m1,i);
	}
    }
  //exit(1);
  i=0;
  while(i<m2->atoms)
    {
      if(!atom_exist(m2->atm[i].name,m2->atm[i].resname,m1))
	{
	  //  printf("NO m2: %s %s %d %d %d\n",m2->atm[i].name,m2->atm[i].residue,i,m2->atm[i].resnum,m1->atm[i].resnum);
	  delete_atom(m2,i);
	  //printf("NO m2: %s %s %d %d %d\n",m2->atm[i].name,m2->atm[i].residue,i,m2->atm[i].resnum,m1->atm[i].resnum);
	  //printf("YES m2: %s %s %d %d %d\n",m2->atm[i].name,m2->atm[i].residue,i,m2->atm[i].resnum,m1->atm[i].resnum);
	}
      else
	{
	  i++;
	}
    }

  //  for(i=0;i<m1->atoms;i++)
  //  {
  //	printf("1: %d %d %s %s %d %lf %lf %lf --- ",i,m1->atm[i].number,m1->atm[i].name,m1->atm[i].residue,m1->atm[i].resnum,m1->atm[i].x,m1->atm[i].y,m1->atm[i].z);
  //	printf("2: %d %d %s %s %d %lf %lf %lf\n",i,m2->atm[i].number,m2->atm[i].name,m2->atm[i].residue,m2->atm[i].resnum,m2->atm[i].x,m2->atm[i].y,m2->atm[i].z);
  //   }
//  while (i<=m1->atoms &&  i<=m2->atoms ){
//    //printf("%d %d %d\n",i,m1->atm[i].resnum,m2->atm[i].resnum);
//    printf("%d %d < %d\n",i,m1->atm[i].resnum,m2->atm[i].resnum);
//  while ( m1->atm[i].resnum < m2->atm[i].resnum && i< m1->atoms )
//    {
//	printf("Deleting m1: %s %s %d %d %d\n",m1->atm[i].name,m1->atm[i].residue,i,m1->atm[i].resnum,m2->atm[i].resnum);
//     delete_atom(m1,i);	
//      }
//    
//    while (m1->atm[i].resnum > m2->atm[i].resnum &&  i< m2->atoms )
//      {
//	printf("Deleting m2: %s %s %d %d %d\n",m2->atm[i].name,m2->atm[i].residue,i,m1->atm[i].resnum,m2->atm[i].resnum);
//	delete_atom(m2,i);	
//      }
//    i++;
//  }
//  while (m1->atoms < m2->atoms){    delete_atom(m2,m2->atoms)  ;}
//  while (m1->atoms > m2->atoms){    delete_atom(m1,m1->atoms)  ;}

//  printf("sizes: %d %d\n",m1->atoms,m2->atoms);
//  for(i=0;i<m1->atoms;i++)
//    {
//	printf("1: %d %d %s %s %d %lf %lf %lf --- ",i,m1->atm[i].number,m1->atm[i].name,m1->atm[i].residue,m1->atm[i].resnum,m1->atm[i].x,m1->atm[i].y,m1->atm[i].z);
//	printf("2: %d %d %s %s %d %lf %lf %lf\n",i,m2->atm[i].number,m2->atm[i].name,m2->atm[i].residue,m2->atm[i].resnum,m2->atm[i].x,m2->atm[i].y,m2->atm[i].z);
//     }
  //for(i=0;i<m2->atoms;i++)
  //  {
      //printf("2: %d %d %s %s %d %lf %lf %lf\n",i,m2->atm[i].number,m2->atm[i].name,m2->atm[i].residue,m2->atm[i].resnum,m2->atm[i].x,m2->atm[i].y,m2->atm[i].z);
  //   }
  
  if (m1->atoms > 0){return(0);}
  else{
    //fprintf(stderr,"no identical atoms in files %s %s \n",m1->filename,m2->filename);
    return(1);
  }
}


int check_molecules_mark_deleted(m1,m2)	
    molecule 	*m1,*m2;		
{
  /* this will mark atoms that are are not in both molecules as deleted*/
  /* Now we should extract the part of the molecules that exists in both */
  int          i,j,minlength;
  int          found;
  int          current_resnum;
  char         current_atom[4];
  minlength=m1->atoms;
  if (minlength>m2->atoms) minlength=m2->atoms;
  i=0;
  while(i<m1->atoms)
    {
      if(!atom_exist(m1->atm[i].name,m1->atm[i].resname,m2) || strcmp(m1->atm[i].name,"CA ")!=0)
	{
	  //delete_atom(m1,i);
	  //printf("NO m1: %d %s %s %d \"%s\" %d\n",
	  //	 i,
	  //	 m1->atm[i].name,
	  //	 m1->atm[i].residue,
	  //	 m1->atm[i].resnum,
	  //	 m1->atm[i].name,strcmp(m1->atm[i].name,"CA "));
		 
	  m1->atm[i].deleted=TRUE;
	  //delete_atom(m1,i);
	}
      i++;
	//      else
	//{
	// i++;
	  //printf("NO\n");
	  //delete_atom(m1,i);
	//}
    }
  //exit(1);
  i=0;
  while(i<m2->atoms)
    {
      if(!atom_exist(m2->atm[i].name,m2->atm[i].resname,m1) || strcmp(m2->atm[i].name,"CA ")!=0)
	{
	  //  printf("NO m2: %s %s %d %d %d\n",m2->atm[i].name,m2->atm[i].residue,i,m2->atm[i].resnum,m1->atm[i].resnum);
	  //delete_atom(m2,i);
	  m2->atm[i].deleted=TRUE;
	  //printf("NO m2: %s %s %d %d %d\n",m2->atm[i].name,m2->atm[i].residue,i,m2->atm[i].resnum,m1->atm[i].resnum);
	  //printf("YES m2: %s %s %d %d %d\n",m2->atm[i].name,m2->atm[i].residue,i,m2->atm[i].resnum,m1->atm[i].resnum);
	}
      // else
      //	{
	  i++;
	  //	}
    }

  //  for(i=0;i<m1->atoms;i++)
  //  {
  //	printf("1: %d %d %s %s %d %lf %lf %lf --- ",i,m1->atm[i].number,m1->atm[i].name,m1->atm[i].residue,m1->atm[i].resnum,m1->atm[i].x,m1->atm[i].y,m1->atm[i].z);
  //	printf("2: %d %d %s %s %d %lf %lf %lf\n",i,m2->atm[i].number,m2->atm[i].name,m2->atm[i].residue,m2->atm[i].resnum,m2->atm[i].x,m2->atm[i].y,m2->atm[i].z);
  //   }
//  while (i<=m1->atoms &&  i<=m2->atoms ){
//    //printf("%d %d %d\n",i,m1->atm[i].resnum,m2->atm[i].resnum);
//    printf("%d %d < %d\n",i,m1->atm[i].resnum,m2->atm[i].resnum);
//  while ( m1->atm[i].resnum < m2->atm[i].resnum && i< m1->atoms )
//    {
//	printf("Deleting m1: %s %s %d %d %d\n",m1->atm[i].name,m1->atm[i].residue,i,m1->atm[i].resnum,m2->atm[i].resnum);
//     delete_atom(m1,i);	
//      }
//    
//    while (m1->atm[i].resnum > m2->atm[i].resnum &&  i< m2->atoms )
//      {
//	printf("Deleting m2: %s %s %d %d %d\n",m2->atm[i].name,m2->atm[i].residue,i,m1->atm[i].resnum,m2->atm[i].resnum);
//	delete_atom(m2,i);	
//      }
//    i++;
//  }
//  while (m1->atoms < m2->atoms){    delete_atom(m2,m2->atoms)  ;}
//  while (m1->atoms > m2->atoms){    delete_atom(m1,m1->atoms)  ;}

//  printf("sizes: %d %d\n",m1->atoms,m2->atoms);
//  for(i=0;i<m1->atoms;i++)
//    {
//	printf("1: %d %d %s %s %d %lf %lf %lf --- ",i,m1->atm[i].number,m1->atm[i].name,m1->atm[i].residue,m1->atm[i].resnum,m1->atm[i].x,m1->atm[i].y,m1->atm[i].z);
//	printf("2: %d %d %s %s %d %lf %lf %lf\n",i,m2->atm[i].number,m2->atm[i].name,m2->atm[i].residue,m2->atm[i].resnum,m2->atm[i].x,m2->atm[i].y,m2->atm[i].z);
//     }
  //for(i=0;i<m2->atoms;i++)
  //  {
      //printf("2: %d %d %s %s %d %lf %lf %lf\n",i,m2->atm[i].number,m2->atm[i].name,m2->atm[i].residue,m2->atm[i].resnum,m2->atm[i].x,m2->atm[i].y,m2->atm[i].z);
  //   }
  
  if (m1->atoms > 0){return(0);}
  else{
    //fprintf(stderr,"no identical atoms in files %s %s \n",m1->filename,m2->filename);
    return(1);
  }
}


int atom_selected(char* name,char* resname,molecule* m)
{
  int i;
  for(i=0;i<m->atoms;i++)
    {
      if(strcmp(name,m->atm[i].name)==0 && strcmp(m->atm[i].resname,resname)==0)
	return m->atm[i].selected;
    }
  return FALSE;
}




int atom_exist(char* name,char* resname,molecule* m)
{
  int i;
  for(i=0;i<m->atoms;i++)
    {
      if(strcmp(name,m->atm[i].name)==0 && strcmp(m->atm[i].resname,resname)==0)
	return TRUE;
    }
  return FALSE;
}

double Levitt_Gerstein(molecule *m1,molecule *m2,double d0)	
{
  int          i,j,last1,last2;
  int          numgap=0;
  double       sum=0.;
  //double       d0=25.; // 5**2
  //  double       d0=5.;  // 2.24**2
  double       M=20.;
  last1=m1->atm[0].resnum-1;
  last2=m2->atm[0].resnum-1;
  j=0;
  for (i=0;i<m1->atoms;i++)
    {
      if (m1->atm[i].selected){
	//printf("%d\n",m1->atm[i].resnum);
	last1++;
	last2++;
	if (m1->atm[i].resnum > last1 || m2->atm[i].resnum > last2 )
	  {
	    	  numgap++;
	  }
	sum+=1/(1+m1->atm[i].rms/d0);
	//	printf("Test>\t%d\t%f\t%f\n",numgap,sum,m1->atm[i].rms);
	last1=m1->atm[i].resnum;
	last2=m2->atm[i].resnum;
	j++;
	//	printf("Test> %d  \t%d\t%d\t%f\n",i,m1->atm[i].selected,sum);
      }
      //else
      //	{
      //  m1->atm[i].Sstr=0;
      //	}
    }
  //printf("TEST-LG>\t%d\t%d\t%f\t%f\n",j,numgap,sum,M*(sum-numgap/2));
  //numgap=0;
  return(M*(sum-numgap/2));
}

double LG_pvalue(int N,double MS)
{
  double min_SD=1.e-20;
  double ln,Mean_MS,SD_MS,Z,expZ;
  double ln60=log(120.);
  double b=2.0*ln60*18.411424-4.501719;
  double a=(ln60*ln60)*18.411424+ln60*(-4.501719)+2.637163-b*ln60;
  if (N<6) {
    expZ=1.;
    return(expZ);
  }
  ln=log((double)N);
  if (N<120){
    Mean_MS=ln*ln*18.411424+ln*(-4.501719)+ 2.637163;
    SD_MS=ln*21.351857 -37.521269;
  }else{
    Mean_MS=ln*b+a;
    SD_MS=ln60*21.351857 -37.521269;
  }
  //  if (SD_MS<min_SD) {SD_MS=min_SD;}
  Z=(MS-Mean_MS)/SD_MS;

  expZ=exp((double)-1.*Z);
  //printf("TEST-P> %d\t%f\t%f\t%e\t%e\t%e\t%e\n",N,MS,Mean_MS,SD_MS,Z,expZ,1-exp(-1*expZ));
  if (Z>20){
    return(expZ);
    //  }else if (Z<-100){
    //expZ=1.;
    //return(expZ);
  }else{
    return(1-exp(-expZ));
  }
}

double superimpose_molecules(m1,m2,s,error_cut)	/* Find RMS superimposition of m1 on m2 */
     molecule 	*m1,*m2;		/* Molecules to be superimposed */
     double 		s[3][3];	/* Final transformation matrix */
     double error_cut;
{
  int 		i,j,k;		/* Counter variables */
  int                 natoms=0;
  double		u[3][3];	/* direct product matrix */
  double 		t[3][3];	/* Temporary storage matrix */
  double		ma[3][3];	/* x axis rotation matrix */
  double		mb[3][3];	/* y axis rotation matrix */
  double		mg[3][3];	/* z axis rotation matrix */
  double 		*d1,*d2;	/* usefule pointers */
  double 		error;		/* Final superimposition error */
  double 		error2;		/* Final superimposition error */
  double		alpha=0.0; 	/* Angle of rotation around x axis */
  double		beta=0.0;	/* Angle of rotation around y axis */
  double		gamma=0.0;	/* Angle of rotation around z axis */
  double		dist,dist2,x,y,z;		/* Temporary coordinate variables */
  for (i=0;i<3;i++)  /* Initialize matrices */
    for (j=0;j<3;j++)
      s[i][j]=u[i][j]=0.0;
  s[0][0]=s[1][1]=s[2][2]=1.0;  /* Initialize S matrix to I */
  for (i=0;i<3;i++)  /* Initialize rotation matrices to I */
    for (j=0;j<3;j++)
      ma[i][j]=mb[i][j]=mg[i][j]=s[i][j];
  for (i=0;i<m1->atoms;i++)  /* Construct U matrix */
    {
      if (m1->atm[i].selected){
	d1= &(m1->atm[i].x);
	for (j=0;j<3;j++)
	  {
	    d2= &(m2->atm[i].x);
	    for (k=0;k<3;k++)
	      {
		u[j][k]+=(*d1)*(*d2);
		d2++;
	      }
	    d1++;
	  }
      }
    }
	  
  do
    {
      error=0.0;
      /* Calculate x axis rotation */
      alpha=atan((u[2][1]-u[1][2])/(u[1][1]+u[2][2]));
      /* Insure we are heading for a minimum, not a maximum */
      if (cos(alpha)*(u[1][1]+u[2][2])+sin(alpha)*(u[2][1]-u[1][2])<0.0)
	alpha+=PI;
      ma[1][1]=ma[2][2]=cos(alpha);
      ma[2][1]=sin(alpha); ma[1][2]= -ma[2][1];
      transpose_matrix(ma);
      multiply_matrix(u,ma,t);
      transpose_matrix(ma);
      copy_matrix(t,u);
      multiply_matrix(ma,s,t);
      copy_matrix(t,s);
      /* Calculate y axis rotation */
      beta=atan((u[0][2]-u[2][0])/(u[0][0]+u[2][2]));
      /* Insure we are heading for a minimum, not a maximum */
      if (cos(beta)*(u[0][0]+u[2][2])+sin(beta)*(u[0][2]-u[2][0])<0.0)
	beta+=PI;
      mb[0][0]=mb[2][2]=cos(beta);
      mb[0][2]=sin(beta); mb[2][0]= -mb[0][2];
      transpose_matrix(mb); 
      multiply_matrix(u,mb,t);
      transpose_matrix(mb);
      copy_matrix(t,u);
      multiply_matrix(mb,s,t);
      copy_matrix(t,s); 
      /* Calculate z axis rotation */
      gamma=atan((u[1][0]-u[0][1])/(u[0][0]+u[1][1]));
      /* Insure we are heading for a minimum, not a maximum */
      if (cos(gamma)*(u[0][0]+u[1][1])+sin(gamma)*(u[1][0]-u[0][1])<0.0)
	gamma+=PI;
      mg[0][0]=mg[1][1]=cos(gamma);
      mg[1][0]=sin(gamma); mg[0][1]= -mg[1][0];
      transpose_matrix(mg);
      multiply_matrix(u,mg,t);
      transpose_matrix(mg);
      copy_matrix(t,u);
      multiply_matrix(mg,s,t);
      copy_matrix(t,s);
      error=fabs(alpha)+fabs(beta)+fabs(gamma); 
      
    }
  while (error>error_cut); 	/* was 0.0001 before optimization Is error low enough to stop? */
  /* Now calculate final RMS superimposition */
  error=0.0;
  error2=0.0;
  natoms=0;
      //printf("testing7b %f %f\n",m2->atm[1].x,error);
      //printf ("s1:\t%f\t%f\t%f\n",s[0][0],s[0][1],s[0][2]);
      //printf ("s2:\t%f\t%f\t%f\n",s[1][0],s[1][1],s[1][2]);
      //printf ("s3:\t%f\t%f\t%f\n",s[2][0],s[2][1],s[2][2]);
  if (! (s[0][0]>0 || s[0][0]<0))  {
    //      printf ("bar\n");
    for (i=0;i<3;i++)  /* Initialize matrices */
      for (j=0;j<3;j++)
	s[i][j]=u[i][j]=0.0;
    s[0][0]=s[1][1]=s[2][2]=1.0;  /* Initialize S matrix to I */
  }
  for (i=0;i<m2->atoms;i++)
    {

      //dist=m2->atm[i].x*m2->atm[i].x+m2->atm[i].y*m2->atm[i].y+m2->atm[i].z*m2->atm[i].z;
      //printf("TEST1-m2\t %d\t%6.3f %6.3f %6.3f\t%e\n",i,m2->atm[i].x,m2->atm[i].y,m2->atm[i].z,dist);
      x=s[0][0]*m2->atm[i].x+s[0][1]*m2->atm[i].y+s[0][2]*m2->atm[i].z;
      y=s[1][0]*m2->atm[i].x+s[1][1]*m2->atm[i].y+s[1][2]*m2->atm[i].z;
      z=s[2][0]*m2->atm[i].x+s[2][1]*m2->atm[i].y+s[2][2]*m2->atm[i].z;
      m2->atm[i].x=x;
      m2->atm[i].y=y;
      m2->atm[i].z=z;
      //dist2=m2->atm[i].x*m2->atm[i].x+m2->atm[i].y*m2->atm[i].y+m2->atm[i].z*m2->atm[i].z;
      //printf("TEST2-m2\t %d\t%6.3f %6.3f %6.3f\t%e\t%e\n",i,m2->atm[i].x,m2->atm[i].y,m2->atm[i].z,dist2,dist2-dist);
      x=m1->atm[i].x-x;
      y=m1->atm[i].y-y;
      z=m1->atm[i].z-z;
      m2->atm[i].rms=m1->atm[i].rms=x*x+y*y+z*z;
      error+=m1->atm[i].rms;
      //printf("TEST:\t %d\t%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f \n",i,m1->atm[i].x,m2->atm[i].x,m1->atm[i].y,m2->atm[i].y,m1->atm[i].z,m2->atm[i].z);
      //printf("TEST1:\t %d\t%f %f %f\n",i,x*x,y*y,z*z);
      //printf("TEST2:\t %d %d\t%f\n",i,m1->atm[i].selected,m1->atm[i].rms);
      if (m1->atm[i].selected){
	error2+=m1->atm[i].rms;
	natoms++;
      }
      /*printf ("TEST1: %f %f %f\n",m1->atm[i].x,m1->atm[i].y,m1->atm[i].z);
	printf ("TEST2: %f %f %f\n",m2->atm[i].x,m2->atm[i].y,m2->atm[i].z);
	printf ("TEST3: %d %f %f\n",i,error,m1->atm[i].rms); */
    }
  //    printf("testing8 %f\n",m2->atm[1].x);
  error/=(double)(m1->atoms);
  error2/=(double)(natoms);
  if(natoms==0)
    error2=99999999;
  //printf("ERROR: %lf %lf %d\n",error2,sqrt(error2),natoms);
  return(sqrt(error2));
}

////////////////////////
//The following is added by Fang
double LG_pvalueF(N,score)
     int N;
     double score;
{
  double Z,s1,mean,SD,sc,expZ;
  s1=(double)N;
  mean=pow(s1,0.3264)/(0.0437*pow(s1,0.0003)+0.0790);
  SD=s1/(0.0417*s1+3.3700);
  Z=(score/10-mean)/SD;
  expZ=exp((double)-1.*Z);
  if (Z>20)s1=expZ;
  else s1=1-exp(-expZ);
  return s1;
}

void copy_coord(m1,m2)
     molecule 	*m1,*m2;
{
  int i;
  for (i=0;i<m1->atoms;i++){
    m2->atm[i].x=m1->atm[i].x;  
    m2->atm[i].y=m1->atm[i].y;  
    m2->atm[i].z=m1->atm[i].z;  
  }
}

void copy_molecule(m1,m2)
     molecule 	*m1,*m2;
{
  int i;
  m2->atoms=m1->atoms;
  m2->xcen=m1->xcen;
  m2->ycen=m1->ycen;
  m2->zcen=m1->zcen;
  strcpy(m2->filename,m1->filename);
  for (i=0;i<m1->atoms;i++){
    copy_res(m1,m2,i,i);
  }
}

void copy_res(m1,m2,i,j)
    molecule 	*m1,*m2;
    int i,j;
{
  m2->atm[j].x=m1->atm[i].x;  
  m2->atm[j].y=m1->atm[i].y;  
  m2->atm[j].z=m1->atm[i].z;  
  m2->atm[j].rms=m1->atm[i].rms;  
  //strcpy(m1->atm[i].name,m2->atm[j].name);
  //strcpy(m1->atm[i].residue,m2->atm[j].residue);
  m2->atm[j].number=m1->atm[i].number;  
  m2->atm[j].resnum=m1->atm[i].resnum;
  m2->atm[j].selected=m1->atm[i].selected;  
  m2->atm[j].bfactor=m1->atm[i].bfactor;
}


void copymolecule(molecule *m1,dyn_molecule *m2)
{
  /* Copy m2 to m1 */
  int i;
  for(i=0;i<m2->atoms;i++)
    {
      m1->atm[i].x=m2->atm[i].x;
      m1->atm[i].y=m2->atm[i].y;
      m1->atm[i].z=m2->atm[i].z;
      //m1->atm[i].rms=0; //m2->atm[i].rms;
      m1->atm[i].resnum=m2->atm[i].resnum;
      m1->atm[i].number=m2->atm[i].number;
      m1->atm[i].rescount=m2->atm[i].rescount;
      m1->atm[i].selected=TRUE;
      strcpy(m1->atm[i].name,m2->atm[i].name);
      strcpy(m1->atm[i].residue,m2->atm[i].residue);
      strcpy(m1->atm[i].resname,m2->atm[i].resname);
      m1->atm[i].selected=TRUE;
      if(strcmp("CA ",m2->atm[i].name) == 0)
	{
	  //printf("%s %s %s %s %d\n",m[0].filename,resname,residue,alt_loc,residues);
	  m1->CA_ref[m2->atm[i].rescount-1]=i;
	}
      //printf("%-30s %lf %lf %lf == %lf %lf %lf\n",m2->filename,m1->atm[i].x,m1->atm[i].y,m1->atm[i].z,m2->atm[i].x,m2->atm[i].y,m2->atm[i].z);
    }
  m1->xcen=m2->xcen;
  m1->ycen=m2->ycen;
  m1->zcen=m2->zcen;
  m1->atoms=m2->atoms;
  m1->residues=m2->residues;
  strcpy(m1->sequence,m2->sequence);
  //printf("%d %d",m1->atoms,m2->atoms);
  strcpy(m1->filename,m2->filename);
}

void reset_rms(molecule *m1)
{
  int i;
  for(i=0;i<m1->atoms;i++)
    {
      m1->atm[i].rms=0; 
    }
}




//End
///////////////////////////
