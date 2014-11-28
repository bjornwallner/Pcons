#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <sys/types.h>
#include <dirent.h>
#include <sys/times.h>
#include <time.h>
#include <linux/limits.h>
#include <math.h>
//#ifdef _OPENMP		      
#include <omp.h>
//#endif
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

void main(int argc,char *argv[])             /* Main routine */
{
  dyn_molecule  *dm;
  int files=10;
  int           ignore_res[2000]={0};
  lgscore *LG;
  int            L=4;
  int            step=2;
  double         minsim=121.0;
  double         d0=2.23606798;  //sqrt(5);
  double         factor=0.5;
  LG=malloc(sizeof(lgscore));
  dm=malloc(sizeof(dyn_molecule)*files); 
  strcpy(dm[0].filename,"large_set/3D-JIGSAW_V5-0_TS1");
  strcpy(dm[1].filename,"large_set/3D-JIGSAW_V5-0_TS2");
  for(int i=0;i<2;i++) {
    if(read_molecules_dynamic(&dm[i],'c',ignore_res)==0) {
      printf("Success reading %s\n",dm[i].filename);
    }
  }
  LGscore_res_pt(&dm[0],&dm[1],LG,d0,minsim,L,factor,step);
}
