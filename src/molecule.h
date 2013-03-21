#ifndef molecule_h
#define molecule_h
#include "types/simple.h"
#include "nets.h"
#include <limits.h>
#define TRUE		1		/* Boolean definitions */
#define FALSE		0
#define	MAXATMS		10000		/* Maximum allowable atoms */
#define	MAXRES		2000	        /* Maximum allowable residues */
#define PI		3.14159265	/* Useful constant */
#undef max
#define max(a,b)    ((a) > (b) ? (a) : (b))
#undef min
#define min(a,b)    ((a) < (b) ? (a) : (b))

//#define SIZE (sizeof(a) / sizeof(a[0]))

typedef struct
 {
   double x,y,z;		/* Atomic coordinates */
   double rms;		/* RMS deviation */
   char residue[8];	/* PDB info for output */
   char name[8];
   int number;
   int resnum;
   char resname[8];
   int rescount;
   int selected;
   int deleted;
   double bfactor;
   char chain[2];
 } atm;

typedef struct {
  atm atm[MAXATMS];
  int CA_ref[MAXRES];
  int res_ref[MAXRES];
  double xcen,ycen,zcen;
  int	atoms;			/* Current # of atoms */
  int   residues;
  char  sequence[MAXRES];
  char  ss[MAXRES];
  double score;
  char method[200];
  int  rank;
  char	filename[PATH_MAX];		/* filename to read molecule from */
  
} molecule;


typedef struct {
  atm   *atm;
  int *CA_ref;//[MAXRES];
  int *res_ref;//[MAXRES];
  double xcen,ycen,zcen;
  int	atoms;			/* Current # of atoms */
  int   residues;
  char *sequence; // [MAXRES];
  char *ss; //[MAXRES];
  double score;
  char *method;
  int  rank;
  char	filename[1000];		/* filename to read molecule from */
  
} dyn_molecule;

//typedef struct {
//  struct
//  {
//    double x,y,z;		/* Atomic coordinates */
//    double rms;		/* RMS deviation */
//    char residue[8];	/* PDB info for output */
//    char name[8];
//    int number;
//    int resnum;
//    char resname[8];
//    int rescount;
//    int selected;
//  } atm[MAXRES];
//  int CA_ref[MAXRES];
//  int res_ref[MAXRES];
//  double xcen,ycen,zcen;
//  int	atoms;			/* Current # of atoms */
//  int   residues;
//  char  sequence[MAXRES];
//  char  ss[MAXRES];
//  double score;
//  char method[200];
//  int  rank;
//  char	filename[1000];		/* filename to read molecule from */
//  
//} molecule_ca;

typedef struct 
{
  double x,y,z;
}my_vector;


/* OBS!!!!!!!!!
   To relace *_ca and _backbone with one routine I
   changed input to read_molecules(molecule *m,char atomflag)
   atomflag == a -> read all atoms (except H)
   atomflag == c -> CA atoms
   atomflag == b -> backbone CA,C,N,O atoms
*/
int    read_molecules_dynamic(dyn_molecule* m,char atomflag, int *ignore_res);
int    free_dyn_molecule(dyn_molecule *m);
int    read_molecules(molecule *m,char atomflag);
int    read_molecules_ca(molecule *m);
int    read_molecules_backbone(molecule *m);
int    read_to_molecule(molecule *m,char** atom_vec,char** residue_vec,int* number_vec,rvec* coord_vec,size_t n);

void   strncpy_NULL(char *dest, char *src, size_t n);
int    get_atomtype3(char *name, char *res);
int    get_atomtype(char *name, char *res);
int    get_res6(char *res);
int    get_res6_no_pointer(char res);
int    get_res(char *res);
double distance(molecule *m,int atomno1, int atomno2);
void   print_type(int type_no, FILE *fp);
int    get_res(char *res);
void   print_res(int res,FILE *fp);
double crd(molecule *m,int atomno1, int atomno2);   /*closest residue distance */
double fatness(molecule *m);
double fatness2(molecule *m);
char   aa321(const char *name);
char*  assign_ss(molecule *m,float cutoff, float angle);
int    calc_index(int size,int row,int col);
int    hbond(molecule *m,int atomno1, int atomno2,float cutoff, float angle);
char*  read_psipred(char* filename);
char*  read_psipred2(char* filename,double* coil,double* helix,double* sheet);
void read_profile(char* filename,double (*prof)[22],char* profseq);
int aa(char aa);
char profile_index_to_aa(int i);
char*  read_stride(char* filename);
double* calculate_parameters(char* pdbfile,char* psipredfile);
double* ProQCA(char* pdbfile);
void ProQ(char** atom_vec,char** residue_vec,int* number_vec,rvec* coord_vec,size_t n,char* psipred,double* quality);
     
//double* ProQ(char** atom_vec,char** residue_vec,int* number_vec,double** coord_vec,size_t n,char* psipred);

#endif
