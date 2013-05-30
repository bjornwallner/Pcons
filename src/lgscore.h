#ifndef lgscore_h
#define lgscore_h
#undef max
#define max(a,b)    ((a) > (b) ? (a) : (b))
#undef min
#define min(a,b)    ((a) < (b) ? (a) : (b))

#define	MAXRES		2000	        /* Maximum allowable residues */

typedef struct
{
  double S[MAXRES];
  double TM[MAXRES];
  double LGscore;
  double TMsum;
  double Ssum;
  char   residue[MAXRES];
  int    resnum[MAXRES];
  int length;
  int maxatoms;
  int    residues;
  double maxrms;
  double maxscore;
  double maxpvalue;
  

} lgscore;

typedef struct
{
  double x;
  double y;
  double z;
} vector;

double     LGscore(char *file1,char *file2, double minsim,int L,double factor);
double     superimpose(char *file1,char *file2,char* file3,double minsim,int L, double factor);
void       LGscore_res(char* file1,char* file2, lgscore *LG,double d0, double minsim,int L,double factor,int step);
void       rms(dyn_molecule *m1,dyn_molecule *m2,lgscore *LG, double d0);
void       LGscore_res_pt(dyn_molecule *m1,dyn_molecule *m2,lgscore *LG, double d0, double minsim,int L,double factor,int step);// ,int *ignore_res);
double	   LG_pvalue(int N,double MS);
double	   Levitt_Gerstein(molecule *m1,molecule *m2);
double	   superimpose_molecules();
double     LG_pvalueF(int N,double score);//added by Fang
//int        read_molecules_ca();	/* Reads in molecules to be superimposed */
//int        read_molecule(molecule *m);	/* Reads in molecules to be superimposed */
vector     center_molecule(molecule *m);	
int        copy_matrix(); /* copy matrix f into matrix t */
int        transpose_matrix(); /* Transpose a 3x3 matrix */
int        delete_atom(molecule *m1,int num);	
int        check_molecules(molecule *m1,molecule *m2);
int        check_molecules_mark_deleted(molecule *m1,molecule *m2);
int        check_molecules2(molecule *m1,molecule *m2);
int        atom_exist(char* name,char* resname,molecule* m);	
void       copymolecule(molecule *m1,dyn_molecule *m2);
void       copy_molecule();
void       copy_res();
void       copy_coord();
void reset_rms(molecule *m1);
int atom_selected(char* name,char* resname,molecule* m);
#endif
